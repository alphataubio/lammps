// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_metadynamics_kokkos.h"

#include "atom_kokkos.h"
#include "update.h"
#include "math_eigen_impl.h"
#include "modify.h"
#include "region.h"
#include "input.h"
#include "variable.h"
#include "math_special_kokkos.h"
#include "memory_kokkos.h"
#include "error.h"
#include "atom_masks.h"
#include "kokkos_base.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathSpecialKokkos;


/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixMetadynamicsKokkos<DeviceType>::FixMetadynamicsKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixMetadynamics(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK;
  datamask_modify = F_MASK;


}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixMetadynamicsKokkos<DeviceType>::~FixMetadynamicsKokkos()
{
  if (copymode) return;

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixMetadynamicsKokkos<DeviceType>::init()
{
  FixMetadynamics::init();

  if (utils::strmatch(update->integrate_style,"^respa"))
    error->all(FLERR,"Cannot (yet) use respa with Kokkos");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixMetadynamicsKokkos<DeviceType>::post_force(int /*vflag*/)
{
  atomKK->sync(execution_space, F_MASK | MASK_MASK);

  d_f = atomKK->k_f.template view<DeviceType>();
  d_mask = atomKK->k_mask.view<DeviceType>();

  int nlocal = atom->nlocal;

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixMetadynamics>(0,nlocal),*this);
  copymode = 0;


  atomKK->modified(execution_space, F_MASK);
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixMetadynamicsKokkos<DeviceType>::operator()(TagFixMetadynamics, const int &i) const {

}

/* ---------------------------------------------------------------------- */

using std::fdim;
using std::sqrt;

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
double FixMetadynamicsKokkos<DeviceType>::rmsd()
{

  typename AT::t_x_array aaXf = atomKK->k_f.template view<DeviceType>();
  typename AT::t_x_array aaXm;
  int N = aaXf.extent(0);

  // Find the center of mass of each object:
  double aCenter_f[3] = {0.0, 0.0, 0.0};
  double aCenter_m[3] = {0.0, 0.0, 0.0};
  for (int n = 0; n < N; n++)
    for (int d = 0; d < 3; d++) {
      aCenter_f[d] += aaXf(n,d);
      aCenter_m[d] += aaXm(n,d);
    }

  for (int d = 0; d < 3; d++) {
    aCenter_f[d] /= N;
    aCenter_m[d] /= N;
  }

  //Subtract the centers-of-mass from the original coordinates for each object
  for (int n = 0; n < N; n++)
    for (int d = 0; d < 3; d++) {
      // shift the coordinates so that the new center of mass is at the origin
      aaXf_shifted[n][d] = aaXf(n,d) - aCenter_f[d];
      aaXm_shifted[n][d] = aaXm(n,d) - aCenter_m[d];
    }

  // Calculate the "M" array from the Diamond paper (equation 16)
  double M[3][3];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) M[i][j] = 0.0;

  for (size_t n = 0; n < N; n++)
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) { M[i][j] += aaXm_shifted[n][i] * aaXf_shifted[n][j]; }
    }

  // Calculate Q (equation 17)
  double traceM = 0.0;
  for (int i = 0; i < 3; i++) traceM += M[i][i];
  double Q[3][3];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      Q[i][j] = M[i][j] + M[j][i];
      if (i == j) Q[i][j] -= 2.0 * traceM;
    }

  // Calculate V (equation 18)
  double V[3];
  V[0] = M[1][2] - M[2][1];
  V[1] = M[2][0] - M[0][2];
  V[2] = M[0][1] - M[1][0];

  // Calculate "P" (equation 22)
  // First we must allocate space for the P matrix.  It's not safe to declare:
  // double P[4][4];
  // ...because most matrix solvers expect arrays in pointer-to-pointer format.
  // (a different format).  Below I create a fixed size matrix P in this format.
  double _PF[4 * 4];             // Contiguous 1D array for storing contents of the 2D P array
  double *P[4];                  // This version of P has has ** (pointer-to-pointer) format.
  for (int i = 0; i < 4; i++)    // We must make sure that
    P[i] = &(_PF[4 * i]);        // P[i] points to the appropriate location in memory

  // Now fill the P array
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) P[i][j] = Q[i][j];
  P[0][3] = V[0];
  P[3][0] = V[0];
  P[1][3] = V[1];
  P[3][1] = V[1];
  P[2][3] = V[2];
  P[3][2] = V[2];
  P[3][3] = 0.0;

  // The vector "p" contains the optimal rotation (backwards quaternion format)
  double p[4] = {0.0, 0.0, 0.0, 1.0};    // default value
  double pPp = 0.0;                      // = p^T * P * p  (zero by default)

  double Evl[4];                 // Store the eigenvalues of P here.
  double *Evc[4];                // Store the eigevectors here. This version has ** format.
  double _Evc[4 * 4];            // Contiguous 1D array for storing contents of "Evc" array
  for (int i = 0; i < 4; i++)    // We must make sure that
    Evc[i] = &(_Evc[4 * i]);     // Evc[i] points to the correct location in memory

  MathEigen::Jacobi<double, double *, double **>(4).Diagonalize(P, Evl, Evc);

  // Note: The eigenvalues are sorted in decreasing order by default.
  pPp = Evl[0];    // = the maximum eigenvalue of P
  for (int i = 0; i < 4; i++)
    p[i] = Evc[0][i];    //copy eigenvector corresponding to this eigenvalue to p

  // Now normalize p
  double pnorm = 0.0;
  for (int i = 0; i < 4; i++) pnorm += p[i] * p[i];
  pnorm = sqrt(pnorm);
  for (int i = 0; i < 4; i++) p[i] /= pnorm;

  // Finally compute the RMSD between the two coordinate sets:
  // First compute E0 from equation 24 of the paper
  double E0 = 0.0;
  for (size_t n = 0; n < N; n++)
    for (int d = 0; d < 3; d++)
      E0 += (square(aaXf_shifted[n][d] - aaXm_shifted[n][d]));

  return sqrt(fdim(E0, 2.0 * pPp)/N);

}



namespace LAMMPS_NS {
template class FixMetadynamicsKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixMetadynamicsKokkos<LMPHostType>;
#endif
}

