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

#include "fix_metadynamics.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "input.h"
#include "math_eigen_impl.h"    //functions to calculate eigenvalues and eigenvectors
#include "math_special.h"
#include "memory.h"
#include "modify.h"
#include "region.h"
#include "respa.h"
#include "update.h"
#include "variable.h"
#include "superpose3d.h"


#include <cmath>
#include <cstring>

#include <iostream>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathSpecial;


/* ---------------------------------------------------------------------- */

FixMetadynamics::FixMetadynamics(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), xstr(nullptr), ystr(nullptr), zstr(nullptr), idregion(nullptr),
    region(nullptr)
{
  //if (narg < 6) utils::missing_cmd_args(FLERR, "fix setforce", error);

  lowerBoundary = 0.0;
  upperBoundary = 30.0;
  width = 0.25;

  hillWeight = 0.1;
  hillWidth = 2.0;
  newHillFrequency = 10;
  outputFreq = 100;

  lowerWalls = 0.5;
  upperWalls = 29.5;
  forceConstant = 10.0;

  N=3;

  MathEigen::Alloc2D(N, 3, &aaXf_shifted);
  MathEigen::Alloc2D(N, 3, &aaXm_shifted);
}

/* ---------------------------------------------------------------------- */

FixMetadynamics::~FixMetadynamics()
{
  if (copymode) return;


}

/* ---------------------------------------------------------------------- */

int FixMetadynamics::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMetadynamics::init()
{

}

/* ---------------------------------------------------------------------- */

void FixMetadynamics::setup(int vflag)
{
  if (utils::strmatch(update->integrate_style, "^verlet"))
    post_force(vflag);
  else
    error->all(FLERR, "fix metadynamics does not support RESPA.");
}


/* ---------------------------------------------------------------------- */

void FixMetadynamics::post_force(int /*vflag*/)
{
  double **x = atom->x;
  //double **f = atom->f;
  //int *mask = atom->mask;
  //int nlocal = atom->nlocal;

  //double b[3][3] = {0.0};

  double **b;
  memory->create(b,3,3,"metadynamics:b");

  Superpose3D<double, double **> superposer(3);
  double tmp = superposer.Superpose(x, b);

  domain->remap(x[0]);
  domain->remap(x[1]);
  domain->remap(x[2]);
  //domain->minimum_image(x[0][0],x[0][1],x[0][2]);
  //domain->minimum_image(x[1][0],x[1][1],x[1][2]);
  //domain->minimum_image(x[2][0],x[2][1],x[2][2]);

  domain->pbc();

  std::cerr << fmt::format("*** x[0] = {} {} {} b[0] = {} {} {} rmsd = {} Superpose3D rmsd {}\n",
    x[0][0],x[0][1],x[0][2],b[0][0],b[0][1],b[0][2],rmsd(x,b,3),tmp );

}

/* ---------------------------------------------------------------------- */

// FIXME: should it be static inline, or inline. do i care...

using std::fdim;
using std::sqrt;

double FixMetadynamics::rmsd( double **aaXf, double **aaXm, int N )
{
  // Find the center of mass of each object:
  double aCenter_f[3] = {0.0, 0.0, 0.0};
  double aCenter_m[3] = {0.0, 0.0, 0.0};
  for (int n = 0; n < N; n++)
    for (int d = 0; d < 3; d++) {
      aCenter_f[d] += aaXf[n][d];
      aCenter_m[d] += aaXm[n][d];
    }

  for (int d = 0; d < 3; d++) {
    aCenter_f[d] /= N;
    aCenter_m[d] /= N;
  }

  //Subtract the centers-of-mass from the original coordinates for each object
  for (int n = 0; n < N; n++)
    for (int d = 0; d < 3; d++) {
      // shift the coordinates so that the new center of mass is at the origin
      aaXf_shifted[n][d] = aaXf[n][d] - aCenter_f[d];
      aaXm_shifted[n][d] = aaXm[n][d] - aCenter_m[d];
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

