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
#include "atom_masks.h"
#include "domain.h"
#include "error.h"
//#include "kokkos_base.h"
#include "kokkos_few.h"
#include "math_eigen_impl.h"
#include "math_extra.h"
//#include "modify.h"
#include "math_special_kokkos.h"
#include "memory_kokkos.h"
#include "update.h"

#include <iostream>


using namespace LAMMPS_NS;
using namespace MathExtra;
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

  memoryKK->create_kokkos(k_group_taglist,group_taglist,group_count,"metadynamics:group_taglist");
  memoryKK->create_kokkos(k_ref_positions,group_count,"metadynamics:ref_positions");
  memoryKK->create_kokkos(d_x_group,group_count,"metadynamics:group_taglist");
  memoryKK->create_kokkos(d_x_group_shifted,group_count,"metadynamics:d_x_group_shifted");
  memoryKK->create_kokkos(d_ref_positions_shifted,group_count,"metadynamics:d_ref_positions_shifted");

  for( int i=0 ; i<group_count ; i++ ) {
    k_ref_positions.h_view(i,0) = ref_positions[i][0];
    k_ref_positions.h_view(i,1) = ref_positions[i][1];
    k_ref_positions.h_view(i,2) = ref_positions[i][2];
  }

  d_group_taglist = k_group_taglist.template view<DeviceType>();
  d_ref_positions = k_ref_positions.template view<DeviceType>();

  k_ref_positions.template modify<LMPHostType>();
  k_ref_positions.template sync<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixMetadynamicsKokkos<DeviceType>::~FixMetadynamicsKokkos()
{
  if (copymode) return;
  memoryKK->destroy_kokkos(k_group_taglist);
  memoryKK->destroy_kokkos(d_x_group);
  memoryKK->destroy_kokkos(d_x_group_shifted);
  memoryKK->destroy_kokkos(d_ref_positions_shifted);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixMetadynamicsKokkos<DeviceType>::init()
{
  k_group_taglist.template modify<LMPHostType>();
  FixMetadynamics::init();
  k_group_taglist.template sync<DeviceType>();

  if (utils::strmatch(update->integrate_style,"^respa"))
    error->all(FLERR,"Cannot (yet) use respa with Kokkos");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixMetadynamicsKokkos<DeviceType>::post_force(int /*vflag*/)
{
  atomKK->sync(execution_space, X_MASK | F_MASK );
  copymode = 1;

  double inverse_quat[4];
  colvar_value = rmsd_grad_gpu(inverse_quat);
  update_hills();

  int grid_index = (int) floor((colvar_value-lower_boundary)/width);

  if( use_grids && grid_index<hills_grid_size ) {
    colvar_energy = hills_grid[grid_index][0];
    colvar_force = hills_grid[grid_index][1];
    std::cerr << fmt::format("   *** hills_grid[{}] {} {}\n", grid_index,hills_grid[grid_index][0],hills_grid[grid_index][1] );
  } else {
    calc_energy_and_force(colvar_value,colvar_energy,colvar_force);
  }

  auto d_f = atomKK->k_f.template view<DeviceType>();
  double drmsddx = (colvar_value>0.0) ? 1.0/(colvar_value*group_count) : 0.0;
  double quat[4];
  qconjugate(inverse_quat,quat);

  Kokkos::parallel_for(group_count, LAMMPS_LAMBDA(const int& n) {
    const int i = atomKK->map(d_group_taglist[n]);


    //double grad[3];
    /*, d_rp_n_rotated[3];
    quatrotvec(quat, &d_ref_positions(n,0), d_rp_n_rotated);
    grad[0] = drmsddx*(d_x_group(n,0) - d_rp_n_rotated[0]);
    grad[1] = drmsddx*(d_x_group(n,1) - d_rp_n_rotated[1]);
    grad[2] = drmsddx*(d_x_group(n,2) - d_rp_n_rotated[2]);
    */

    d_f(i,0) += colvar_force * d_x_group(n,0);
    d_f(i,1) += colvar_force * d_x_group(n,1);
    d_f(i,2) += colvar_force * d_x_group(n,2);

    //std::cerr << fmt::format(" *** x[{}] {:.6} {:.6} {:.6} grad {:.6} {:.6} {:.6} f[{}] {:.6} {:.6} {:.6}\n", i,d_x_group(n,0),d_x_group(n,1),d_x_group(n,2),grad[0],grad[1],grad[2], i,d_f(i,0),d_f(i,1),d_f(i,2) );

    std::cerr << fmt::format(" *** grad[{}] {:.6} {:.6} {:.6} f[{}] {:.6} {:.6} {:.6}\n", n,d_x_group(n,0),d_x_group(n,1),d_x_group(n,2), i,d_f(i,0),d_f(i,1),d_f(i,2) );

  });

  copymode = 0;

  atomKK->modified(execution_space, F_MASK);

  std::cerr << fmt::format(" *** FixMetadynamicsKK colvar_value {:.6} colvar_energy {:.6} colvar_force {:.6} quat {:.6} {:.6} {:.6} {:.6}\n", colvar_value,colvar_energy,colvar_force, quat[0],quat[1],quat[2],quat[3]);


}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
double FixMetadynamicsKokkos<DeviceType>::rmsd_grad_gpu(double inverse_quat[4])
{
  auto d_x = atomKK->k_x.template view<DeviceType>();
  auto d_image = atomKK->k_image.template view<DeviceType>();

  Kokkos::parallel_for(group_count, LAMMPS_LAMBDA(const int& n) {
      const int i = atomKK->map(d_group_taglist[n]);
      domain->unmap(&d_x(i,0),d_image[i],&d_x_group(n,0));
      //std::cerr << fmt::format("d_group_taglist[{}] {}\n", n,d_group_taglist[n]);
    });

  return gpu_q_j(d_x_group,d_ref_positions,inverse_quat);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
double FixMetadynamicsKokkos<DeviceType>::gpu_q_j(typename AT::t_x_array p, typename AT::t_x_array q, double inverse_quat[4])
{

  int N = p.extent_int(0);

  double px_qx, px_qy, px_qz;
  double py_qx, py_qy, py_qz;
  double pz_qx, pz_qy, pz_qz;
  double px_px, py_py, pz_pz;
  double qx_qx, qy_qy, qz_qz;
  double px, py, pz, qx, qy, qz;

  Kokkos::parallel_reduce(N, LAMMPS_LAMBDA(const int& n,
    double& l_px_qx, double& l_px_qy, double& l_px_qz,
    double& l_py_qx, double& l_py_qy, double& l_py_qz,
    double& l_pz_qx, double& l_pz_qy, double& l_pz_qz,
    double& l_px_px, double& l_py_py, double& l_pz_pz,
    double& l_qx_qx, double& l_qy_qy, double& l_qz_qz,
    double& l_px, double& l_py, double& l_pz,
    double& l_qx, double& l_qy, double& l_qz ) {

      l_px_qx+=p(n,0)*q(n,0); l_px_qy+=p(n,0)*q(n,1); l_px_qz+=p(n,0)*q(n,2);
      l_py_qx+=p(n,1)*q(n,0); l_py_qy+=p(n,1)*q(n,1); l_py_qz+=p(n,1)*q(n,2);
      l_pz_qx+=p(n,2)*q(n,0); l_pz_qy+=p(n,2)*q(n,1); l_pz_qz+=p(n,2)*q(n,2);

      l_px_px+=p(n,0)*p(n,0); l_py_py+=p(n,1)*p(n,1); l_pz_pz+=p(n,2)*p(n,2);
      l_qx_qx+=q(n,0)*q(n,0); l_qy_qy+=q(n,1)*q(n,1); l_qz_qz+=q(n,2)*q(n,2);

      l_px+=p(n,0); l_py+=p(n,1); l_pz+=p(n,2);
      l_qx+=q(n,0); l_qy+=q(n,1); l_qz+=q(n,2);

    },
    px_qx, px_qy, px_qz,
    py_qx, py_qy, py_qz,
    pz_qx, pz_qy, pz_qz,
    px_px, py_py, pz_pz,
    qx_qx, qy_qy, qz_qz,
    px, py, pz, qx, qy, qz);

  double r11 = px_qx - px*qx/N;
  double r12 = px_qy - px*qy/N;
  double r13 = px_qz - px*qz/N;

  double r21 = py_qx - py*qx/N;
  double r22 = py_qy - py*qy/N;
  double r23 = py_qz - py*qz/N;

  double r31 = pz_qx - pz*qx/N;
  double r32 = pz_qy - pz*qy/N;
  double r33 = pz_qz - pz*qz/N;

  double P[4][4], Evl[4], Evc[4][4];
  P[0][0]=r11+r22+r33; P[0][1]=r23-r32; P[0][2]=r31-r13; P[0][3]=r12-r21;
  P[1][0]=r23-r32; P[1][1]=r11-r22-r33; P[1][2]=r12+r21; P[1][3]=r13+r31;
  P[2][0]=r31-r13; P[2][1]=r12+r21; P[2][2]=-r11+r22-r33; P[2][3]=r23+r32;
  P[3][0]=r12-r21; P[3][1]=r13+r31; P[3][2]=r23+r32; P[3][3]=-r11-r22+r33;

  MathEigen::Jacobi<double,double*,double(*)[4],double const(*)[4]> ecalc(4);
  int ierror = ecalc.Diagonalize(P,Evl,Evc);

  if(ierror)
    error->all(FLERR, "fix metadynamics: Too many iterations in jacobi diagonalization.\n"
      "This is usually the result of an ill-defined set of atoms for "
      "rotational alignment (RMSD, rotateReference, etc).\n");

  inverse_quat[0] = Evc[0][0];
  inverse_quat[1] = -Evc[0][1];
  inverse_quat[2] = -Evc[0][2];
  inverse_quat[3] = -Evc[0][3];

  double mat[3][3];
  quat_to_mat(inverse_quat,mat);
  write3(mat);

  double d_x_group_center[3] = { px/N, py/N, pz/N };
  double d_rp_center[3] = { qx/N, qy/N, qz/N };

  double e = px_px+py_py+pz_pz + qx_qx+qy_qy+qz_qz
    - ( px*px+py*py+pz*pz + qx*qx+qy*qy+qz*qz)/N;

  double d_rp_center_rotated[3];
  quatrotvec(inverse_quat, d_rp_center, d_rp_center_rotated);
  double rmsd = Kokkos::sqrt(Kokkos::fdim(e, 2.0 * Evl[0])/N);


  Kokkos::parallel_for(N, LAMMPS_LAMBDA(const int& n) {

    d_x_group(n,0) -= d_x_group_center[0];
    d_x_group(n,1) -= d_x_group_center[1];
    d_x_group(n,2) -= d_x_group_center[2];

    std::cerr << fmt::format(" *** AFTER_CG d_x_group[{}] {:.6} {:.6} {:.6}\n", n,d_x_group(n,0),d_x_group(n,1),d_x_group(n,2) );

    double tmp[3];
    quatrotvec(inverse_quat, &d_x_group(n,0), tmp);

    std::cerr << fmt::format(" *** AFTER_ROT inverse_quat {:.6} {:.6} {:.6} {:.6}  tmp {:.6} {:.6} {:.6}\n", inverse_quat[0],inverse_quat[1],inverse_quat[2],inverse_quat[3], tmp[0],tmp[1],tmp[2]);

    d_x_group(n,0) = tmp[0] + d_x_group_center[0];
    d_x_group(n,1) = tmp[1] + d_x_group_center[1];
    d_x_group(n,2) = tmp[2] + d_x_group_center[2];

    std::cerr << fmt::format(" *** AFTER_TRANSLATION d_x_group[{}] {:.6} {:.6} {:.6}\n", n,x_group[n][0],x_group[n][1],x_group[n][2]);

  });
  return rmsd;

}























































/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
double FixMetadynamicsKokkos<DeviceType>::rmsd(double inverse_quat[4])
{

  double aCenter_f[3] = {0.0};
  double aCenter_m[3] = {0.0};

  // Find the center-of-geometry of each object:
  auto d_x = atomKK->k_x.template view<DeviceType>();
  auto d_image = atomKK->k_image.template view<DeviceType>();

  Kokkos::parallel_reduce(group_count,
    LAMMPS_LAMBDA(const int& n, double& fx, double& fy, double& fz, double& mx, double& my, double& mz ) {
      const int i = atomKK->map(d_group_taglist[n]);
      domain->unmap(&d_x(i,0),d_image[i],&d_x_group(n,0));
      fx += d_x_group(n,0);
      fy += d_x_group(n,1);
      fz += d_x_group(n,2);
      mx += d_ref_positions(n,0);
      my += d_ref_positions(n,1);
      mz += d_ref_positions(n,2);
      std::cerr << fmt::format("d_group_taglist[{}] {}\n", n,d_group_taglist[n]);
    },aCenter_f[0],aCenter_f[1],aCenter_f[2],aCenter_m[0],aCenter_m[1],aCenter_m[2]);

  std::cerr << fmt::format(" *** aCenter_f {:.6} {:.6} {:.6} aCenter_m {:.6} {:.6} {:.6}\n", aCenter_f[0],aCenter_f[1],aCenter_f[2], aCenter_m[0],aCenter_m[1],aCenter_m[2]);


    //std::cerr << fmt::format(" *** BEFORE_CG x_group[{}] {:.6} {:.6} {:.6}\n", n,x_group[n][0],x_group[n][1],x_group[n][2]);

  for (int d = 0; d < 3; d++) {
    aCenter_f[d] /= group_count;
    aCenter_m[d] /= group_count;
  }

  // Calculate the "M" array from the Diamond paper (equation 16)
  double M[3][3];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {

    M[i][j] = 0.0;

    Kokkos::parallel_reduce(group_count,
      LAMMPS_LAMBDA(const int& n, double& mij ) {

        // x_group_shifted[n][d] = x_group[n][d] - aCenter_f[d];
        // ref_positions_shifted[n][d] = ref_positions[n][d] - aCenter_m[d];
        // M[i][j] += ref_positions_shifted[n][i] * x_group_shifted[n][j]

        d_ref_positions_shifted(n,i) = d_ref_positions(n,i) - aCenter_m[i];
        d_x_group_shifted(n,j) = d_x_group(n,j) - aCenter_f[j];
        mij += d_ref_positions_shifted(n,i)*d_x_group_shifted(n,j);

    }, M[i][j] );
  }

  write3(M);

  // Calculate Q (equation 17)
  double traceM = 0.0;
  for (int i = 0; i < 3; i++) traceM += M[i][i];
  double Q[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      Q[i][j] = M[i][j] + M[j][i];
      if (i == j) Q[i][j] -= 2.0 * traceM;
    }
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
  double p[4] = {0.0, 0.0, 0.0, 1.0};

  double Evl[4];                 // Store the eigenvalues of P here.
  double *Evc[4];                // Store the eigevectors here. This version has ** format.
  double _Evc[4 * 4];            // Contiguous 1D array for storing contents of "Evc" array
  for (int i = 0; i < 4; i++)    // We must make sure that
    Evc[i] = &(_Evc[4 * i]);     // Evc[i] points to the correct location in memory

  //for (unsigned i = 0; i < 4; i++) {
  //  printf("[ ");
  //  for (unsigned j = 0; j < 4; j++) printf("%g ",P[i][j]);
  //  printf("]\n");
  //}

  int ierror = MathEigen::Jacobi<double, double *, double **>(4).Diagonalize(P, Evl, Evc);

  if(ierror)
    error->all(FLERR, "fix metadynamics: Too many iterations in jacobi diagonalization.\n"
      "This is usually the result of an ill-defined set of atoms for "
      "rotational alignment (RMSD, rotateReference, etc).\n");

  for (int i = 0; i < 4; i++)
    p[i] = Evc[0][i];    //copy eigenvector corresponding to this eigenvalue to p

  // Now normalize p
  double pnorm = 0.0;
  for (int i = 0; i < 4; i++) pnorm += p[i] * p[i];
  pnorm = Kokkos::sqrt(pnorm);
  for (int i = 0; i < 4; i++) p[i] /= pnorm;


  // Note: The "p" variable is not a quaternion in the
  //       conventional sense because its elements
  //       are in the wrong order.  I correct for that here.
  //       "q" is the quaternion correspond to rotation R.
  //       [Andrew Jewett (Scripps Research)]
  //q[0] = p[3];
  //q[1] = p[0];
  //q[2] = p[1];
  //q[3] = p[2];

  // BUT... we actually need the INVERSE ROTATION
  // [alphataubio (2024/08)]

  inverse_quat[0] = p[3];
  inverse_quat[1] = -p[0];
  inverse_quat[2] = -p[1];
  inverse_quat[3] = -p[2];

  // Finally compute the RMSD between the two coordinate sets:
  // First compute E0 from equation 24 of the paper
  double E0 = 0.0;

  Kokkos::parallel_reduce(group_count,
    LAMMPS_LAMBDA(const int& n, double& e0 ) {
      double tmp[3];
      // quaternion rotation of vector: c = a*b*conj(a)
      quatrotvec(inverse_quat, &d_x_group_shifted(n,0), tmp);
      //std::cerr << fmt::format(" *** AFTER_ROT x_group[{}] {:.6} {:.6} {:.6}\n", n,tmp[0],tmp[1],tmp[2]);
      for (int d = 0; d < 3; d++) {
        e0 += square(d_x_group_shifted(n,d) - d_ref_positions_shifted(n,d));
        d_x_group(n,d) = tmp[d]+aCenter_m[d];
      }
      std::cerr << fmt::format(" *** AFTER_TRANSLATION x_group[{}] {:.6} {:.6} {:.6}\n", n,d_x_group(n,0),d_x_group(n,1),d_x_group(n,2));
    },E0);

  return Kokkos::sqrt(Kokkos::fdim(E0, 2.0 * Evl[0])/group_count); // Evl[0] = the maximum eigenvalue of P

}

namespace LAMMPS_NS {
template class FixMetadynamicsKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixMetadynamicsKokkos<LMPHostType>;
#endif
}

