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

  double quat[4];
  colvar_value = rmsd_grad_gpu(quat);
  update_hills();

  int grid_index = (int) floor((colvar_value-lower_boundary)/width);

  if( use_grids && grid_index<hills_grid_size ) {
    colvar_energy = hills_grid[grid_index][0];
    colvar_force = hills_grid[grid_index][1];
    //std::cerr << fmt::format("   *** hills_grid[{}] {} {}\n", grid_index,hills_grid[grid_index][0],hills_grid[grid_index][1] );
  } else {
    calc_energy_and_force(colvar_value,colvar_energy,colvar_force);
  }

  auto d_f = atomKK->k_f.template view<DeviceType>();
  double drmsddx = (colvar_value>0.0) ? 1.0/(colvar_value*group_count) : 0.0;
  double inverse_quat[4];
  qconjugate(quat,inverse_quat);

  Kokkos::parallel_for(group_count, LAMMPS_LAMBDA(const int& n) {

    const int i = atomKK->map(d_group_taglist[n]);

    double grad[3], grad_rot[3];
    grad[0] = drmsddx*(d_x_group(n,0) - d_ref_positions(n,0));
    grad[1] = drmsddx*(d_x_group(n,1) - d_ref_positions(n,1));
    grad[2] = drmsddx*(d_x_group(n,2) - d_ref_positions(n,2));

    // quaternion rotation of vector: c = a*b*conj(a)
    quatrotvec(inverse_quat, grad, grad_rot);

    d_f(i,0) += colvar_force * grad_rot[0];
    d_f(i,1) += colvar_force * grad_rot[1];
    d_f(i,2) += colvar_force * grad_rot[2];

    //std::cerr << fmt::format(" *** x[{}] {:.6} {:.6} {:.6} grad {:.6} {:.6} {:.6} f[{}] {:.6} {:.6} {:.6}\n", i,d_x_group(n,0),d_x_group(n,1),d_x_group(n,2),grad[0],grad[1],grad[2], i,d_f(i,0),d_f(i,1),d_f(i,2) );

    //std::cerr << fmt::format(" *** f[{}] {:.6} {:.6} {:.6}\n", i,d_f(i,0),d_f(i,1),d_f(i,2) );

  });

  copymode = 0;

  atomKK->modified(execution_space, F_MASK);

  //std::cerr << fmt::format(" *** FixMetadynamicsKK colvar_value {:.6} colvar_energy {:.6} colvar_force {:.6} quat {:.6} {:.6} {:.6} {:.6}\n", colvar_value,colvar_energy,colvar_force, quat[0],quat[1],quat[2],quat[3]);


}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
double FixMetadynamicsKokkos<DeviceType>::rmsd_grad_gpu(double quat[4])
{
  auto d_x = atomKK->k_x.template view<DeviceType>();
  auto d_image = atomKK->k_image.template view<DeviceType>();

  Kokkos::parallel_for(group_count, LAMMPS_LAMBDA(const int& n) {
      const int i = atomKK->map(d_group_taglist[n]);
      domain->unmap(&d_x(i,0),d_image[i],&d_x_group(n,0));
      //std::cerr << fmt::format("d_group_taglist[{}] {}\n", n,d_group_taglist[n]);
    });

  return gpu_q_j(d_x_group,d_ref_positions,quat);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
double FixMetadynamicsKokkos<DeviceType>::gpu_q_j(typename AT::t_x_array p, typename AT::t_x_array q, double quat[4])
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

  quat[0] = Evc[0][0];
  quat[1] = Evc[0][1];
  quat[2] = Evc[0][2];
  quat[3] = Evc[0][3];

  //double mat[3][3];
  //quat_to_mat(quat,mat);
  //write3(mat);

  double d_x_group_center[3] = { px/N, py/N, pz/N };
  double d_rp_center[3] = { qx/N, qy/N, qz/N };

  double e = px_px+py_py+pz_pz + qx_qx+qy_qy+qz_qz
    - ( px*px+py*py+pz*pz + qx*qx+qy*qy+qz*qz)/N;

  double d_rp_center_rotated[3];
  quatrotvec(quat, d_rp_center, d_rp_center_rotated);
  double rmsd = Kokkos::sqrt(Kokkos::fdim(e, 2.0 * Evl[0])/N);

  Kokkos::parallel_for(N, LAMMPS_LAMBDA(const int& n) {

    p(n,0) -= d_x_group_center[0];
    p(n,1) -= d_x_group_center[1];
    p(n,2) -= d_x_group_center[2];

    //std::cerr << fmt::format(" *** AFTER_CG d_x_group[{}] {:.6} {:.6} {:.6} d_x_group_center {:.6} {:.6} {:.6} \n", n,p(n,0),p(n,1),p(n,2), d_x_group_center[0],d_x_group_center[1],d_x_group_center[2] );

    double tmp[3];
    quatrotvec(quat, &p(n,0), tmp);

    //std::cerr << fmt::format(" *** AFTER_ROT tmp {:.6} {:.6} {:.6} quat {:.6} {:.6} {:.6} {:.6}\n", tmp[0],tmp[1],tmp[2], quat[0],quat[1],quat[2],quat[3]);

    p(n,0) = tmp[0] + d_rp_center[0];
    p(n,1) = tmp[1] + d_rp_center[1];
    p(n,2) = tmp[2] + d_rp_center[2];

    //std::cerr << fmt::format(" *** AFTER_TRANSLATION d_x_group[{}] {:.6} {:.6} {:.6}\n", n,p(n,0),p(n,1),p(n,2) );

  });
  return rmsd;

}



















































namespace LAMMPS_NS {
template class FixMetadynamicsKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixMetadynamicsKokkos<LMPHostType>;
#endif
}

