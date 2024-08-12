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
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "group.h"
#include "input.h"
#include "math_eigen_impl.h"
#include "math_extra.h"
#include "math_special.h"
#include "memory.h"
#include "modify.h"
#include "region.h"
#include "respa.h"
#include "update.h"
#include "variable.h"
#include "superpose3d.h"


#include <algorithm>
#include <cmath>
#include <cstring>

#include <iostream>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathExtra;
using namespace MathSpecial;

static constexpr int HILLS_DELTA = 16;

/* ---------------------------------------------------------------------- */

FixMetadynamics::FixMetadynamics(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg)
{
  //if (narg < 6) utils::missing_cmd_args(FLERR, "fix setforce", error);

  lower_boundary = 0.0;
  upper_boundary = 10.0;
  width = 0.1;

  hill_weight = 1.0;
  hill_width = 2.0;
  new_hill_freq = 1000;
  output_freq = 100;

  lower_walls = 0.5;
  upper_walls = 29.5;
  force_constant = 10.0;

  if (hill_width > 0.0)
    colvar_sigma = width * hill_width / 2.0;

  hills_grid_size = floor(upper_boundary-lower_boundary)/width;
  memory->create(hills_grid,hills_grid_size,2,"metadynamics:hills_grid");

  count = group->count(igroup);

  MathEigen::Alloc2D(count, 3, &x_group);
  MathEigen::Alloc2D(count, 3, &aaXf_shifted);
  MathEigen::Alloc2D(count, 3, &aaXm_shifted);
  memory->create(refPositions,count,3,"metadynamics:refPositions");
  memory->create(group_index,count,"metadynamics:group_index");

  // FIXME
  // error->all(FLERR, "Error: the number {} of reference positions provided does not match the number {} of atoms of group {}", count, group->names[igroup]);

  refPositions[0][0] = 1.0;
  refPositions[0][1] = 2.0;
  refPositions[0][2] = 3.0;

  refPositions[1][0] = 4.0;
  refPositions[1][1] = 5.0;
  refPositions[1][2] = 6.0;

  refPositions[2][0] = 7.0;
  refPositions[2][1] = 8.0;
  refPositions[2][2] = 9.0;

}

/* ---------------------------------------------------------------------- */

FixMetadynamics::~FixMetadynamics()
{
  if (copymode) return;
  MathEigen::Dealloc2D(&x_group);
  MathEigen::Dealloc2D(&aaXf_shifted);
  MathEigen::Dealloc2D(&aaXm_shifted);
  memory->destroy(refPositions);
  memory->destroy(group_index);
}

/* ---------------------------------------------------------------------- */

int FixMetadynamics::setmask()
{
  int mask = 0;
  //mask |= INITIAL_INTEGRATE;
  mask |= POST_FORCE;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMetadynamics::setup(int vflag)
{
  if (utils::strmatch(update->integrate_style, "^verlet"))
    post_force(vflag);
  else
    error->all(FLERR, "fix metadynamics does not support RESPA.");
}

void FixMetadynamics::initial_integrate(int /*vflag*/)
{
  std::cerr << "-------- STEP " << update->ntimestep << " --------\n";
}

/* ---------------------------------------------------------------------- */

void FixMetadynamics::post_force(int /*vflag*/)
{
  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;

  int j=0;



  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      const tagint k = atom->map(i+1);
      domain->unmap(x[k],image[k],x_group[j]);
      group_index[j]=k;
      j++;
    }

  double inverse_quat[4];

  colvar_value = rmsd(x_group,refPositions,count,inverse_quat);

  // add new biasing energy/forces
  update_bias();

  if (use_grids) {
    // update grid content to reflect new bias
    update_grid_data();
  }

  int grid_index;

  if( use_grids && (grid_index=(colvar_value-lower_boundary)/width-0.5)<hills_grid_size ) {
    bias_energy = hills_grid[grid_index][0];
    colvar_force = hills_grid[grid_index][1];
    //std::cerr << fmt::format("   *** hills_grid[{}] {} {}\n", grid_index,hills_grid[grid_index][0],hills_grid[grid_index][1] );
  } else {
    bias_energy = calc_energy(colvar_value);
    colvar_force = calc_force(colvar_value);
  }

  double drmsddx = (colvar_value>0.0) ? 1.0/(colvar_value*count) : 0.0;

  double quat[4];
  qconjugate(inverse_quat,quat);

  double mat[3][3],inverse_mat[3][3];
  quat_to_mat(quat,mat);
  //write3(mat);
  //quat_to_mat(inverse_quat,inverse_mat);
  //write3(inverse_mat);


  for( int i=0 ; i<count ; i++ ) {

    double grad[3], grad_rot[3];
    j=group_index[i];
    grad[0] = drmsddx*(x_group[i][0]-refPositions[i][0]);
    grad[1] = drmsddx*(x_group[i][1]-refPositions[i][1]);
    grad[2] = drmsddx*(x_group[i][2]-refPositions[i][2]);

    // quaternion rotation of vector: c = a*b*conj(a)
    quatrotvec(quat, grad, grad_rot);

    f[j][0] += colvar_force * grad_rot[0];
    f[j][1] += colvar_force * grad_rot[1];
    f[j][2] += colvar_force * grad_rot[2];

    //std::cerr << fmt::format(" *** FixMetadynamics x_group[{}] {:.6} {:.6} {:.6} grad {:.6} {:.6} {:.6} f[{}] {:.6} {:.6} {:.6}\n", i,x_group[i][0],x_group[i][1], x_group[i][2],grad[0],grad[1],grad[2], j,f[j][0],f[j][1],f[j][2] );

  }

  //std::cerr << fmt::format(" *** FixMetadynamics colvar_value {:.6} bias_energy {:.6} colvar_force {:.6} quat {:.6} {:.6} {:.6} {:.6}\n", colvar_value,bias_energy,colvar_force, quat[0],quat[1],quat[2],quat[3]);

}

/* ---------------------------------------------------------------------- */

void FixMetadynamics::end_of_step()
{

  colvar_history.push_back(std::make_pair(update->ntimestep,colvar_value));

  if((update->ntimestep % output_freq == 0)&&(update->ntimestep > update->firststep)) {

    char traj[] = "output.traj";
    FILE *fp = fopen(traj, "w");
    fmt::print(fp,"{:>18} {:>6}\n", "STEP", "RMSD");

    if (fp == nullptr)
      error->one(FLERR, "Cannot open traj file {}: {}", traj, utils::getsyserror());

    std::for_each(colvar_history.begin(), colvar_history.end(),
      [fp](const std::pair<bigint, double> p) {
        fmt::print(fp,"{:18} {:6f}\n", std::get<0>(p), std::get<1>(p)); });

    if (comm->me == 0) fclose(fp);

    char pmf[] = "output.pmf";

    if( update->ntimestep == output_freq ) {
      fp = fopen(pmf, "w");
      fmt::print(fp,"{:>18} {:>6} {:6}\n", "STEP", "RMSD", "PMF");
    } else
      fp = fopen(pmf, "a");

    if (fp == nullptr)
      error->one(FLERR, "Cannot open pmf file {}: {}", pmf, utils::getsyserror());

    double energy_max = 0.0;

    for( int i=0; i<hills_grid_size ; i++ )
      if( hills_grid[i][0]>energy_max )
        energy_max = hills_grid[i][0];

    for( int i=0; i<hills_grid_size ; i++ ) {
      double v = lower_boundary+((double)i+0.5)*width;
      fmt::print(fp,"{:>18} {:>6f} {:6f}\n", update->ntimestep, v, hills_grid[i][0]-energy_max);
    }

    if (comm->me == 0) fclose(fp);

  }

}

/* ---------------------------------------------------------------------- */

// FIXME: should it be static inline, or inline. do i care...

using std::fdim;
using std::sqrt;

double FixMetadynamics::rmsd( double **aaXf, double **aaXm, int N, double inverse_quat[4])
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

  double Evl[4];                 // Store the eigenvalues of P here.
  double *Evc[4];                // Store the eigevectors here. This version has ** format.
  double _Evc[4 * 4];            // Contiguous 1D array for storing contents of "Evc" array
  for (int i = 0; i < 4; i++)    // We must make sure that
    Evc[i] = &(_Evc[4 * i]);     // Evc[i] points to the correct location in memory


 //for (int i = 0; i < 4; i++) std::cerr << fmt::format("{} {} {} {}\n", P[i][0],P[i][1],P[i][2],P[i][3]);

  int ierror = MathEigen::Jacobi<double, double *, double **>(4).Diagonalize(P, Evl, Evc);

  if(ierror)
    error->all(FLERR, "fix metadynamics: Too many iterations in jacobi diagonalization.\n"
      "This is usually the result of an ill-defined set of atoms for "
      "rotational alignment (RMSD, rotateReference, etc).\n");

  // Note: The eigenvalues are sorted in decreasing order by default.
  double pPp = Evl[0];    // = the maximum eigenvalue of P


  for (int i = 0; i < 4; i++)
    p[i] = Evc[0][i];    //copy eigenvector corresponding to this eigenvalue to p

  // Now normalize p
  double pnorm = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2] + p[3]*p[3]);
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

  inverse_quat[0] =  p[3]; // =  q[0]
  inverse_quat[1] = -p[0]; // = -q[1]
  inverse_quat[2] = -p[1]; // = -q[2]
  inverse_quat[3] = -p[2]; // = -q[3]

  // Finally compute the RMSD between the two coordinate sets:
  // First compute E0 from equation 24 of the paper
  double E0 = 0.0;
  for (size_t n = 0; n < N; n++)
    for (int d = 0; d < 3; d++)
      E0 += (square(aaXf_shifted[n][d] - aaXm_shifted[n][d]));

  for (int n = 0; n < N; n++) {

    //std::cerr << fmt::format(" * STEP0 aaXf[{}] {:.6} {:.6} {:.6}\n", n,aaXf[n][0],aaXf[n][1],aaXf[n][2]);

    //std::cerr << fmt::format(" * STEP1 aaXf[{}] {:.6} {:.6} {:.6}\n", n,aaXf_shifted[n][0],aaXf_shifted[n][1],aaXf_shifted[n][2]);

    double tmp[3];
    // quaternion rotation of vector: c = a*b*conj(a)
    quatrotvec(inverse_quat, aaXf_shifted[n], tmp);

    //std::cerr << fmt::format(" * STEP2 aaXf[{}] {:.6} {:.6} {:.6}\n", n,tmp[0],tmp[1],tmp[2]);

    aaXf[n][0] = tmp[0]+aCenter_m[0];
    aaXf[n][1] = tmp[1]+aCenter_m[1];
    aaXf[n][2] = tmp[2]+aCenter_m[2];

    //std::cerr << fmt::format(" * STEP3 aaXf[{}] {:.6} {:.6} {:.6}\n", n,aaXf[n][0],aaXf[n][1],aaXf[n][2]);

  }

  return sqrt(fdim(E0, 2.0 * pPp)/N);
}

bool FixMetadynamics::can_accumulate_data()
{
  return( update->ntimestep > update->firststep );
}

void FixMetadynamics::update_bias()
{
  // add a new hill if the required time interval has passed
  if ((update->ntimestep % new_hill_freq == 0) && can_accumulate_data() ) {

    double hills_scale=1.0;

    hill h = hill(
      update->ntimestep,
      hill_weight*hills_scale,
      colvar_value,
      colvar_sigma);

    hills.push_back(h);

    // FIXME: if there is more than one replica,
    // communicate it to the others

  }
}

void FixMetadynamics::update_grid_data()
{
  if ((update->ntimestep % new_hill_freq) == 0) {

    // memset is for C coders back in 1990s
    std::fill_n(hills_grid[0],hills_grid_size*2,0.0);

    for( int i=0; i<hills_grid_size ; i++ ) {
      double x = lower_boundary+((double)i+0.5)*width;
      hills_grid[i][0] = calc_energy(x);
      hills_grid[i][1] = calc_force(x);
      //std::cerr << fmt::format("  *** hills_grid[{}] {} {}\n", i,hills_grid[i][0],hills_grid[i][1] );
    }
  }
}

double FixMetadynamics::calc_energy(double x)
{
  double energy = 0.0;

  for (hill_iter h = hills.begin(); h != hills.end(); h++) {

    // compute the gaussian exponent
    double cv_sqdev = square(x - h->center) / square(h->sigma);

    // compute the gaussian
    if (cv_sqdev > 23.0) {
      // set it to zero if the exponent is more negative than log(1.0E-06)
      h->value(0.0);
    } else {
      h->value(exp(-0.5*cv_sqdev));
    }
    energy += h->energy();
    //std::cout << fmt::format("   *** colvar_value {:.6} h->center {:.6} h->sigma {:.6} cv_sqdev {:.6} h->value {:.6} h->energy {:.6} bias_energy {:.6}\n",colvar_value,h->center,h->sigma,cv_sqdev,h->value(),h->energy(),energy);
  }

  return energy;
}

double FixMetadynamics::calc_force(double x)
{
  double force = 0.0;

  for (hill_iter h = hills.begin(); h != hills.end(); h++) {
    if (h->value() == 0.0) continue;
    force += ( h->weight() * h->value() * (0.5 / square(h->sigma)) *
        (2.0 * (x - h->center)) );
  }

  return force;
}
