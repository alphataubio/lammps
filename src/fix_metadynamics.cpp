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
#include "math_eigen_impl.h"
#include "math_extra.h"
#include "math_special.h"
#include "memory.h"
#include "text_file_reader.h"
#include "update.h"

#include <cmath>

#include <iostream>


using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathExtra;
using namespace MathSpecial;

/* ---------------------------------------------------------------------- */

FixMetadynamics::FixMetadynamics(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg)
{

  scalar_flag = 1;
  extscalar = 0;

  if (narg < 4) utils::missing_cmd_args(FLERR, "fix metadynamics", error);

  // xyz_file
  group_count = group->count(igroup);
  read_xyz(arg[3]);

  // traj and pmf output files

  char traj[] = "output.traj";
  char pmf[] = "output.pmf";

  fp_traj = fopen(traj, "w");
  //fmt::print(fp,"{:>18} {:>6}\n", "STEP", "RMSD");

  if (fp_traj == nullptr)
    error->one(FLERR, "Cannot open traj file {}: {}", traj, utils::getsyserror());

  fp_pmf = fopen(pmf, "w");
  //fmt::print(fp,"{:>18} {:>6} {:6}\n", "STEP", "RMSD", "PMF");

  if (fp_pmf == nullptr)
    error->one(FLERR, "Cannot open pmf file {}: {}", pmf, utils::getsyserror());

  // default options

  lower_boundary = 0.0;
  upper_boundary = 10.0;
  width = 0.1;

  hill_weight = 1.0;
  hill_width = 2.0;
  new_hill_freq = 1000;
  output_freq = 1000;

  upper_wall_force_constant = 2.0;

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "lowerBoundary") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix metadynamics lowerBoundary", error);
      lower_boundary = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      //if (nevery <= 0) error->all(FLERR, "Invalid fix addforce every argument: {}", nevery);
      iarg += 2;
    } else if (strcmp(arg[iarg], "upperBoundary") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix metadynamics upperBoundary", error);
      upper_boundary = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      //if (nevery <= 0) error->all(FLERR, "Invalid fix addforce every argument: {}", nevery);
      iarg += 2;
    } else if (strcmp(arg[iarg], "width") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix metadynamics width", error);
      width = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      //if (nevery <= 0) error->all(FLERR, "Invalid fix addforce every argument: {}", nevery);
      iarg += 2;
    } else if (strcmp(arg[iarg], "hillWeight") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix metadynamics hillWeight", error);
      hill_weight = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      //if (nevery <= 0) error->all(FLERR, "Invalid fix addforce every argument: {}", nevery);
      iarg += 2;
    } else if (strcmp(arg[iarg], "hillWidth") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix metadynamics hillWidth", error);
      hill_width = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      //if (nevery <= 0) error->all(FLERR, "Invalid fix addforce every argument: {}", nevery);
      iarg += 2;
    } else if (strcmp(arg[iarg], "newHillFrequency") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix metadynamics newHillFrequency", error);
      new_hill_freq = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      //if (nevery <= 0) error->all(FLERR, "Invalid fix addforce every argument: {}", nevery);
      iarg += 2;
    } else if (strcmp(arg[iarg], "upperWallForceConstant") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix metadynamics upperWallForceConstant", error);
      upper_wall_force_constant = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      //if (nevery <= 0) error->all(FLERR, "Invalid fix addforce every argument: {}", nevery);
      iarg += 2;
    }
  }

  number_hills = 0;
  hill_sigma = width * hill_width / 2.0;
  hills_grid_size = floor(upper_boundary-lower_boundary)/width;
  memory->create(hills_grid,hills_grid_size,2,"metadynamics:hills_grid");

  MathEigen::Alloc2D(group_count, 3, &x_group);
  MathEigen::Alloc2D(group_count, 3, &x_group_shifted);
  MathEigen::Alloc2D(group_count, 3, &ref_positions_shifted);
  memory->create(group_taglist,group_count,"metadynamics:group_taglist");
}

/* ---------------------------------------------------------------------- */

FixMetadynamics::~FixMetadynamics()
{
  if (copymode) return;
  MathEigen::Dealloc2D(&x_group);
  MathEigen::Dealloc2D(&x_group_shifted);
  MathEigen::Dealloc2D(&ref_positions_shifted);
  memory->destroy(hill_centers);
  memory->destroy(ref_positions);
  memory->destroy(group_taglist);
  if (comm->me == 0) {
    fclose(fp_traj);
    fclose(fp_pmf);
  }
}

/* ---------------------------------------------------------------------- */

int FixMetadynamics::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  //mask |= POST_NEIGHBOR;
  mask |= POST_FORCE;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMetadynamics::init()
{
  utils::logmesg(lmp,"FixMetadynamics:\n  firststep {} laststep {} beginstep {} endstep {} ntimestep {} nsteps {}\n  lowerBoundary {} upperBoundary {} width {}\n  hillWeight {} hillWidth {} newHillFreq {} outputFreq {}\n",
      update->firststep, update->laststep, update->beginstep, update->endstep, update->ntimestep, update->nsteps,
      lower_boundary, upper_boundary, width, hill_weight, hill_width, new_hill_freq, output_freq );

  memory->create(hill_centers,ceil(update->nsteps / new_hill_freq),"metadynamics:hill_centers");

  int *mask = atom->mask;
  for( int i=0, j=0 ; i < atom->nlocal ; i++ )
    if (mask[i] & groupbit) {
      group_taglist[j] = atom->tag[i];
      //std::cerr << fmt::format("group_taglist[{}] {}\n", j,group_taglist[j]);
      j++;
    }

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

void FixMetadynamics::post_force(int /*vflag*/)
{

  double inverse_quat[4];
  colvar_value = rmsd(inverse_quat);
  update_hills();

  int grid_index = (int) floor((colvar_value-lower_boundary)/width);

  if( use_grids && grid_index<hills_grid_size ) {
    colvar_energy = hills_grid[grid_index][0];
    colvar_force = hills_grid[grid_index][1];
    //std::cerr << fmt::format("   *** hills_grid[{}] {} {}\n", grid_index,hills_grid[grid_index][0],hills_grid[grid_index][1] );
  } else {
    calc_energy_and_force(colvar_value,colvar_energy,colvar_force);
  }

  double **f = atom->f;
  double drmsddx = (colvar_value>0.0) ? 1.0/(colvar_value*group_count) : 0.0;
  double quat[4];
  qconjugate(inverse_quat,quat);

  for( int n=0 ; n<group_count ; n++ ) {

    double grad[3], grad_rot[3];
    const int i = atom->map(group_taglist[n]);
    grad[0] = drmsddx*(x_group[n][0]-ref_positions[n][0]);
    grad[1] = drmsddx*(x_group[n][1]-ref_positions[n][1]);
    grad[2] = drmsddx*(x_group[n][2]-ref_positions[n][2]);

    // quaternion rotation of vector: c = a*b*conj(a)
    quatrotvec(quat, grad, grad_rot);

    f[i][0] += colvar_force * grad_rot[0];
    f[i][1] += colvar_force * grad_rot[1];
    f[i][2] += colvar_force * grad_rot[2];

    std::cerr << fmt::format(" *** x_group[{}] {:.6} {:.6} {:.6} grad {:.6} {:.6} {:.6} f[{}] {:.6} {:.6} {:.6}\n", n,x_group[n][0],x_group[n][1], x_group[n][2],grad[0],grad[1],grad[2], i,f[i][0],f[i][1],f[i][2] );

  }

  //double mat[3][3],inverse_mat[3][3];
  //quat_to_mat(quat,mat);
  //write3(mat);
  //quat_to_mat(inverse_quat,inverse_mat);
  //write3(inverse_mat);

  std::cerr << fmt::format(" *** FixMetadynamics colvar_value {:.6} colvar_energy {:.6} colvar_force {:.15} quat {:.6} {:.6} {:.6} {:.6} |quat| {}\n", colvar_value,colvar_energy,colvar_force, quat[0],quat[1],quat[2],quat[3],sqrt(quat[0]*quat[0]+quat[1]*quat[1]+quat[2]*quat[2]+quat[3]*quat[3]));

}


/* ---------------------------------------------------------------------- */

// FIXME: should it be static inline, or inline. do i care...

using std::fdim;
using std::sqrt;

double FixMetadynamics::rmsd( double inverse_quat[4] )
{

  // Find the center-of-geometry of each object:
  double **x = atom->x;
  imageint *image = atom->image;
  double aCenter_f[3] = {0.0};
  double aCenter_m[3] = {0.0};
  for (int n = 0; n < group_count; n++) {
    const int i = atom->map(group_taglist[n]);
    domain->unmap(x[i],image[i],x_group[n]);
    //std::cerr << fmt::format(" *** BEFORE_CG x_group[{}] {:.6} {:.6} {:.6}\n", n,x_group[n][0],x_group[n][1],x_group[n][2]);
    for (int d = 0; d < 3; d++) {
      aCenter_f[d] += x_group[n][d];
      aCenter_m[d] += ref_positions[n][d];
    }
  }

  for (int d = 0; d < 3; d++) {
    aCenter_f[d] /= group_count;
    aCenter_m[d] /= group_count;
  }

  //Subtract the centers-of-geometry from the original coordinates for each object
  for (int n = 0; n < group_count; n++) {
    for (int d = 0; d < 3; d++) {
      // shift the coordinates so that the new center of geometry is at the origin
      x_group_shifted[n][d] = x_group[n][d] - aCenter_f[d];
      ref_positions_shifted[n][d] = ref_positions[n][d] - aCenter_m[d];
    }
        std::cerr << fmt::format(" *** AFTER_CG x_group[{}] {:.6} {:.6} {:.6}\n", n,x_group_shifted[n][0],x_group_shifted[n][1],x_group_shifted[n][2]);
  }

  // Calculate the "M" array from the Diamond paper (equation 16)
  double M[3][3];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) M[i][j] = 0.0;

  for (size_t n = 0; n < group_count; n++) {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) { M[i][j] += ref_positions_shifted[n][i] * x_group_shifted[n][j]; }
    }
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

  for (unsigned i = 0; i < 4; i++) {
    printf("[ ");
    for (unsigned j = 0; j < 4; j++) printf("%g ",P[i][j]);
    printf("]\n");
  }

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
  pnorm = sqrt(pnorm);
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

  for (size_t n = 0; n < group_count; n++) {
    double tmp[3];
    // quaternion rotation of vector: c = a*b*conj(a)
    quatrotvec(inverse_quat, x_group_shifted[n], tmp);
    std::cerr << fmt::format(" *** AFTER_ROT inverse_quat {:.6} {:.6} {:.6} {:.6} x_group[{}] {:.6} {:.6} {:.6}\n", inverse_quat[0],inverse_quat[1],inverse_quat[2],inverse_quat[3], n,tmp[0],tmp[1],tmp[2]);
    for (int d = 0; d < 3; d++) {
      E0 += (square(x_group_shifted[n][d] - ref_positions_shifted[n][d]));
      x_group[n][d] = tmp[d]+aCenter_m[d];
    }
    std::cerr << fmt::format(" *** AFTER_TRANSLATION x_group[{}] {:.6} {:.6} {:.6}\n", n,x_group[n][0],x_group[n][1],x_group[n][2]);
  }
  return sqrt(fdim(E0, 2.0 * Evl[0])/group_count); // Evl[0] = the maximum eigenvalue of P
}













































































/* ---------------------------------------------------------------------- */

void FixMetadynamics::end_of_step()
{
  fmt::print(fp_traj,"{:18} {:6f}\n", update->ntimestep, colvar_value);

  if((update->ntimestep % output_freq == 0)&&(update->ntimestep > update->firststep)) {

    double energy_max = 0.0;

    for( int i=0; i<hills_grid_size ; i++ )
      if( hills_grid[i][0]>energy_max )
        energy_max = hills_grid[i][0];

    for( int i=0; i<hills_grid_size ; i++ ) {
      double v = lower_boundary+((double)i+0.5)*width;
      fmt::print(fp_pmf,"{:>18} {:>6f} {:6f}\n", update->ntimestep, v, hills_grid[i][0]-energy_max);
    }
  }
}

/* ---------------------------------------------------------------------- */

double FixMetadynamics::compute_scalar() { return colvar_value; }

/* ---------------------------------------------------------------------- */

void FixMetadynamics::read_xyz(char *filename)
{

  memory->create(ref_positions,group_count,3,"metadynamics:ref_positions");
  if (comm->me != 0) return;
  FILE *fp = fopen(filename, "r");

  if (fp == nullptr)
    error->one(FLERR, "FixMetadynamics: Cannot open file {}: {}", filename, utils::getsyserror());

  TextFileReader reader(fp, "xyz");
  try {
    char *line = reader.next_line(1);
    int xyz_count = ValueTokenizer(line).next_int();

    if( xyz_count != group_count )
      error->all(FLERR, "FixMetadynamics: number {} of reference positions in {} does not match the number {} of atoms of group {}", xyz_count, filename, group_count, group->names[igroup]);

    reader.skip_line();
    for( int i=0 ; i<xyz_count ; i++ ) {
      double buffer[4];
      reader.next_dvector(buffer,4);
      ref_positions[i][0] = buffer[1];
      ref_positions[i][1] = buffer[2];
      ref_positions[i][2] = buffer[3];
    }
  } catch (std::exception &e) {
    error->all(FLERR, "FixMetadynamics: error reading xyz file {}: {}", filename, e.what());
  }
  fclose(fp);
}

/* ---------------------------------------------------------------------- */

void FixMetadynamics::update_hills()
{
  // add a new hill if the required time interval has passed
  if((update->ntimestep % new_hill_freq == 0)&&(update->ntimestep > update->firststep)) {
    hill_centers[number_hills] = colvar_value;
    number_hills++;

    // FIXME: if there is more than one replica,
    // communicate it to the others

    if (use_grids)
      for( int i=0; i<hills_grid_size ; i++ ) {
        double x = lower_boundary+((double)i+0.5)*width;
        calc_energy_and_force(x,hills_grid[i][0],hills_grid[i][1]);
      //std::cerr << fmt::format(" *** hills_grid[{}] {} {}\n", i,hills_grid[i][0],hills_grid[i][1] );
      }
  }
}

/* ---------------------------------------------------------------------- */

void FixMetadynamics::calc_energy_and_force(double xi, double &energy, double &force)
{
  energy = force = 0.0;

  for( int i=0 ; i<number_hills ; i++ ) {

    // compute the gaussian exponent
    double hill_center = hill_centers[i];
    double cv_sqdev = square(xi - hill_center) / square(hill_sigma);

    // compute the gaussian, only if exponent > log(1.0E-06)
    if (cv_sqdev < 23.0) {
      double gaussian = exp(-0.5*cv_sqdev);
      energy += gaussian;
      force += gaussian * (xi-hill_center) / square(hill_sigma);
    }
  }

  energy *= hill_weight;
  force *= hill_weight;

  if( xi > upper_boundary ) {
    force -= upper_wall_force_constant*(xi - upper_boundary)/square(width);
    std::cerr << fmt::format(" *** xi {:.6} energy {:.6} force {:.6} \n", xi,0.0,
      upper_wall_force_constant*( xi - upper_boundary)/square(width));

  }
  //std::cerr << fmt::format(" *** calc_energy_and_force xi {:.6} energy {:.6} force {:.6} \n", xi,energy,force);

}
