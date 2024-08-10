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
#include "math_eigen_impl.h"
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
using namespace MathSpecial;

static constexpr int HILLS_DELTA = 16;

/* ---------------------------------------------------------------------- */

FixMetadynamics::FixMetadynamics(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg)
{
  //if (narg < 6) utils::missing_cmd_args(FLERR, "fix setforce", error);

  lower_boundary = 0.0;
  upper_boundary = 10.0;
  width = 1.0;

  hill_weight = 0.1;
  hill_width = 2.0;
  new_hill_freq = 5;
  output_freq = 100;

  lower_walls = 0.5;
  upper_walls = 29.5;
  force_constant = 10.0;

  if (hill_width > 0.0)
    colvar_sigma = width * hill_width / 2.0;

  hills_grid_size = floor(upper_boundary-lower_boundary)/width;
  memory->create(hills_grid,hills_grid_size,2,"metadynamics:hills_grid");

  N=3;

  MathEigen::Alloc2D(N, 3, &aaXf_shifted);
  MathEigen::Alloc2D(N, 3, &aaXm_shifted);

  memory->create(refPositions,3,3,"metadynamics:refPositions");

}

/* ---------------------------------------------------------------------- */

FixMetadynamics::~FixMetadynamics()
{
  if (copymode) return;

  memory->destroy(refPositions);

}

/* ---------------------------------------------------------------------- */

int FixMetadynamics::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
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


/* ---------------------------------------------------------------------- */

void FixMetadynamics::post_force(int /*vflag*/)
{
  double **x = atom->x;
  //double **f = atom->f;
  //int *mask = atom->mask;
  //int nlocal = atom->nlocal;

  // Update the cached colvar values
  //for (size_t i = 0; i < num_variables(); i++) {
  //  colvar_values[i] = colvars[i]->value();
  //}

  //domain->remap(x[0]);
  //domain->remap(x[1]);
  //domain->remap(x[2]);
  domain->minimum_image_big(x[0][0],x[0][1],x[0][2]);
  domain->minimum_image_big(x[1][0],x[1][1],x[1][2]);
  domain->minimum_image_big(x[2][0],x[2][1],x[2][2]);
  //domain->pbc();

  colvar_value = rmsd(x,refPositions,3);

  // add new biasing energy/forces
  update_bias();

  if (use_grids) {
    // update grid content to reflect new bias
    update_grid_data();
    int i = 2.0*(colvar_value-lower_boundary)-1.0;
    bias_energy = hills_grid[i][0];
    colvar_force = hills_grid[i][1];
  } else {
    bias_energy = calc_energy(colvar_value);
    colvar_force = calc_force(colvar_value);
  }

  std::cout << fmt::format("   *** x[0] {:.6} {:.6} {:.6} colvar_value {:.6} bias_energy {:.6} colvar_force {:.6}\n", x[0][0],x[0][1],x[0][2],colvar_value,bias_energy,colvar_force );

}

bool FixMetadynamics::can_accumulate_data()
{
  return( (update->ntimestep - update->firststep) > 0 );
}

bool FixMetadynamics::debug() { return true; }

void FixMetadynamics::update_bias()
{

  // add a new hill if the required time interval has passed
  if (((update->ntimestep % new_hill_freq) == 0) && can_accumulate_data() ) {

    double hills_scale=1.0;

    hill h = hill(
      update->ntimestep,
      hill_weight*hills_scale,
      colvar_value,
      colvar_sigma);

    hills.push_back(h);

    // FIXME: if there is more than one replica,
    // communicate it to the others

    if (debug())
      utils::logmesg(lmp,"   *** Metadynamics bias \"{}\": adding a new hill at step {} W {:.6} center {:.6} sigma {:.6}.\n", this->name,update->ntimestep,h.W,h.center,h.sigma);

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
      std::cerr << fmt::format("   *** hills_grid[{}] {} {} {}\n",
        i,x,hills_grid[i][0],hills_grid[i][1]);
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
