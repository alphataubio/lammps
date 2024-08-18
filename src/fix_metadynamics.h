/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(metadynamics,FixMetadynamics);
// clang-format on
#else

#ifndef LMP_FIX_METADYNAMICS_H
#define LMP_FIX_METADYNAMICS_H

#include "fix.h"

namespace LAMMPS_NS {

class FixMetadynamics : public Fix {
 public:
  FixMetadynamics(class LAMMPS *, int, char **);
  ~FixMetadynamics() override;
  int setmask() override;
  void init() override;
  void initial_integrate(int) override;
  void post_force(int) override;
  void end_of_step() override;
  double compute_scalar() override;

  double **x_group, **x_group_shifted, **ref_positions_shifted;

  double rmsd(double *);

 protected:

  tagint *group_taglist;

  static int idcompare(const int, const int, void *);

  // -------- RMSD --------

  bigint group_count;
  double **ref_positions;

  // -------- COLVAR --------

  // Location of the lower boundary
  double lower_boundary;

  // Location of the upper boundary
  double upper_boundary;

  // Resolution of colvar
  double width;

  // Up to date value of colvar
  double colvar_value;

  // Current energy of colvar
  double colvar_energy;

  // Current force from this colvar
  double colvar_force;

  // -------- HARMONIC UPPER WALL --------

  double upper_wall_force_constant;

  // -------- HILLS --------

  // Height of new hills
  double hill_weight;

  // The local width of each collective variable, multiplied by this
  // number, provides the hill width along that direction
  double hill_width;

  // The sigma parameters of the Gaussian hills
  double hill_sigma;

  int number_hills;
  double *hill_centers;
  
  // Number of simulation steps between two hills
  int new_hill_freq;

  // Bin the hills on grids of energy and forces, and use them
  // to force the colvars (as opposed to deriving the hills analytically)
  bool use_grids = true;

  // Hills grid cache *[0] = energy, *[1]=force
  int hills_grid_size;
  double **hills_grid;

  // -------- MULTIPLE REPLICAS --------
  
  // Frequency at which data the "mirror" biases are updated
  int replica_update_freq;
  
  // MPI comm with 1 root proc from each world
  MPI_Comm roots;
  
  // -------- OUTPUT FILES --------

  // trajectory and pmf output
  FILE *fp_traj, *fp_pmf;

  // Frequency for writing pmf file
  int output_freq;

  // -------- PRIVATE IMPLEMENTATION METHODS --------

  void read_xyz(char *);
  void update_hills();
  void update_replicas();
  void calc_energy_and_force(double, double &, double &);

};
}    // namespace LAMMPS_NS


#endif
#endif
