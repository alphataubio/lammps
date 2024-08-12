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

#include <array>
#include <list>
#include <vector>

namespace LAMMPS_NS {

class hill;
typedef std::list<hill>::iterator hill_iter;

class FixMetadynamics : public Fix {
 public:
  FixMetadynamics(class LAMMPS *, int, char **);
  ~FixMetadynamics() override;
  int setmask() override;
  void setup(int) override;
  void initial_integrate(int) override;
  void post_force(int) override;
  void end_of_step() override;

  // rmsd
  bigint count;
  double **x_group, **aaXf_shifted, **aaXm_shifted;

  double rmsd(double **,double **,int,double *);

 protected:

  // colvar
  double lower_boundary, upper_boundary, width;
  double **refPositions;

  std::vector<std::pair<bigint, double>> colvar_history;

  // metadynamics
  int new_hill_freq, output_freq;
  int *group_index;
  double hill_weight;

  // harmonic walls
  double lower_walls, upper_walls, force_constant;

 private:

  bool use_grids = true;

  // The local width of each collective variable, multiplied by this
  // number, provides the hill width along that direction
  double hill_width;

  // The sigma parameters of the Gaussian hills
  double colvar_sigma;

  // Up to date value of each colvar
  double colvar_value;

  // Current energy of this bias (colvar_forces should be obtained by deriving this)
  double bias_energy;

  // Current forces from this bias to the variables
  double colvar_force;

  // List of hills used on this bias (total); if a grid is
  // employed, these don't need to be updated at every time step
  std::list<hill> hills;

  // colvarsbias_meta
  void update_bias();
  void update_grid_data();

  // Hill energy, cached on a grid *[0]
  // Hill forces, cached on a grid *[1]
  int hills_grid_size;
  double **hills_grid;

  double calc_energy(double);
  double calc_force(double);

  // Returns true if the current step represent a valid increment, whose data
  // can be recorded (as opposed to e.g. a continuation step from a restart)
  bool can_accumulate_data();

  bool debug();

  char name[5] = "RMSD";

  //  Add a new hill; if a .hills trajectory is written,
  // write it there; 
  void add_hill(hill const &h);


};


// A hill for the metadynamics bias
class hill {

protected:

  /// Time step at which this hill was added
  bigint it;

  /// Value of the hill function (ranges between 0 and 1)
  double hill_value;

  /// Scale factor, which could be modified at runtime (default: 1)
  double sW;

  /// Maximum height in energy of the hill
  double W;

  /// Centers of the hill in the collective variable space
  double center;

  /// Half-widths of the hill in the collective variable space
  double sigma;

  /// Identity of the replica who added this hill
  std::string replica;
  
public:

  friend class FixMetadynamics;

  /// Constructor of a hill object
  /// \param it Step number at which the hill was added
  /// \param W Weight of the hill (energy units)
  /// \param cv_values Array of collective variable values
  /// \param cv_sigmas Array of collective variable values
  /// \param replica ID of the replica that creates the hill (optional)
  hill(bigint it_in, double W_in, double center_in, double sigma_in) :
    it(it_in), hill_value(0.0), sW(1.0), W(W_in), center(center_in), sigma(sigma_in) {}

  /// Copy constructor
  hill(hill const &h) : it(h.it), hill_value(0.0), sW(1.0), W(h.W), center(h.center),
    sigma(h.sigma) {}

  /// Destructor
  ~hill() {}

  /// Assignment operator
  hill & operator = (hill const &h);

  /// Get the energy
  inline double energy()
  {
    return W * sW * hill_value;
  }

  /// Get the energy using another hill weight
  inline double energy(double const &new_weight)
  {
    return new_weight * sW * hill_value;
  }

  /// Get the current hill value
  inline double const &value()
  {
    return hill_value;
  }

  /// Set the hill value as specified
  inline void value(double const &new_value)
  {
    hill_value = new_value;
  }

  /// Get the weight
  inline double weight()
  {
    return W * sW;
  }

  /// Scale the weight with this factor (by default 1.0 is used)
  inline void scale(double const &new_scale_fac)
  {
    sW = new_scale_fac;
  }

  /// Comparison operator
  inline friend bool operator < (hill const &h1, hill const &h2)
  {
    if (h1.it < h2.it) return true;
    else return false;
  }

  /// Comparison operator
  inline friend bool operator <= (hill const &h1, hill const &h2)
  {
    if (h1.it <= h2.it) return true;
    else return false;
  }

  /// Comparison operator
  inline friend bool operator > (hill const &h1, hill const &h2)
  {
    if (h1.it > h2.it) return true;
    else return false;
  }

  /// Comparison operator
  inline friend bool operator >= (hill const &h1, hill const &h2)
  {
    if (h1.it >= h2.it) return true;
    else return false;
  }

  /// Comparison operator
  inline friend bool operator == (hill const &h1, hill const &h2)
  {
    if ( (h1.it >= h2.it) && (h1.replica == h2.replica) ) return true;
    else return false;
  }

  /// Represent the hill ina string suitable for a trajectory file
  std::string output_traj();

};

}    // namespace LAMMPS_NS


#endif
#endif
