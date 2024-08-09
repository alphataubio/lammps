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
  void setup(int) override;
  void post_force(int) override;

  // rmsd
  int N;
  double **aaXf_shifted, **aaXm_shifted;

  double rmsd(double **,double **,int);

 protected:

  // colvar
  double lowerBoundary, upperBoundary, width;

  // metadynamics
  double hillWeight, hillWidth;
  int newHillFrequency, outputFreq;

  double *hill_centers;

  // harmonic walls
  double lowerWalls, upperWalls, forceConstant;



  // fix setforce
  char *xstr, *ystr, *zstr;
  char *idregion;
  class Region *region;

/*
  double xvalue, yvalue, zvalue;
  int varflag;
  int xvar, yvar, zvar, xstyle, ystyle, zstyle;
  double foriginal[3], foriginal_all[3], foriginal_saved[3];
  int force_flag;
  int nlevels_respa, ilevel_respa;

  int maxatom;
  double **sforce;
*/

};

}    // namespace LAMMPS_NS

#endif
#endif
