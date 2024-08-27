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
FixStyle(bond/react/sbml,FixBondReactSbml);
// clang-format on
#else

#ifndef LMP_FIX_BOND_REACT_SBML_H
#define LMP_FIX_BOND_REACT_SBML_H

#include "fix_bond_react.h"

namespace LAMMPS_NS {

//class FixBondReactSbml : public FixBondReact {
class FixBondReactSbml : public Fix {
 public:
  FixBondReactSbml(class LAMMPS *, int, char **);
  ~FixBondReactSbml() override;

  int setmask() override;

 protected:


};

}    // namespace LAMMPS_NS

#endif
#endif
