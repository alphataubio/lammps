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

/* ----------------------------------------------------------------------
   Contributing author: Mitch Murphy (alphataubio at gmail)
------------------------------------------------------------------------- */

#include "bond_scan.h"

#include "atom.h"
#include "input.h"
#include "variable.h"


using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

BondScan::BondScan(LAMMPS *lmp) : Command(lmp) {}

/* ---------------------------------------------------------------------- */

void BondScan::command(int narg, char **arg)
{
  if (narg < 1) utils::missing_cmd_args(FLERR, "bond_scan", error);
  if (atom->tag_enable == 0) error->all(FLERR, "Cannot use bond_scan unless atoms have IDs");

}


