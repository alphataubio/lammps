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

#include "write_reaxff.h"

#include "pair.h"

#include <cctype>
#include <cstring>

using namespace LAMMPS_NS;

enum { REGULAR_MODE, CLASS2_MODE };

static constexpr int BUF_SIZE = 256;

/* ----------------------------------------------------------------------
   called as write_coeff command in input script
------------------------------------------------------------------------- */

void WriteReaxff::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all(FLERR, "Write_reaxff command before simulation box is defined");

  if (narg != 1) utils::missing_cmd_args(FLERR, "write_reaxff", error);

  char *file = utils::strdup(fmt::format("{}.tmp", arg[0]));

  // initialize relevant styles
  lmp->init();

  if (comm->me == 0)
    Write_Force_Field(arg[0], &(api->system->reax_param), api->control, world);

}
