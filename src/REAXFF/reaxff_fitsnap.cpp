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

#include "reaxff_api.h"

#include "error.h"
#include "memory.h"
#include "text_file_reader.h"
#include "utils.h"

#include <cmath>
#include <cstring>
#include <exception>
#include <string>

using LAMMPS_NS::utils::open_potential;
using LAMMPS_NS::utils::getsyserror;
using LAMMPS_NS::utils::strmatch;
using LAMMPS_NS::utils::uppercase;
using LAMMPS_NS::EOFException;
using LAMMPS_NS::ValueTokenizer;

namespace ReaxFF {



#undef THROW_ERROR

}
