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
Contributing Author: Mitch Murphy (alphataubio at gmail com)
------------------------------------------------------------------------- */

#include "fix_bond_react_sbml.h"

#include "sbml/SBMLTypes.h"

#include <filesystem>
#include <iostream>

using namespace LAMMPS_NS;
using namespace FixConst;

using namespace std;
using namespace libsbml;



/* ---------------------------------------------------------------------- */
// clang-format off

//FixBondReactSbml::FixBondReactSbml(LAMMPS *lmp, int narg, char **arg) :
//  FixBondReact(lmp, narg, arg)

FixBondReactSbml::FixBondReactSbml(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{

  const char* filename   = arg[3];

  cerr << " *** filename " << filename <<"\n";

  SBMLDocument* document = readSBML(filename);

  if (document->getNumErrors() > 0)
  {
    cerr << "Encountered the following SBML errors:" << endl;
    document->printErrors(cerr);
    //return;
  }

  unsigned int level   = document->getLevel  ();
  unsigned int version = document->getVersion();

  cout << endl
       << "File: " << filename
       << " (Level " << level << ", version " << version << ")" << endl;

  Model* model = document->getModel();

  if (model == 0)
  {
    cout << "No model present." << endl;
    return;
  }

  cout << "               "
       << (level == 1 ? "name: " : "  id: ")
       << (model->isSetId() ? model->getId() : std::string("(empty)")) << endl;

  if (model->isSetSBOTerm())
    cout << "      model sboTerm: " << model->getSBOTerm() << endl;

  cout << "functionDefinitions: " << model->getNumFunctionDefinitions() << endl;
  cout << "    unitDefinitions: " << model->getNumUnitDefinitions    () << endl;
  cout << "   compartmentTypes: " << model->getNumCompartmentTypes   () << endl;
  cout << "        specieTypes: " << model->getNumSpeciesTypes       () << endl;
  cout << "       compartments: " << model->getNumCompartments       () << endl;
  cout << "            species: " << model->getNumSpecies            () << endl;
  cout << "         parameters: " << model->getNumParameters         () << endl;
  cout << " initialAssignments: " << model->getNumInitialAssignments () << endl;
  cout << "              rules: " << model->getNumRules              () << endl;
  cout << "        constraints: " << model->getNumConstraints        () << endl;
  cout << "          reactions: " << model->getNumReactions          () << endl;
  cout << "             events: " << model->getNumEvents             () << endl;
  cout << endl;

  delete document;


}

/* ---------------------------------------------------------------------- */

FixBondReactSbml::~FixBondReactSbml()
{
  if (copymode) return; // needed for KOKKOS [alphataubio,2024/08]

}

/* ---------------------------------------------------------------------- */

int FixBondReactSbml::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  mask |= POST_INTEGRATE_RESPA;
  mask |= POST_FORCE;
  return mask;
}
