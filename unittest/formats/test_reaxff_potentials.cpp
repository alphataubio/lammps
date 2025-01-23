/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */


#include "info.h"
#include "input.h"
#include "potential_file_reader.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "../testing/core.h"

#include <cstring>
//#include <filesystem>
#include <iostream>
#include <mpi.h>
#include <vector>
#include <dirent.h>



using namespace LAMMPS_NS;
using utils::split_words;
//namespace fs = std::filesystem;

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;


class ReaxffTest : public LAMMPSTest {

protected:

  std::string first_atom(char *filename) {

    PotentialFileReader reader(lmp, filename, "potential", "reaxff");
    reader.skip_line();
    int num_general_parameters = reader.next_int();

    for( int i=0 ; i<num_general_parameters+4 ; i++ )
      reader.skip_line();

    return reader.next_string();

  }

};

// open for native units
TEST_F(ReaxffTest, ForceFields)
{

    if (!Info::has_package("REAXFF")) GTEST_SKIP();

    //BEGIN_HIDE_OUTPUT();

    command("echo none");
    command("units real");
    command("atom_style charge");
    command("region box block 0 2 0 2 0 2");
    command("create_box 1 box");
    command("mass 1 1.0");

    command("pair_style reaxff NULL checkqeq no");

    const char *folder = getenv("LAMMPS_POTENTIALS");
    fprintf(stderr, "LAMMPS_POTENTIALS = %s\n", folder);
/*
    if (folder != nullptr)
      for (auto const& dir_entry : fs::directory_iterator(folder)) {
        fprintf(stderr, "%s\n", dir_entry.path());
      }
*/

DIR *dir;
struct dirent *ent;
if ((dir = opendir (folder)) != NULL) {
  /* print all the files and directories within directory */
  while ((ent = readdir (dir)) != NULL)
    if( strncmp(ent->d_name,"reaxff-",7)==0 ){

      try {
        fprintf (stderr, "%s ", ent->d_name);
        std::string atom = first_atom(ent->d_name);
        fprintf (stderr, "%s\n", atom.c_str());
        auto cmd = fmt::format("pair_coeff * * {} {}", ent->d_name, atom);
        command(cmd);
      //} catch (LAMMPSException &e) {
      } catch (...) {



      }

  }
  closedir (dir);
} else {
  /* could not open directory */
  perror ("");
  //return EXIT_FAILURE;
}

    //END_HIDE_OUTPUT();

}



int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleMock(&argc, argv);

    // handle arguments passed via environment variable
    if (const char *var = getenv("TEST_ARGS")) {
        std::vector<std::string> env = split_words(var);
        for (auto arg : env) {
            if (arg == "-v") {
                verbose = true;
            }
        }
    }
    if ((argc > 1) && (strcmp(argv[1], "-v") == 0)) verbose = true;

    int rv = RUN_ALL_TESTS();
    MPI_Finalize();
    return rv;
}
