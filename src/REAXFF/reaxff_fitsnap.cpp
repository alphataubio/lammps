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

void Write_Force_Field(const char *filename, reax_interaction *reax,
                        control_params *control, MPI_Comm world)
{
  char ****tor_flag;
  auto error = control->error_ptr;
  auto lmp = control->lmp_ptr;
  auto memory = control->lmp_ptr->memory;


#define THROW_ERROR(txt) throw ffield_parser_error(fmt::format("{}:{}: {}", filename, lineno, txt))

  // write force field only on rank 0
  if (control->me != 0) return;

  FILE *fp = fopen(filename, "w");
  if (!fp) error->one(FLERR, "Cannot open reaxff {} file", filename);

  auto &gp = reax->gp;

  try {
    int i,j,k,l,m,n,lineno = 0;

    // header comment line
    fmt::print(fp, "Reactive MD-force field: second-generation water force field-2017\n");

    // number of global parameters
    n = gp.n_global;
    fmt::print(fp, " {}       ! Number of general parameters\n", n);

    //memory->create(gp.l,n,"reaxff:gp.l");
    // see reaxff_types.h for mapping between l[i] and the lambdas used in ff

    //for (i = 0; i < n; ++i) {
    //      values = reader.next_values(0);
    //      ++lineno;
    //      gp.l[i] = values.next_double();
    //    }


    fmt::print(fp, "{:10.4f} !p_boc1 Eq(4c): Overcoordination parameter\n", gp.l[0]);
    fmt::print(fp, "{:10.4f} !p_boc2 Eq(4d): Overcoordination parameter\n", gp.l[1]);


    fmt::print(fp, "{:10.4f} !p_coa2 Eq(15): Valency angle conjugation\n", gp.l[2]);
    fmt::print(fp, "{:10.4f} !p_trip4 Eq(20): Triple bond stabilisation\n", gp.l[3]);
    fmt::print(fp, "{:10.4f} !p_trip3 Eq(20): Triple bond stabilisation\n", gp.l[4]);
    fmt::print(fp, "{:10.4f} !k_c2 Eq(19): C2-correction\n", gp.l[5]);
    fmt::print(fp, "{:10.4f} !p_ovun6 Eq(12): Undercoordination\n", gp.l[6]);
    fmt::print(fp, "{:10.4f} !p_trip2 Eq(20): Triple bond stabilisation\n", gp.l[7]);
    fmt::print(fp, "{:10.4f} !p_ovun7 Eq(12): Undercoordination\n", gp.l[8]);
    fmt::print(fp, "{:10.4f} !p_ovun8 Eq(12): Undercoordination\n", gp.l[9]);
    fmt::print(fp, "{:10.4f} !p_trip1 Eq(20): Triple bond stabilization\n", gp.l[10]);
    fmt::print(fp, "{:10.4f} !Lower Taper-radius (must be 0)\n", gp.l[11]);
    fmt::print(fp, "{:10.4f} !R_cut Eq(21): Upper Taper-radius\n", gp.l[12]);
    fmt::print(fp, "{:10.4f} !p_fe1 Eq(6a): Fe dimer correction\n", gp.l[13]);
    fmt::print(fp, "{:10.4f} !p_val6 Eq(13c): Valency undercoordination\n", gp.l[14]);
    fmt::print(fp, "{:10.4f} !p_lp1 Eq(8): Lone pair param\n", gp.l[15]);
    fmt::print(fp, "{:10.4f} !p_val9 Eq(13f): Valency angle exponent\n", gp.l[16]);
    fmt::print(fp, "{:10.4f} !p_val10 Eq(13g): Valency angle parameter\n", gp.l[17]);
    fmt::print(fp, "{:10.4f} !p_fe2 Eq(6a): Fe dimer correction\n", gp.l[18]);
    fmt::print(fp, "{:10.4f} !p_pen2 Eq(14a): Double bond/angle param\n", gp.l[19]);
    fmt::print(fp, "{:10.4f} !p_pen3 Eq(14a): Double bond/angle param\n", gp.l[20]);
    fmt::print(fp, "{:10.4f} !p_pen4 Eq(14a): Double bond/angle param\n", gp.l[21]);
    fmt::print(fp, "{:10.4f} !p_fe3 Eq(6a): Fe dimer correction\n", gp.l[22]);
    fmt::print(fp, "{:10.4f} !p_tor2 Eq(16b): Torsion/BO parameter\n", gp.l[23]);
    fmt::print(fp, "{:10.4f} !p_tor3 Eq(16c): Torsion overcoordination\n", gp.l[24]);
    fmt::print(fp, "{:10.4f} !p_tor4 Eq(16c): Torsion overcoordination\n", gp.l[25]);
    fmt::print(fp, "{:10.4f} !p_elho Eq(26a): electron-hole interaction\n", gp.l[26]);
    fmt::print(fp, "{:10.4f} !p_cot2 Eq(17b): Conjugation if tors13=0\n", gp.l[27]);
    fmt::print(fp, "{:10.4f} !p_vdW1 Eq(23b): vdWaals shielding\n", gp.l[28]);
    fmt::print(fp, "{:10.4f} !Cutoff for bond order (*100)\n", gp.l[29]);
    fmt::print(fp, "{:10.4f} !p_coa4 Eq(15): Valency angle conjugation\n", gp.l[30]);
    fmt::print(fp, "{:10.4f} !p_ovun4 Eq(11b): Over/Undercoordination\n", gp.l[31]);
    fmt::print(fp, "{:10.4f} !p_ovun3 Eq(11b): Over/Undercoordination\n", gp.l[32]);
    fmt::print(fp, "{:10.4f} !p_val8 Eq(13d): Valency/lone pair param\n", gp.l[33]);
    fmt::print(fp, "{:10.4f} !X_soft Eq(25): ACKS2 softness for X_ij\n", gp.l[34]);
    fmt::print(fp, "{:10.4f} !d Eq(23d): Scale factor in lg-dispersion\n", gp.l[35]);
    fmt::print(fp, "{:10.4f} !p_val Eq(27): Gauss exponent for electrons\n", gp.l[36]);
    fmt::print(fp, "{:10.4f} !1 Eq(13e): disable undecoord in val angle\n", gp.l[37]);
    fmt::print(fp, "{:10.4f} !p_coa3 Eq(15): Valency angle conjugation\n", gp.l[38]);

    //fmt::print(fp, "{:10.4f} !1: triple bond stabilization for all bonds
    //fmt::print(fp, "{:10.4f} !R_cut,e Eq(26): eReax taper radius
    //fmt::print(fp, "{:10.4f} !1: electron correction on valence angle

    // next line is number of atom types followed by 3 lines of comments
    const int ntypes = reax->num_atom_types;
    fmt::print(fp, "{:<4}   ! Nr of atoms; cov.r; valency;a.m;Rvdw;Evdw;gammaEEM;cov.r2;#\n"
                   "         alfa;gammavdW;valency;Eunder;Eover;chiEEM;etaEEM;n.u.\n"
                   "         cov r3;Elp;Heat inc.;bo131;bo132;bo133;softcut;n.u.\n"
                   "         ov/un;val1;n.u.;val3,vval4\n", ntypes);

    // ffield data
    auto &sbp = reax->sbp;
    auto &tbp = reax->tbp;
    auto &thbp = reax->thbp;
    auto &hbp = reax->hbp;
    auto &fbp = reax->fbp;

    // atomic parameters
    // four lines per atom type, or 5 if lgvdw != 0
    // the first starts with the symbol and has 9 words
    // the next three have 8 words
    // the fifth will have 2 words, if present

    const int lgflag = control->lgflag;

     for (i = 0; i < ntypes; ++i) {

      // line one

          // FIXME
          //if ((values.count() < 8) && !lgflag)
          //  THROW_ERROR("This force field file requires using 'lgvdw yes'");

      fmt::print(fp, " {:2} {:8.4f} {:8.4f} {:8.4f} {:8.4f} {:8.4f} {:8.4f} {:8.4f} {:8.4f}\n",
          sbp[i].name, sbp[i].r_s, sbp[i].valency, sbp[i].mass, sbp[i].r_vdw, sbp[i].epsilon, sbp[i].gamma, sbp[i].r_pi, sbp[i].valency_e );

      // line two

      fmt::print(fp, "    {:8.4f} {:8.4f} {:8.4f} {:8.4f} {:8.4f} {:8.4f} {:8.4f} {:8.4f}\n",
          sbp[i].alpha, sbp[i].gamma_w, sbp[i].valency_boc, sbp[i].p_ovun5, 0.0, sbp[i].chi, sbp[i].eta/2.0, (double)sbp[i].p_hbond );

      // line three

      fmt::print(fp, "    {:8.4f} {:8.4f} {:8.4f} {:8.4f} {:8.4f} {:8.4f} {:8.4f} {:8.4f}\n",
          sbp[i].r_pi_pi, sbp[i].p_lp2, 0.0, sbp[i].b_o_131, sbp[i].b_o_132, sbp[i].b_o_133, sbp[i].bcut_acks2, 0.0 );

      // line four

      fmt::print(fp, "    {:8.4f} {:8.4f} {:8.4f} {:8.4f} {:8.4f} {:8.4f} {:8.4f} {:8.4f}\n",
          sbp[i].p_ovun2, sbp[i].p_val3, 0.0, sbp[i].valency_val, sbp[i].p_val5, sbp[i].rcore2, sbp[i].ecore2, sbp[i].acore2 );

          // FIXME
          // read line five only when lgflag != 0

          //if (lgflag) {
          //  values = reader.next_values(0);
          //  ++lineno;
          //  CHECK_COLUMNS(2);
          //  sbp[i].lgcij    = values.next_double();
          //  sbp[i].lgre     = values.next_double();
          //} else sbp[i].lgcij = sbp[i].lgre = 0.0;

    }

    // next line is number of two body parameters followed by 1 comment line

    int nbondtypes = 0;
    for (j = 0; j < ntypes; ++j)
      for (k = 0; k <= j; ++k)
        if ( tbp[j][k].p_bo2 > 0.0 || tbp[j][k].p_bo4 > 0.0 || tbp[j][k].p_bo6 > 0.0 )
          nbondtypes++;

    fmt::print(fp, " {:<3}   ! Nr of bonds; Edis1;LPpen;n.u.;pbe1;pbo5;13corr;pbo6\n"
                   "                          pbe2;pbo3;pbo4;n.u.;pbo1;pbo2;ovcorr\n", nbondtypes);

    for (j = 0; j < ntypes; ++j)
      for (k = 0; k <= j; ++k)
        if ( tbp[j][k].p_bo2 > 0.0 || tbp[j][k].p_bo4 > 0.0 || tbp[j][k].p_bo6 > 0.0 )
        {

          // first line

          fmt::print(fp, " {:2}  {:2} {:8.4f} {:8.4f} {:8.4f} {:8.4f} {:8.4f} {:8.4f} {:8.4f} {:8.4f}\n",  k+1, j+1, tbp[k][j].De_s, tbp[k][j].De_p, tbp[k][j].De_pp, tbp[k][j].p_be1, tbp[k][j].p_bo5, tbp[k][j].v13cor, tbp[k][j].p_bo6, tbp[k][j].p_ovun1 );

          // second line

          fmt::print(fp, "        {:8.4f} {:8.4f} {:8.4f} {:8.4f} {:8.4f} {:8.4f} {:8.4f} {:8.4f}\n", tbp[k][j].p_be2, tbp[k][j].p_bo3, tbp[k][j].p_bo4, 0.0, tbp[k][j].p_bo1, tbp[k][j].p_bo2, tbp[k][j].ovc, 0.0 );

    }

/*
        // next line is number of two body off-diagonal parameters

        values = reader.next_values(0);
        n = values.next_int();
        ++lineno;

        double val;
        for (i = 0; i < n; ++i) {
          values = reader.next_values(0);
          ++lineno;
          CHECK_COLUMNS(8 + lgflag);

          j = values.next_int() - 1;
          k = values.next_int() - 1;

          if ((j < 0) || (k < 0))
            THROW_ERROR("Inconsistent force field file");

          if ((j < ntypes) && (k < ntypes)) {
            val = values.next_double();
            if (val > 0.0) tbp[j][k].D = tbp[k][j].D = val;

            val = values.next_double();
            if (val > 0.0) tbp[j][k].r_vdW = tbp[k][j].r_vdW = 2*val;

            val = values.next_double();
            if (val > 0.0) tbp[j][k].alpha = tbp[k][j].alpha = val;

            val = values.next_double();
            if (val > 0.0) tbp[j][k].r_s = tbp[k][j].r_s = val;

            val = values.next_double();
            if (val > 0.0) tbp[j][k].r_p = tbp[k][j].r_p = val;

            val = values.next_double();
            if (val > 0.0) tbp[j][k].r_pp = tbp[k][j].r_pp = val;

            if (lgflag) {
              val = values.next_double();
              if (val >= 0.0) tbp[j][k].lgcij = tbp[k][j].lgcij = val;
            }
          }
        }

        // next line is number of three body parameters

        values = reader.next_values(0);
        n = values.next_int();
        ++lineno;

        int cnt;
        for (i = 0; i < n; ++i) {
          values = reader.next_values(0);
          ++lineno;
          CHECK_COLUMNS(10);

          j = values.next_int() - 1;
          k = values.next_int() - 1;
          l = values.next_int() - 1;

          if ((j < 0) || (k < 0) || (l < 0))
            THROW_ERROR("Inconsistent force field file");

          if ((j < ntypes) && (k < ntypes) && (l < ntypes)) {

            cnt = thbp[j][k][l].cnt;
            thbp[j][k][l].cnt++;
            thbp[l][k][j].cnt++;

            val = values.next_double();
            thbp[j][k][l].prm[cnt].theta_00 = val;
            thbp[l][k][j].prm[cnt].theta_00 = val;

            val = values.next_double();
            thbp[j][k][l].prm[cnt].p_val1 = val;
            thbp[l][k][j].prm[cnt].p_val1 = val;

            val = values.next_double();
            thbp[j][k][l].prm[cnt].p_val2 = val;
            thbp[l][k][j].prm[cnt].p_val2 = val;

            val = values.next_double();
            thbp[j][k][l].prm[cnt].p_coa1 = val;
            thbp[l][k][j].prm[cnt].p_coa1 = val;

            val = values.next_double();
            thbp[j][k][l].prm[cnt].p_val7 = val;
            thbp[l][k][j].prm[cnt].p_val7 = val;

            val = values.next_double();
            thbp[j][k][l].prm[cnt].p_pen1 = val;
            thbp[l][k][j].prm[cnt].p_pen1 = val;

            val = values.next_double();
            thbp[j][k][l].prm[cnt].p_val4 = val;
            thbp[l][k][j].prm[cnt].p_val4 = val;
          }
        }

        // next line is number of four body parameters

        values = reader.next_values(0);
        n = values.next_int();
        ++lineno;

        for (i = 0; i < n; ++i) {
          values = reader.next_values(0);
          ++lineno;
          CHECK_COLUMNS(9);

          j = values.next_int() - 1;
          k = values.next_int() - 1;
          l = values.next_int() - 1;
          m = values.next_int() - 1;

          if ((j < -1) || (k < 0) || (l < 0) || (m < -1))
            THROW_ERROR("Inconsistent force field file");

          const double val1 = values.next_double();
          const double val2 = values.next_double();
          const double val3 = values.next_double();
          const double val4 = values.next_double();
          const double val5 = values.next_double();

          if ((j >= 0) && (m >= 0)) { // this means the entry is not in compact form

            if ((j < ntypes) && (k < ntypes) && (l < ntypes) && (m < ntypes)) {

              tor_flag[j][k][l][m] = tor_flag[m][l][k][j] = 1;
              fbp[j][k][l][m].cnt  = fbp[m][l][k][j].cnt  = 1;

              fbp[j][k][l][m].prm[0].V1     = fbp[m][l][k][j].prm[0].V1     = val1;
              fbp[j][k][l][m].prm[0].V2     = fbp[m][l][k][j].prm[0].V2     = val2;
              fbp[j][k][l][m].prm[0].V3     = fbp[m][l][k][j].prm[0].V3     = val3;
              fbp[j][k][l][m].prm[0].p_tor1 = fbp[m][l][k][j].prm[0].p_tor1 = val4;
              fbp[j][k][l][m].prm[0].p_cot1 = fbp[m][l][k][j].prm[0].p_cot1 = val5;
            }

          } else { // This means the entry is of the form 0-X-Y-0

            if ((k < ntypes) && (l < ntypes)) {
              for (int p = 0; p < ntypes; ++p) {
                for (int o = 0; o < ntypes; ++o) {

                  if (tor_flag[p][k][l][o] == 0) {
                    fbp[p][k][l][o].cnt = 1;
                    fbp[p][k][l][o].prm[0].V1 = val1;
                    fbp[p][k][l][o].prm[0].V2 = val2;
                    fbp[p][k][l][o].prm[0].V3 = val3;
                    fbp[p][k][l][o].prm[0].p_tor1 = val4;
                    fbp[p][k][l][o].prm[0].p_cot1 = val5;
                  }

                  if (tor_flag[o][l][k][p] == 0) {
                    fbp[o][l][k][p].cnt = 1;
                    fbp[o][l][k][p].prm[0].V1 = val1;
                    fbp[o][l][k][p].prm[0].V2 = val2;
                    fbp[o][l][k][p].prm[0].V3 = val3;
                    fbp[o][l][k][p].prm[0].p_tor1 = val4;
                    fbp[o][l][k][p].prm[0].p_cot1 = val5;
                  }
                }
              }
            }
          }
        }

        // next line is number of hydrogen bond parameters. that block may be missing

        for (i = 0; i < ntypes; ++i)
          for (j = 0; j < ntypes; ++j)
            for (k = 0; k < ntypes; ++k)
              hbp[i][j][k].r0_hb = -1.0;

        auto thisline = reader.next_line();
        if (!thisline) throw EOFException("ReaxFF parameter file has no hydrogen bond parameters");

        values = ValueTokenizer(thisline);
        n = values.next_int();
        ++lineno;

        for (i = 0; i < n; ++i) {
          values = reader.next_values(0);
          ++lineno;
          CHECK_COLUMNS(7);

          j = values.next_int() - 1;
          k = values.next_int() - 1;
          l = values.next_int() - 1;

          if ((j < 0) || (k < 0) || (l < 0))
            THROW_ERROR("Inconsistent force field file");

          if ((j < ntypes) && (k < ntypes) && (l < ntypes)) {
            hbp[j][k][l].r0_hb = values.next_double();
            hbp[j][k][l].p_hb1 = values.next_double();
            hbp[j][k][l].p_hb2 = values.next_double();
            hbp[j][k][l].p_hb3 = values.next_double();
          }
        }

        memory->destroy(tor_flag);
*/
      } catch (EOFException &e) {
        error->warning(FLERR, e.what());
      } catch (std::exception &e) {
        error->one(FLERR,e.what());
      }


  fclose(fp);

/*
  // apply global parameters to global control settings
  control->bo_cut    = 0.01 * reax->gp.l[29];
  control->nonb_low  = reax->gp.l[11];
  control->nonb_cut  = reax->gp.l[12];
*/

}
#undef THROW_ERROR

}
