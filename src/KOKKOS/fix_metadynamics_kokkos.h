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
FixStyle(metadynamics/kk,FixMetadynamicsKokkos<LMPDeviceType>);
FixStyle(metadynamics/kk/device,FixMetadynamicsKokkos<LMPDeviceType>);
FixStyle(metadynamics/kk/host,FixMetadynamicsKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_METADYNAMICS_KOKKOS_H
#define LMP_FIX_METADYNAMICS_KOKKOS_H

#include "fix_metadynamics.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

struct TagFixMetadynamics{};

template<class DeviceType>
class FixMetadynamicsKokkos : public FixMetadynamics {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  FixMetadynamicsKokkos(class LAMMPS *, int, char **);
  ~FixMetadynamicsKokkos() override;
  void init() override;
  void post_force(int) override;

  KOKKOS_INLINE_FUNCTION
  double rmsd(double *);

  KOKKOS_INLINE_FUNCTION
  double gpu_q_j(typename AT::t_x_array, typename AT::t_x_array, double*);

  KOKKOS_INLINE_FUNCTION
  double rmsd_grad_gpu(double*);

 private:
  typename AT::t_x_array d_x_group;

  typename AT::tdual_tagint_1d k_group_taglist;
  typename AT::t_tagint_1d_randomread d_group_taglist;

  typename AT::tdual_x_array k_ref_positions;
  typename AT::t_x_array d_ref_positions;

};

}

#endif
#endif

