Reactive Force Field (ReaxFF)
=============================

The `Reactive Force Field (ReaxFF) <https://doi.org/10.1038/npjcompumats.2015.11>`_ replaces fixed bond topologies of classical force fields with the concept of bond order to simulate bond breaking/formation of chemical reactions. Originally conceived for hydrocarbons in the gas phase :footcite:p:`vanduin2001`, ReaxFF has been extended to a wide range of applications :footcite:p:`senftle2016`. ReaxFF in LAMMPS :footcite:p:`aktulga2012` supports both the QEq charge equilibration :footcite:p:`rappe1991,nakano1997` and atom-condensed Kohn-Sham DFT to second order (ACKS2) :footcite:p:`verstraelen2013,ohearn2020` methods to represent the dynamics of electron density while fixed partial charges in classical force fields (eg. CHARMM) do not.


ReaxFF Potential Functions
--------------------------

.. math::

  E_{system} & = E_{bond} + E_{lp} + E_{over} + E_{under} + E_{val} + E_{pen} + E_{coa} + E_{C2}\\[.6em]
  & \qquad + E_{triple} + E_{tors} + E_{conj} + E_{H-bond} + E_{vdWaals} + E_{Coulomb}


ReaxFF LAMMPS commands
----------------------

* :doc:`pair_style <pair_reaxff>` reaxff (/kk)
* :doc:`fix <fix_qeq_reaxff>` qeq/reaxff (/kk)
* :doc:`fix <fix_acks2_reaxff>` acks2/reaxff (/kk)
* :doc:`fix <fix_qtpie_reaxff>` qtpie/reaxff (/kk)
* :doc:`pair_style <fix_reaxff_bonds>` reaxff/bonds (/kk)
* :doc:`pair_style <fix_reaxff_species>` reaxff/species (/kk)
* :doc:`compute <compute_reaxff_atom>` reaxff/atom (/kk)


ReaxFF Force Fields
-------------------

.. list-table:: Historical serial Fortran 77 force fields (no longer compatible and not available)
   :widths: 10 10 10 70
   :header-rows: 1
   :align: center

   * - Branch
     - Atoms
     - Filename
     - Source
   * - combustion
     - C / H
     - *n/a*
     - :footcite:t:`vanduin2001`

Combustion Branch
^^^^^^^^^^^^^^^^^

.. list-table:: Available COMBUSTION force fields in LAMMPS
   :widths: 10 10 10 10 60
   :header-rows: 1
   :align: center

   * - Branch
     - Atoms
     - Filename (LAMMPS)
     - `Filename (SCM) <https://www.scm.com/doc/ReaxFF/Included_Forcefields.html>`_
     - Source
   * - combustion
     - H/O/N/B
     - reaxff-weismiller2010.ff
     - Ab.ff
     - :footcite:t:`weismiller2010`
   * - FIXME
     - H/O/Au
     - reaxff-joshi2010.ff
     - *n/a*
     - :footcite:t:`joshi2010`
   * - FIXME
     - C/H/O/N
     - reaxff-budzien2009.ff
     - *n/a*
     - :footcite:t:`budzien2009`
   * - combustion
     - C/H/O
     - reaxff-chenoweth2008a.ff
     - CHO.ff
     - :footcite:t:`chenoweth2008a`
   * - FIXME
     - C/H/O/N/S/F/Pt/Cl/Ni/X
     - reaxff-singh2013.ff
     - *n/a*
     - :footcite:t:`singh2013`
   * - FIXME
     - C/H/O/N/S
     - reaxff-mattsson2010.ff
     - *n/a*
     - :footcite:t:`mattsson2010`
   * - combustion
     - C/H/O/N
     - reaxff-strachan2003.ff
     - *n/a*
     - :footcite:t:`strachan2003`
   * - combustion
     - V/O/C/H
     - reaxff-chenoweth2008b.ff
     - VOCH.ff
     - :footcite:t:`chenoweth2008b`
   * - combustion
     - C/H/O/N/S/Si
     - reaxff-liu2011.ff
     - dispersion/CHONSSi-lg.ff
     - :footcite:t:`liu2011`
   * - combustion
     - C/H/Na
     - reaxff-hjertenaes2016.ff
     - CHNa.ff
     - :footcite:t:`hjertenaes2016`
   * - combustion
     - C/H/O/S
     - reaxff-muller2016.ff
     - Mue2016.ff
     - :footcite:t:`muller2016`
   * - combustion
     - C/H/O
     - reaxff-ashraf2017.ff
     - CHO-2016.ff
     - :footcite:t:`ashraf2017`
   * - combustion
     - C/H
     - reaxff-mao2017.ff
     - CH_aromatics.ff
     - :footcite:t:`mao2017`
   * - combustion
     - C/H/O/S
     - reaxff-komissarov2021.ff
     - *n/a*
     - :footcite:t:`komissarov2021`

Independent Branch
^^^^^^^^^^^^^^^^^^

.. list-table:: Available INDEPENDENT force fields in LAMMPS
   :widths: 10 10 10 10 60
   :header-rows: 1
   :align: center

   * - Branch
     - Atoms
     - Filename (LAMMPS)
     - `Filename (SCM) <https://www.scm.com/doc/ReaxFF/Included_Forcefields.html>`_
     - Source
   * - independent
     - Zr/Y/O/H
     - reaxff-mayernick2010.ff
     - ZrYOHVac.ff
     - :footcite:t:`mayernick2010`
   * - independent
     - C/H/Fe
     - reaxff-islam2016.ff
     - CHFe.ff
     - :footcite:t:`islam2016`
   * - independent
     - H/O/Si/Al/Li
     - reaxff-ostadhossein2016.ff
     - HOSiAlLi.ff
     - :footcite:t:`ostadhossein2016`
   * - independent
     - C/H/Ar/He/Ne/Kr
     - reaxff-yoon2016.ff
     - CHArHeNeKr.ff
     - :footcite:t:`yoon2016`
   * - independent
     - H/S/Mo
     - reaxff-ostadhossein2017.ff
     - HSMo.ff
     - :footcite:t:`ostadhossein2017`
   * - independent
     - Ru/N/H
     - reaxff-kim2018.ff
     - RuNH.ff
     - :footcite:t:`kim2018`
   * - independent
     - C/H/O/Ge
     - reaxff-nayir2018.ff
     - CHOGe.ff
     - :footcite:t:`nayir2018`
   * - independent
     - C/H/O/Li/Al/Ti/P
     - reaxff-shin2018.ff
     - CHOLiAlTiP.ff
     - :footcite:t:`shin2018`
   * - independent
     - C/H/O/N/S/Si
     - reaxff-soria2018.ff
     - CHONSSi.ff
     - :footcite:t:`soria2018`
   * - independent
     - C/H/O/S/Cu/Cl/X
     - reaxff-yeon2018.ff
     - CuSCH.ff
     - :footcite:t:`yeon2018`
   * - independent
     - Cu/Zr
     - reaxff-huang2019.ff
     - CuZr.ff
     - :footcite:t:`huang2019`
   * - independent
     - Zr/Y/O/Ni/H
     - reaxff-liu2019.ff
     - ZrYONiH.ff
     - :footcite:t:`liu2019`
   * - independent
     - Si/O/H
     - reaxff-nayir2019.ff
     - SiOHv2.ff
     - :footcite:t:`nayir2019`
   * - independent
     - C/H/O/N/S/Zr
     - reaxff-dwivedi2020.ff
     - CHONSZr.ff
     - :footcite:t:`dwivedi2020`
   * - independent
     - C/H/O/N/Si
     - reaxff-dwivedi2020.ff
     - CHONSi.ff
     - :footcite:t:`wang2020`
   * - independent
     - | C/H/Ga
       | C/H/In
     - | reaxff-rajabpour2021a.ff
       | reaxff-rajabpour2021b.ff
     - | GaCH-2020.ff
       | InCH-2020.ff
     - :footcite:t:`rajabpour2021`
   * - independent
     - Ni/Cr
     - reaxff-shin2021.ff
     - NiCr.ff
     - :footcite:t:`shin2021`
   * - independent
     - W/S/H/Al/O
     - reaxff-nayir2021.ff
     - WSHAlO.ff
     - :footcite:t:`nayir2021`
   * - independent
     - Si/Al/Mg/O
     - reaxff-yeon2021.ff
     - SiAlMgO.ff
     - :footcite:t:`yeon2021`
   * - independent
     - I/Pb/Cs/X
     - reaxff-pols2021.ff
     - CsPbI.ff
     - :footcite:t:`pols2021`
   * - independent
     - H/O/N/Si/F
     - reaxff-kim2021.ff
     - HONSiF.ff
     - :footcite:t:`kim2021`
   * - independent
     - C/H/O/S/Mo/Ni/Au/Ti
     - reaxff-mao2022.ff
     - CHOSMoNiAuTi.ff
     - :footcite:t:`mao2022`
   * - independent
     - Ru/H
     - reaxff-onwudinanti2022.ff
     - RuH.ff
     - :footcite:t:`onwudinanti2022`
   * - independent
     - C/H/O/N/S/Mg/P/Na/Cu/Cl/Ti/X
     - reaxff-hou2022.ff
     - CHONSMgPNaCuClTi.ff
     - :footcite:t:`hou2022`
   * - independent
     - Mg/O
     - reaxff-fiesinger2023.ff
     - MgO.ff
     - :footcite:t:`fiesinger2023`
   * - independent
     - Li/Si/C
     - reaxff-olou2023.ff
     - LiSiC.ff
     - :footcite:t:`olou2023`
   * - independent
     - Ni/Al
     - reaxff-du2023.ff
     - NiAl.ff
     - :footcite:t:`du2023`
   * - independent
     - C/H/O/N/B/Al/Si/Cl
     - reaxff-uene2024.ff
     - CHONBAlSiCl.ff
     - :footcite:t:`uene2024`
   * - independent
     - C/H/O/N/S/Si/Ge/Ga/Ag
     - reaxff-niefind2024.ff
     - CHONSSiGeGaAg.ff
     - :footcite:t:`niefind2024`
   * - independent
     - I/Br/Pb/Cs
     - reaxff-pols2024.ff
     - IBrPbCs.ff
     - :footcite:t:`pols2024`




Water Branch
^^^^^^^^^^^^

.. list-table:: Available WATER force fields in LAMMPS
   :widths: 10 10 10 10 60
   :header-rows: 1
   :align: center

   * - Branch
     - Atoms
     - Filename (LAMMPS)
     - `Filename (SCM) <https://www.scm.com/doc/ReaxFF/Included_Forcefields.html>`_
     - Source
   * - water
     - C/H/O/Fe
     - reaxff-aryanpour2010.ff
     - FeOCHCl.ff
     - :footcite:t:`aryanpour2010`
   * - water
     - C/H/O/Zn
     - reaxff-han2010.ff
     - CHOZn.ff
     - :footcite:t:`han2010`
   * - water
     - Zn/O/H
     - reaxff-raymand2010.ff
     - ZnOH.ff
     - :footcite:t:`raymand2010`
   * - water
     - C/H/O/N
     - reaxff-rahaman2011.ff
     - Glycine.ff
     - :footcite:t:`rahaman2011`
   * - water
     - C/H/O/Ca/Si/X
     - reaxff-manzano2012.ff
     - CaSiOH.ff
     - :footcite:t:`manzano2012`
   * - water
     - C/H/O/N/S/Mg/P/Na/Cu
     - reaxff-huang2013.ff
     - CuBTC.ff
     - :footcite:t:`huang2013`
   * - water
     - C/H/O/N/S/Mg/P/Na/Cu/Cl
     - reaxff-monti2013.ff
     - CHONSMgPNaCuCl_v2.ff
     - :footcite:t:`monti2013`
   * - water
     - Al/C/H/O
     - reaxff-hong2016.ff
     - AlCHO.ff
     - :footcite:t:`hong2016`
   * - water
     - C/H/O/N/Si/Cu/Ag/Zn
     - reaxff-lloyd2016.ff
     - AgZnO.ff
     - :footcite:t:`lloyd2016`
   * - water
     - C/H/O/N/S/Mg/P/Na/Ti/Cl/F/Au
     - reaxff-monti2016.ff
     - CHONSMgPNaTiClFAu.ff
     - :footcite:t:`monti2016`
   * - water
     - C/H/O/Si/Na
     - reaxff-hahn2018.ff
     - CHOSiNa.ff
     - :footcite:t:`hahn2018`
   * - water
     - H/O/X
     - reaxff-zhang2017.ff
     - Water2017.ff
     - :footcite:t:`zhang2017`
   * - water
     - C/H/O/Al/Ge/X
     - reaxff-zheng2017.ff
     - CHOAlGeX.ff: ()
     - :footcite:t:`zheng2017`
   * - water
     - C/H/O/N
     - reaxff-trnka2018.ff
     - *n/a*
     - :footcite:t:`trnka2018`
   * - water
     - C/H/O/N/S/Mg/P/Na/Cu/Cl/X
     - reaxff-zhang2018.ff
     - CHON2017_weak.ff
     - :footcite:t:`zhang2018`
   * - water
     - C/H/O/Cs/K/Na/Cl/I/F/Li
     - reaxff-fedkin2019.ff
     - CHOCsKNaClIFLi.ff
     - :footcite:t:`fedkin2019`
   * - water
     - C/H/O/N
     - reaxff-kowalik2019.ff
     - CHON-2019.ff
     - :footcite:t:`kowalik2019`
   * - water
     - C/H/O/N/S/Mg/P/Na/Ti/Cl/F/K/Li
     - reaxff-ganeshan2020.ff
     - CHONSMgPNaTiClFKLi.ff
     - :footcite:t:`ganeshan2020`
   * - water
     - C/H/O/S/Mo/Ni/Li/B/F/P/N
     - reaxff-liu2021.ff
     - CHOSMoNiLiBFPN-2.ff
     - :footcite:t:`liu2021`









----


FitSNAP-ReaxFF
--------------

If a parameter set is not available for your intented application, then you can fit a new parameter set with `FitSNAP <https://fitsnap.github.io/>`_ from DFT training data.

.. table:: Parameters that can be optimized
  :widths: auto
  :align: center

  ===== ========= ====================================
  Block Parameter Description
  ===== ========= ====================================
  ATM   r_s       Sigma bond covalent radius
  ATM   r_pi      Pi bond covalent radius
  ATM   r_pi2     Double pi bond covalent radius
  BND   p_bo1     Sigma bond order
  BND   p_bo2     Sigma bond order
  BND   p_bo3     Pi bond order parameter
  BND   p_bo4     Pi bond order parameter
  BND   p_bo5     Double pi bond order parameter
  BND   p_bo6     Double pi bond order parameter
  BND   p_be1     Bond energy parameter
  BND   p_be2     Bond energy parameter
  BND   De_s      Sigma-bond dissociation energy
  BND   De_p      Pi-bond dissociation energy
  BND   De_pp     Double pi-bond dissociation energy
  BND   p_ovun1   Overcoordination penalty
  OFD   r_s       Sigma bond length
  OFD   r_pi      Pi bond length
  OFD   r_pi2     PiPi bond length
  ANG   theta_00  180o-(equilibrium angle)
  ANG   p_val1    Valence angle parameter
  ANG   p_val2    Valence angle parameter
  TOR   V1        V1-torsion barrier
  TOR   V2        V2-torsion barrier
  TOR   V3        V3-torsion barrier
  TOR   p_tor1    Torsion angle parameter
  HBD   r0_hb     Hydrogen bond equilibrium distance
  HBD   p_hb1     Hydrogen bond energy
  ===== ========= ====================================

----

ReaxFF Bibliography
-------------------

  :download:`download reaxff.bib<reaxff.bib>`

.. footbibliography::

