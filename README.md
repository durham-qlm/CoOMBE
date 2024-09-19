# CoOMBE
**Co**mbined **O**ptical and **M**axwell-**B**loch **E**quations

CoOMBE is a collection of Fortran programs with two general aims:

1. integrating the optical Bloch equations for atoms interacting with one or several laser fields (i.e., integrating the equations governing the time evolution of the density matrix representing the quantum state of these atoms),
2. integrating the Maxwell-Bloch equations for one or two laser fields (co)- propagating in an atomic vapour (i.e., integrating the equations describing how these fields evolve when propagating in such a medium).

This cillection of programs has been published open access in Computer Physics Communications: 

R M Potvliege and S A Wrathmall, CoOMBE: A suite of open-source programs for the integration of the optical Bloch equations and Maxwell-Bloch equations, Comput. Phys. Commun. 306 (2025) 109374, https://doi.org/10.1016/j.cpc.2024.109374

This folder contains the following items:

1. Documentation
----------------

user_manual.pdf : The User Manual, giving a detailed description
                  of the present codes. This document also contains
                  a short tutorial on writing the Hamiltonian within
                  the rotating wave approximation. 

2. Program files 
----------------

driveall.f90         : The driveall program described in Section 4.4 of
                       the article and Section 4 of the User Manual.
general_settings.f90 : The general_settings module.
ldbl.f90             : The ldbl and ldblstore modules and the fcn_dummy and
                       solout_dummy subroutines.
mbe.f90              : The mbe module.
obe.f90              : The obe and obe_constants modules and the ext_setsys
                       subroutine.


3. Examples
-----------

examples : A directory containing copies of the codes described in the
           appendices E and F of the article, as well as additional examples
           and associated documentation.
