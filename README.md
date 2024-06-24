# CoOMBE
**Co**mbined **O**ptical and **M**axwell-**B**loch **E**quations

This repository contains the following items:

## Documentation

- *doc/user_manual.pdf* : The User Manual, giving a detailed description
                  of the present codes. This document also contains
                  a short tutorial on writing the Hamiltonian within
                  the rotating wave approximation. 

## Program files 

source/driveall.f90         : The driveall program described in Section 4.4 of
                       the article and Section 4 of the User Manual.
source/general_settings.f90 : The general_settings module.
source/ldbl.f90             : The ldbl and ldblstore modules and the fcn_dummy and
                       solout_dummy subroutines.
source/mbe.f90              : The mbe module.
source/obe.f90              : The obe and obe_constants modules and the ext_setsys
                       subroutine.


## Examples

examples/ : A directory containing copies of the codes described in the
           appendices E and F of the article, as well as additional examples and associated documentation.
