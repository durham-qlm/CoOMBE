# CoOMBE
**Co**mbined **O**ptical and **M**axwell-**B**loch **E**quations

CoOMBE is a collection of Fortran programs with two general aims:

1. integrating the optical Bloch equations for atoms interacting with one or several laser fields (i.e., integrating the equations governing the time evolution of the density matrix representing the quantum state of these atoms),
2. integrating the Maxwell-Bloch equations for one or two laser fields (co)- propagating in an atomic vapour (i.e., integrating the equations describing how these fields evolve when propagating in such a medium).

## Repository contents

This repository contains the following items:

### Documentation

- **doc/user_manual.pdf** : The User Manual, giving a detailed description
                  of the present codes. This document also contains
                  a short tutorial on writing the Hamiltonian within
                  the rotating wave approximation. 

### Program files 

- **source/driveall.f90** : The driveall program described in Section 4.4 of the article and Section 4 of the User Manual.
- **source/general_settings.f90** : The general_settings module.
- **source/ldbl.f90** : The ldbl and ldblstore modules and the fcn_dummy and solout_dummy subroutines.
- **source/mbe.f90**  : The mbe module.
- **source/obe.f90** : The obe and obe_constants modules and the ext_setsys subroutine.


### Examples

**examples/** : A directory containing copies of the codes described in the appendices E and F of the article, as well as additional examples and associated documentation.

## Prerequisites

The standard installation requires a Fortran 90 compiler. The LAPACK and BLAS numerical libraries also need to be installed. We also provide an alternative method of running the programs using podman (REF). This require no installation of a compiler or numerical librariers on the user's machine.

## License

The files distributed with this program are provided subject to the GNU General Public License, Version 3.0. A Copy of the license is provided.

## Change Log
v 1.0 (XX/XX/2024)
- Program uploaded to repoitory.
