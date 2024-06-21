      program example
 
!  Modules directly used by this program:
      use general_settings
      use obe
 
!  Declare all the variables. The type obecfield is defined in
!  the obe module. The variable nst (the number of states) is
!  defined in the general_settings module.
      implicit none
      type(obecfield) :: coupling_field, probe_field
      double precision, dimension(nst,nst) :: Gamma_decay_f
      double precision, dimension(nst*nst) :: rhovec
      double precision, dimension(nst) :: energ_f
      integer :: ioption, iRabi, iweakprb, mim, mre, nfields
 
!  Properties of the probe field
      probe_field%detuning = 5.0d0
      probe_field%detuning_fact(1) =  0.0d0
      probe_field%detuning_fact(2) = -1.0d0
      probe_field%detuning_fact(3) = -1.0d0
      probe_field%Rabif = (0.0d0 , 0.0d0)
      probe_field%Rabif(1,2) = (5.0d0 , 0.0d0) 
 
!  Properties of the coupling field
      coupling_field%detuning = 0.0d0
      coupling_field%detuning_fact(1) =  0.0d0
      coupling_field%detuning_fact(2) =  0.0d0
      coupling_field%detuning_fact(3) = -1.0d0
      coupling_field%Rabif = (0.d0 , 0.0d0)
      coupling_field%Rabif(2,3) = (10.0d0 , 0.0d0)
 
!  Frequency offset of each of the states. Here they are zero
!  as the three states are assumed to coincide in energy with
!  their respective reference energy level.
      energ_f(1) = 0.0d0
      energ_f(2) = 0.0d0
      energ_f(3) = 0.0d0
 
!  State 2 decays to state 1 and state 3 decays to state 2.
!  Corresponding decay rates:
      Gamma_decay_f = 0.0d0
      Gamma_decay_f(1,2) = 5.0d0 
      Gamma_decay_f(2,3) = 1.0d0
 
!  Initialisation
      nfields = 2    ! Number of fields
      iweakprb = 0   ! 0 means that the weak probe approximation
                     ! is not made
      iRabi = 1      ! 1 means that the Rabi frequencies are
                     ! provided! directly rather than through
                     ! dipole matrix elements and field amplitudes.
      call obe_setcsts(energ_f,Gamma_decay_f,nfields,iweakprb,iRabi)
      call obe_setfields(1,probe_field)    ! Declare what is field 1 
      call obe_setfields(2,coupling_field) ! Declare what is field 2 
 
!  Calculation of the steady state density matrix. The result is
!  returned by obe_steadystate trough the 1D array rhovec.
      ioption = 1    ! See the description of obe_steadystate for
                     ! that option.
      call obe_steadystate(rhovec,ioption)
 
!  The coherence rho(1,2) is printed out.
      call obe_coher_index(1,2,mre,mim)
!  Given the values of mre and mim returned by obe_coher_index,
!  rhovec(mre) and rhovec(mim) contain, respectively, Re rho_12 
!  and Im rho_12.
      print 1000,rhovec(mre),rhovec(mim)
 1000 format(1x,'rho(1,2) = ',2(1pe12.5,2x))
 
      end program example
