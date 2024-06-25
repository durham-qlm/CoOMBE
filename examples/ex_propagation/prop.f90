      program prop
 
!  Modules directly used by this program:
      use general_settings
      use obe
      use mbe
 
!  Declare all the variables. The type obecfield is defined in
!  the obe module. The variable nst (the number of states) is
!  defined in the general_settings module.
      implicit none
      type(obecfield) :: coupling_field,probe_field
      double precision, dimension(nst,nst) :: Gamma_decay_f
      double precision, dimension(nst*nst) :: rhovec
      double precision, dimension(nst) :: energ_f, popinit
      double complex :: Field_c_0, Field_p_0
      double precision :: density,tmax,tmin,tw,t0,zmax
      integer :: iDoppler,iinterp,imethod,istart,iRabi,iweakprb,izrule,   &
                 nfields,nsubsteps,n_time_steps,n_z_steps,nt_writeout,    &
                 nz_writeout
      logical :: force0
      external :: my_output_pr

 
!  Properties of the probe field
      probe_field%wavelength = 794.979d0
      probe_field%idir = 1 ! The field propagates in the positive direction.
      probe_field%detuning = 0.0d0
      probe_field%detuning_fact(1) =  0.0d0
      probe_field%detuning_fact(2) = -1.0d0
      probe_field%detuning_fact(3) =  0.0d0
      probe_field%dip_mom = 0.0d0
      probe_field%dip_mom(2,1) = 1.465d-29
 
!  Properties of the coupling field
      coupling_field%wavelength = 780.241d0
      coupling_field%idir = 1 ! The field propagates in the positive direction.
      coupling_field%detuning = 0.0d0
      coupling_field%detuning_fact(1) =  0.0d0
      coupling_field%detuning_fact(2) =  0.0d0
      coupling_field%detuning_fact(3) = -1.0d0
      coupling_field%dip_mom = 0.0d0
      coupling_field%dip_mom(3,1) = 2.06937d-29
 
!  Frequency offset of each of the states. Here they are zero
!  as the three states are assumed to coincide in energy with
!  their respective reference energy level.
      energ_f(1) = 0.0d0
      energ_f(2) = 0.0d0
      energ_f(3) = 0.0d0
 
!  States 2 and 3 both decay to state 1. Corresponding decay rates:
      Gamma_decay_f = 0.0d0
      Gamma_decay_f(1,2) = 5.746d0 
      Gamma_decay_f(1,3) = 6.0666d0

!  Parameters of the Gaussian pulse for the coupling field:
      tw = 1.131d-03  ! 0.008 times sqrt(2), for 0.8 ns FWHM in power
      t0 = 0.d0
      force0 = .false. ! No forcing of the field to zero at the initial time.
!
!  The electric field amplitudes of the two fields
      Field_p_0 = (8.68021d0,0.d0)  ! 10 muW/cm^2 constant intensity
      Field_c_0 = (86802.1d0,0.d0)  ! 1 kW/cm^2 peak intensity
!
!  Initial and final times (in mus), and number of time steps
      tmin = -0.002d0
      tmax =  0.002d0
      n_time_steps = 100
 
!  Initialisation
      nfields = 2    ! Number of fields
      iweakprb = 0   ! 0 means that the weak probe approximation
                     ! is not made
      iRabi = 0      ! 0 means that the Rabi frequencies are not provided.
      call obe_setcsts(energ_f,Gamma_decay_f,nfields,iweakprb,iRabi)
      call obe_setfields(1,probe_field)    ! Declare what is field 1 
      call obe_setfields(2,coupling_field) ! Declare what is field 2 

!  Define the fields
      call mbe_set_envlp('probe','cw')
      call mbe_set_envlp('coupl','Gs',tw=tw,t0=t0,force0=force0)
      iinterp = 1
      call mbe_set_tdfields_A(tmin,tmax,n_time_steps,nfields,   &
                              iinterp,Field_p_0,Field_c_0)
 
!  Prepare the propagation calculation.
      istart = 2
      popinit = 0.d0
      popinit(1) = 1.d0
      imethod = 5
      nsubsteps = 2
      izrule = 3
      zmax = 16.d0
      n_z_steps = 1600
      density = 1.96d21
      iDoppler = 0
      nz_writeout = 20
      nt_writeout = 1

!  Call mbe_propagate_2, asking this subroutine to use the external
!  subroutine my_output_pr for writing out the results.
      call mbe_propagate_2(istart,popinit,imethod,nsubsteps,     &
                           izrule,zmax,n_z_steps,probe_field,    &
                           coupling_field,density,iDoppler,      &
                           nz_writeout,nt_writeout,my_output_pr)
      
      end program prop

      subroutine my_output_pr(i,k,iv,z,tp,rhovec,F_p_re,F_p_im,   &
                              F_c_re,F_c_im)
!  Used by mbe_propagate_2 to write out the results of the calculation

      use general_settings
      real(kd), dimension(nst*nst), intent(in) :: rhovec
      real(kd), intent(in) :: F_p_re,F_p_im,F_c_re,F_c_im,tp,z
      integer, intent(in) :: i,iv,k

!  Write tp, z, and the real and imaginary parts of each amplitude:
      print 1000,tp,z,F_p_re,F_p_im,F_c_re,F_c_im
 1000 format(1x,6(1pe12.5,1x))

      end subroutine my_output_pr
