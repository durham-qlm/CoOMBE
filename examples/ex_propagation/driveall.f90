      program driveall
!
!  A program based on the obe and mbe modules for integrating
!  the optical Bloch equations or the Maxwell-Bloch equations.
!
      use general_settings
      use obe
      use obe_constants, only : hbar
      use mbe
!
      implicit none
!
      integer, parameter :: nmx = nmn + nst - 1
!
!  Key parameters (read from the keyparams file)
!
      integer :: icmplxfld,nfields,nmin,nstates
      character(len=50) :: filename_controlparams,filename_defaultdata
      namelist/keyparams/filename_controlparams,filename_defaultdata,      &
                  icmplxfld,nfields,nmin,nstates
!
!  Other control parameters, atomic data, field data (read from the
!  controparams file and/or the defaultdata file)
!
      double complex, allocatable, dimension(:,:,:) :: cdip_mom,cRabif
      double complex, allocatable, dimension(:) :: camplitude
      double precision, allocatable, dimension(:,:,:) :: dip_mom,Rabif
      double precision, allocatable, dimension(:,:) :: detuning_fact
      double precision, allocatable, dimension(:) :: amplitude,detuning,   &
                                          detuningsarr,t_mesh,wavelength
      double precision, dimension(nmn:nmx,nmn:nmx) :: Gamma_decay_f,       &
                                                      add_dephas
      double precision, dimension(nmn:nmx) :: energ_f,popinit
      double precision, dimension(2) :: tw,t0,t1
      double precision :: atol,density,detuning_min,detuning_max,          &
                          detuning_step,rtol,ti,tf,urms,vmax,zmax
      integer, allocatable, dimension(:) :: idir
      integer, dimension(2) :: iforce0
      integer :: iAorB,iappend,icalc,iDoppler,iDoppler_numer_st,iinterp,   &
                 iladder_wkprb,imethod,index_field,inoncw,ioption,iRabi,   &
                 irate,iprintrho,irule,istart,isuscept,itdfieldsAorB,      &
                 iunformatted,ivarydetuning,iweakprb,izrule,nsubsteps,     &
                 n_time_steps,nt_writeout,n_v_values,n_z_steps,nz_writeout
      character(len=50), allocatable, dimension(:,:) :: filename_rho_out
      character(len=50) :: filename_chi_out,filename_Dopplerquad,          &
                           filename_tdamps_in,filename_rhoall_out,         &
                           filename_tdamps_out
      character(len=2), dimension(2) :: pulse_type
      namelist/controlparams/add_dephas,amplitude,atol,camplitude,         &
                  cdip_mom,cRabif,density,detuning,detuning_fact,          &
                  detuning_min,detuning_max,detuning_step,dip_mom,energ_f, &
                  filename_chi_out,filename_Dopplerquad,filename_rho_out,  &
                  filename_tdamps_in,filename_rhoall_out,                  &
                  filename_tdamps_out,Gamma_decay_f,iAorB,iappend,icalc,   &
                  idir,iDoppler,iDoppler_numer_st,iforce0,iinterp,         &
                  iladder_wkprb,imethod,index_field,inoncw,ioption,        &
                  iprintrho,iRabi,irate,irule,istart,isuscept,             &
                  itdfieldsAorB,iunformatted,ivarydetuning,iweakprb,       &
                  izrule,nsubsteps,n_time_steps,nt_writeout,n_v_values,    &
                  n_z_steps,nz_writeout,popinit,pulse_type,Rabif,rtol,     &
                  ti,tf,tw,t0,t1,urms,vmax,wavelength,zmax
      namelist/defaultdata/add_dephas,amplitude,camplitude,cdip_mom,       &
                  cRabif,density,detuning,detuning_fact,dip_mom,energ_f,   &
                  Gamma_decay_f,idir,popinit,Rabif,urms,wavelength
!
!  Other variables
!
      type(obecfield), allocatable, dimension(:) :: cfields
      type(obefield), allocatable, dimension(:) :: fields
      real(kd), dimension(nst*nst) :: rhovec
      double complex, allocatable, dimension(:) :: cplfield_ampl_arr,      &
                                                   prbfield_ampl_arr
      double complex :: chi,cplfield_ampl,prbfield_ampl
      double precision, allocatable, dimension(:) :: tmesh,vmesh,vweight
      double precision :: alpha,cplim,cplre,pi,prbim,prbre,refr_index
      integer, allocatable, dimension(:,:), save :: iunits_arr
      integer :: i,imltpls,istat,istore,iunit,iunit_ampl,j,jRabi,k,m,mre,  &
                 mim,nd,ndmax,nmax
      integer :: iprep
      logical :: fl_exist,fl_no_chi_outfile,fl_no_rho_outfile,             &
                 fl_Doppler_numer,force0
      external :: extsub_outpr
!
      common /common_driveall/ iDoppler,iunformatted,iunit,iunit_ampl,     &
                               nfields
!
      pi = 4.d0*datan(1.d0)
!
!  Read the key parameters from the standard input file
!  and check they are suitable.
!
      call sub_keyparams
!
!  Prepare the array containing the parameters of the field(s)
!  and initialise what needs to be initialised at this point.
!
      call sub_preparrays
      call sub_initialise

!  Get the rest of the necessary information from the controlparams file
!  and optionally from the defaultdata file.
!
      call sub_otherparams
!
!  If the input files specify real fields (read into the variable fields
!  of type obefield) rather than complex fields (read into the variable
!  cfields of type obecfield), copy these data into the obecfield
!  variable cfields, converting their type as necessary. Only the content
!  of cfields is used thereafter.
!
      if(icmplxfld.eq.0)then
         do k = 1,nfields
            call obe_fieldtocfield(fields(k),cfields(k))
         enddo
      endif
!
!  Prepare the output files and the numerical quadrature over atomic
!  velocities as necessary, and initialise obe.
!
      call sub_prepareoutputfiles_rho
      call sub_prepareDoppler
      call sub_convertfromRabif
      call sub_initialiseobe
!
!  Calculate the density matrix and, optionally, the susceptibility,
!  and output the results as required.
!
      if(icalc.eq.1)then
!
!  Time-dependent calculation
!
         call sub_openoutputfiles_td
         call sub_calcrho_td
         call sub_outputrho_td
         call sub_closeoutputfiles_td
!
      elseif(icalc.eq.2)then
!
!  Steady-state calculation
!
         call sub_preparedetunings_steady
         call sub_preparechiout
         do nd = 1,ndmax 
            if(nd.eq.1)then
               iprep = 1
            else
               iprep = 0
            endif
            call sub_calcrho_steady
            call sub_outputrho_steady
            call sub_suscept
         enddo
!
      elseif(icalc.eq.3)then
!
!  Integrate the Maxwell-Bloch equations to propagate the field(s).
!
         inoncw = 1   !  The fields are always non-cw in propagation
                      !  calculations.
         call sub_openoutputfiles_td
         call sub_calcrho_mb
         call sub_closeoutputfiles_td
!
      else
!
         call printerror
         print*,'Wrong value of icalc. ',icalc
         stop
!
      endif
!
      contains
 
!!!!!!!!!!!!!!!!!!!

      subroutine sub_keyparams
!
!  Read the key parameters from the standard input file
!  and check they are suitable.
!
      implicit none
!
      filename_defaultdata = 'undefined'
      icmplxfld = 0
!
      read(5,keyparams)
      if(nstates.ne.nst)then
         call printerror
         print*,'The number of states is inconsistent with the number'
         print*,'set in the general_settings module.',nstates,nst
         stop
      endif
      if(nfields.le.0)then
         call printerror
         print*,'The number of fields is nonsensical.',nfields
         stop
      endif
      if(nmin.ne.nmn)then
         call printerror
         print*,'The starting value of the index labelling the states is'
         print*,'inconsistent with the value set in the general_settings'
         print*,'module.',nmin,nmn
         stop
      endif
      if(icmplxfld.lt.0 .or. icmplxfld.gt.1)then
         call printerror
         print*,'The variable icmplxfld has been given an illegal value'
         print*,'in the keyparams file. ',icmplxfld
         stop
      endif
!
      nmax = nmin + nstates - 1   ! Largest state index
!
      end subroutine sub_keyparams

!!!!!!!!!!!!!!!!!!!

      subroutine sub_preparrays
!
!  Prepare the array containing the parameters of the field(s).
!
      implicit none
!
      allocate(amplitude(1:nfields),camplitude(1:nfields),      &
               detuning(1:nfields),idir(1:nfields),             &
               wavelength(1:nfields),stat=istat)
      if(istat.gt.0)then
         call printerror
         print*,'Error 1 in sub_preparrays.'
         stop
      endif
      allocate(detuning_fact(nmin:nmax,1:nfields),stat=istat)
      if(istat.gt.0)then
         call printerror
         print*,'Error 2 in sub_preparrays.'
         stop
      endif
      allocate(cdip_mom(nmin:nmax,nmin:nmax,1:nfields),         &
               cRabif(nmin:nmax,nmin:nmax,1:nfields),           &
               dip_mom(nmin:nmax,nmin:nmax,1:nfields),          &
               Rabif(nmin:nmax,nmin:nmax,1:nfields),stat=istat)
      if(istat.gt.0)then
         call printerror
         print*,'Error 3 in sub_preparrays.'
         stop
      endif
      allocate(cfields(1:nfields),fields(1:nfields),stat=istat)
      if(istat.gt.0)then
         call printerror
         print*,'Error 4 in sub_preparrays.'
         stop
      endif
!
!  Prepare the array containing the names of the output files.
!
      allocate(filename_rho_out(nmin:nmax,nmin:nmax),stat=istat)
      if(istat.gt.0)then
         call printerror
         print*,'Error 5 in sub_preparrays.'
         stop
      endif
!
      end subroutine sub_preparrays

!!!!!!!!!!!!!!!!!!!

      subroutine sub_initialise
!
!  Initialise the variables which need to be initialised
!  at this stage.
!
      implicit none
!
!  Default values of variables forming part of the namelists
!  controlparams or defaultdata. Any of these default values is 
!  superseded by the value read from either the controlparams
!  file or the defaultdata file if a value of the corresponding
!  variable is specified in either file.
!
!  The following variables have no default values, or the default
!  value set below will make the execution stop or give wrong results
!  if a suitable value is not specified in the controlparams file:
!  atol, density, detuning_min, detuning_max, detuning_step,
!  the amplitude and wavelength components of the cfields and fields
!  variables, iAorB, icalc, imethod, irule, itdfieldsAorB,
!  n_time_steps, n_v_values, n_z_steps, popinit, pulse_type,
!  rtol, ti, tf, tw, t0, t1, urms, vmax, zmax.
!  
!  Components of the fields and cfields variables:
!
      detuning = 0.d0
      amplitude = 0.d0
      camplitude = (0.d0,0.d0)
      wavelength = -1.d0
      detuning_fact = 0.d0
      idir = 1
      dip_mom = 0.d0
      cdip_mom = (0.d0,0.d0)
      Rabif = 0.d0
      cRabif = (0.d0,0.d0)
!
!  Other input data:
!
      add_dephas = 0.d0
      atol = 0.d0
      energ_f = 0.d0
      filename_chi_out = 'undefined'
      filename_Dopplerquad = 'undefined'
      filename_rho_out = 'undefined'
      filename_tdamps_in = 'undefined'
      filename_rhoall_out = 'undefined'
      filename_tdamps_out = 'undefined'
      Gamma_decay_f = 0.d0
      iAorB = -1 
      iappend = 1 
      icalc = 0
      iDoppler = 0
      iDoppler_numer_st = 0
      iforce0 = 0
      iinterp = 1
      iladder_wkprb = 0
      imethod = -1
      index_field = 1
      inoncw = 0
      ioption = 0
      iprintrho = 1
      iRabi = -1
      irate = 0
      irule = -1
      istart = 0
      isuscept = 0
      itdfieldsAorB = 0
      iunformatted = 0
      ivarydetuning = 0
      iweakprb = 0
      izrule = 3
      nsubsteps = 1
      n_time_steps = 0
      nt_writeout = 1
      n_v_values = -1
      n_z_steps = 0
      nz_writeout = 1
      popinit = 0.d0
      rtol = 0.d0
      urms = -1.d0
      vmax = 1.1d99
!
      end subroutine sub_initialise

!!!!!!!!!!!!!!!!!!!

      subroutine sub_otherparams
!
!  Get the rest of the necessary information from the controlparams file,
!  and also from the defaultdata file if there is one (starting with the
!  latter).
!
      implicit none
!
!  Get the default atomic and field data from the defaultdata file
!  if a name for this file was specified in the keyparams file.
!  Data specified in the defaultdata file will be superseded
!  by values specified in the controlparams file for the same
!  quantities, if any are specified in the latter.
!
      if(filename_defaultdata.ne.'undefined')then
         inquire(file=filename_defaultdata,exist=fl_exist)
         if(.not.fl_exist)then
            call printerror
            print*,'The defaultdata file ',filename_defaultdata
            print*,'is not found.'
            stop
         endif
         open(1,file=filename_defaultdata,status='OLD')
         read(1,defaultdata)
         close(1)
         print*,'Source file for the default data: ',   &
                filename_defaultdata
      endif
!
!  Get the control parameters and the rest of the input data from the
!  controlparams file specified in the keyparams file.
!
      inquire(file=filename_controlparams,exist=fl_exist)
      if(.not.fl_exist)then
         call printerror
         print*,'The controlparams file ',filename_controlparams
         print*,'is not found.'
         stop
      endif
      open(1,file=filename_controlparams,status='OLD')
      read(1,controlparams)
      close(1)
      print*,'Source file for the control parameters: ',  &
             filename_controlparams
!
!  Check that the value of iprintrho is legal 
!
      if(iprintrho.ne.0 .and. iprintrho.ne.1 .and. icalc.ne.3)then
         if(.not.(icalc.eq.2 .and. ivarydetuning.ge.1))then
            call printerror
            print*,'The control parameter iprintrho should be either'
            print*,'either 0 or 1.',iprintrho
            stop
         endif
      endif
!
!
!  Check that the value of ivarydetuning is consistent with the
!  type of calculation required.
!
      if(icalc.ne.2 .and. ivarydetuning.ge.1)then
         call printerror
         print*,'Calculations for a range of detunings are'
         print*,'possible only for steady state calculations.'
         stop
      endif
!
!  Pack the field data into either the obefield array fields or
!  the cobefield array cfields for later use.
!
      if(icmplxfld.eq.1)then
         do k = 1,nfields
            cfields(k)%detuning = detuning(k)
            cfields(k)%amplitude = camplitude(k)
            cfields(k)%wavelength = wavelength(k)
            cfields(k)%detuning_fact = detuning_fact(:,k)
            cfields(k)%idir = idir(k)
            cfields(k)%dip_mom = cdip_mom(:,:,k)
            cfields(k)%Rabif = cRabif(:,:,k)
         enddo
      else
         do k = 1,nfields
            fields(k)%detuning = detuning(k)
            fields(k)%amplitude = amplitude(k)
            fields(k)%wavelength = wavelength(k)
            fields(k)%detuning_fact = detuning_fact(:,k)
            fields(k)%idir = idir(k)
            fields(k)%dip_mom = dip_mom(:,:,k)
            fields(k)%Rabif = Rabif(:,:,k)
         enddo
      endif
!
      end subroutine sub_otherparams

!!!!!!!!!!!!!!!!!!!

      subroutine sub_prepareoutputfiles_rho
!
!  Prepare the files used to output specific elements of the
!  density matrix, as necessary. These files replace any
!  existing ones of the same name if a value of 0 is specified for
!  the control variable iappend in the input file.
!
      implicit none
!
      if(icalc.eq.1 .and. iprintrho.eq.0)then
         print*,'The density matrix is not written to the standard ', &
                'output stream (iprintrho = 0).'
      elseif(icalc.eq.2)then
         if(ivarydetuning.eq.0 .and. iprintrho.eq.0)then
            print*,'The density matrix is not written to the standard ', &
                   'output stream (iprintrho = 0).'
         elseif(ivarydetuning.ge.1)then
            print*,'The density matrix is not written to the standard ', &
                   'output stream (ivarydetuning > 0).'
         endif
      endif
!
      fl_no_rho_outfile = .true.
!
      if(any(filename_rho_out(:,:).ne.'undefined'))then
!
         fl_no_rho_outfile = .false.
!
         if(iappend.eq.0)then
            do j = nmin,nmax
               do i = nmin,nmax
                  if(filename_rho_out(i,j).ne.'undefined')then
                     inquire(file=filename_rho_out(i,j),exist=fl_exist)
                     if(fl_exist)then
                        open(1,file=filename_rho_out(i,j),status='REPLACE')
                        close(1)
                     endif
                  endif
               enddo
            enddo
         elseif(iappend.ne.1)then
            call printerror
            print*,'The value of iappend specified in the input file'
            print*,'is illegal. ',iappend
            print*,iappend
            stop
         endif
!
         print*,'Selected elements of the density matrix will be written to'
         print*,'the files specified for this in ',filename_controlparams
!
      endif
!
      end subroutine sub_prepareoutputfiles_rho

!!!!!!!!!!!!!!!!!!!

      subroutine sub_convertfromRabif
!
!  Re-define the dipole couplings in terms of complex electric field
!  amplitudes and complex dipole moments, if these couplings were
!  given in terms of complex Rabi frequencies only and this would not
!  be suitable for the calculation to be done (i.e., because mbe_tdint_1
!  mbe_tdint_2 are to be used). It is assumed that the values of the 
!  dipole moments is irrelevant, as long as in combination with the 
!  electric field amplitude they give the correct Rabi frequencies. For
!  simplicity, it is assumed that all the electric field amplitudes are
!  100 V/m and calculate the dipole moments accordingly.
!
!  This conversion is done to go round the impossibility of using
!  certain subroutines of obe or mbe when iRabi = 1. The difficulty is
!  avoided in this program by calling obe_setcsts with the variable
!  jRabi, not iRabi, in its list of arguments, where jRabi is normally
!  equal to iRabi unless pseudo complex amplitudes and dipole moments 
!  producing the correct Rabi frequencies have been prepared by the
!  present subroutine, in which case jRabi is set to 0.
!
!  Note that the dipole moments generated by this procedure may
!  considerably differ from their true values, as, of course, the
!  electric field amplitude. However, multiplying the dipole moments by
!  the amplitude will generate the correct Rabi frequencies.
!
      implicit none
!
!  The conversion is necessary only if the option iRabi = 1 has been
!  selected and the calculation required is to integrate the optical Bloch
!  equations for a non-cw field:
!
      if(icalc.eq.1 .and. inoncw.eq.1 .and. iRabi.eq.1)then
         jRabi = 0
         do k = 1,nfields
            cfields(k)%amplitude = (100.d0,0.d0)
            do j = nmin,nmax
               do i = j+1,nmax
                  if(cfields(k)%Rabif(i,j).ne.(0.d0,0.d0))then
                        cfields(k)%dip_mom(i,j) =                            &
                           cfields(k)%Rabif(i,j)*2.d0*pi*1.d6*hbar/          &
                                                     cfields(k)%amplitude
                  elseif(cfields(k)%Rabif(j,i).ne.(0.d0,0.d0))then
                        cfields(k)%dip_mom(j,i) =                            &
                           cfields(k)%Rabif(j,i)*2.d0*pi*1.d6*hbar/          &
                                                     cfields(k)%amplitude
                  endif
               enddo
            enddo
         enddo
!
      else
!
         jRabi = iRabi
!
      endif
!
      end subroutine sub_convertfromRabif

!!!!!!!!!!!!!!!!!!!

      subroutine sub_initialiseobe
!
!  Initialise obe
!
      implicit none
!
      imltpls = 0
      call obe_setcsts(energ_f,Gamma_decay_f,nfields,iweakprb,jRabi, &
                       imltpls,add_dephas)
      do k = 1,nfields
         call obe_setfields(k,cfields(k))
      enddo
!
      end subroutine sub_initialiseobe

!!!!!!!!!!!!!!!!!!!

      subroutine sub_prepareDoppler
!
!  Check that the relevant parameters have been specified in the
!  controlparams file if inhomogeneous broadening needs to be taken
!  into account. Also, prepare a numerical quadrature over atomic
!  velocities if this is required.
!
      implicit none
      integer :: k
!
!  First, check whether inhomogenous broadening needs to be taken
!  into account. If not, return to the calling program without doing
!  anything else.
!
      if(iDoppler.ne.1)return
!
!  From this point onwards inhomogenous broadening is assumed to be
!  necessary. Check that the necessary parameters have been provided
!  (the default values of idir(k) are used if not provided):
      if(urms.lt.0.d0)then
         call printerror
         print*,'No value or an illegal value of urms is specified'
         print*,'in the input file although a Doppler averaging is'
         print*,'required.'
         stop
      endif
      do k = 1,nfields
         if(wavelength(k).lt.0.d0)then
            call printerror
            print*,'No value or an illegal value of wavelength(k)'
            print*,'is specified in the input file for field ',k
            print*,'although a Doppler averaging is required.'
            stop
         endif
      enddo
!
!  Find out whether a numerical quadrature is required.
      fl_Doppler_numer = .false.
!
      if(icalc.eq.1)fl_Doppler_numer=.true.
      if(icalc.eq.2 .and. iladder_wkprb.ge.1 .and.        &
                             iDoppler_numer_st.eq.1)then
         call printerror
         print*,'No numerical Doppler averaging is possible under'
         print*,'the options iladder_wkprb.eq.1, 2 or 3.'
         stop
      endif
      if(icalc.eq.2 .and. iladder_wkprb.eq.0 .and.        &
                          iDoppler_numer_st.eq.1)fl_Doppler_numer=.true.
      if(icalc.eq.3)fl_Doppler_numer=.true.
!           
!  Prepare the numerical integration if this is necessary.
      if(fl_Doppler_numer)then
!
         if(irule.lt.0)then
            call printerror
            print*,'No value or an illegal value of irule is specified'
            print*,'in the input file although a numerical Doppler'
            print*,'averaging is required.'
            print*,irule
            stop
         elseif(irule.eq.0)then
            if(filename_Dopplerquad.eq.'undefined')then
               call printerror
               print*,'The name of the file containing the quadrature'
               print*,'abscissas and weights for the Doppler'
               print*,'averaging is not specified in the input file.'
               stop
            endif
            if(n_v_values.le.0)then
               call printerror
               print*,'No value or a wrong value n_v_values is'
               print*,'specified in the input file. ',n_v_values
               stop
            endif
            if(vmax.gt.1.d99)then  ! 1.1d99 is the default value of vmax.
               call printerror
               print*,'No value or a wrong value vmax is'
               print*,'specified in the input file. ',vmax
               stop
            endif
            allocate(vmesh(1:n_v_values),vweight(1:n_v_values),stat=istat)
            if(istat.gt.0)then
               call printerror
               print*,'Error at the allocation of vmesh and vweight.'
               stop
            endif
            inquire(file=filename_Dopplerquad,exist=fl_exist)
            if(.not.fl_exist)then
               call printerror
               print*,'The Dopplerquad file ',filename_Dopplerquad
               print*,'is not found.'
               stop
            endif
            open(1,file=filename_Dopplerquad,status='OLD')
            do k=1,n_v_values
               read(1,*,end=999,err=999)vmesh(k),vweight(k)
            enddo
            close(1)
            print*,'Source file for the numerical quadrature: ', &
                    filename_Dopplerquad
         endif
!
         call obe_set_Doppler(urms,n_v_values,irule,vmax,vmesh,vweight)
!
      endif
!
      return
!
 999  continue
!
      call printerror
      print*,'Error or premature end-of-file when reading the file'
      print*,'containing the quadrature mesh and weights for the'
      print*,'numerical Doppler averaging.'
      stop
!
      end subroutine sub_prepareDoppler

!!!!!!!!!!!!!!!!!!!

      subroutine sub_openoutputfiles_td
!
!  Open the files used to output specific elements of the
!  density matrix, if necessary, and assign a unit number to each
!  of these files. The unit numbers used for this range start at 10
!  and go upwards up to 99, which makes it possible to open up
!  to 90 files in that way. Then pass these unit numbers to obe.
!
!  If necessary, also open the file used for writing out the field(s)
!  at each time step.
!
      implicit none
!
      if(any(filename_rho_out(:,:).ne.'undefined'))then
!
         if(allocated(iunits_arr))deallocate(iunits_arr)
         allocate(iunits_arr(nmin:nmax,nmin:nmax),stat=istat)
         if(istat.gt.0)then
            call printerror
            print*,'Error at the allocation of iunits_arr.'
            stop
         endif
         iunits_arr = -1
!
         m = 9
         do j = nmin,nmax
            do i = nmin,nmax
               if(filename_rho_out(i,j).ne.'undefined')then
                  m = m + 1
                  if(m.gt.99)then
                     call printerror
                     print*,'Not more than 90 output files can be opened.'
                     stop
                  endif
                  open(m,file=filename_rho_out(i,j),status='UNKNOWN',      &
                       position='APPEND')
                  iunits_arr(i,j) = m
               endif
            enddo
         enddo
!
         call obe_setoutputfiles(iunits_arr)
!
      endif
!
!  Also, if necessary, create (if not existing previously) and open the
!  files used to dump the whole density matrix and/or write out the
!  amplitude(s) of the field(s) at each time step.
!
      if(filename_rhoall_out.ne.'undefined')then
!
         if(any(filename_rho_out(:,:).ne.'undefined'))then
            call printerror
            print*,'Specifying a rhoall_out file is incompatible with'
            print*,'specifying rho_out files.'
            stop
         endif
         if(iunformatted.ne.0 .and. iunformatted.ne.1)then
            call printerror
            print*,'iunformatted has been given an illegal value in the'
            print*,'input file.'
            stop
         endif
!
         iunit = 2
         inquire(file=filename_rhoall_out,exist=fl_exist)
         if(fl_exist .and. iappend.eq.0)then
            if(iunformatted.eq.1)then
               open(iunit,file=filename_rhoall_out,status='REPLACE', &
                                               form='UNFORMATTED')
            else
               open(iunit,file=filename_rhoall_out,status='REPLACE')
            endif
         elseif(fl_exist .and. iappend.eq.1)then
            if(iunformatted.eq.1)then
               open(iunit,file=filename_rhoall_out,status='UNKNOWN', &
                               form='UNFORMATTED',position='APPEND')
            else
               open(iunit,file=filename_rhoall_out,status='UNKNOWN', &
                                                  position='APPEND')
            endif
         else
            if(iunformatted.eq.1)then
               open(iunit,file=filename_rhoall_out,status='NEW', &
                                           form='UNFORMATTED')
            else
               open(iunit,file=filename_rhoall_out,status='NEW')
            endif
         endif
!
         print*,'The whole density matrix will be written at each step'
         print*,'to ',filename_rhoall_out
!
      else
!
         if(any(filename_rho_out(:,:).ne.'undefined'))then
            iunit = -1
         else
            iunit = 0
         endif
!
      endif
!
      if(inoncw.eq.1 .and. filename_tdamps_out.ne.'undefined')then
!
         iunit_ampl = 3
         inquire(file=filename_tdamps_out,exist=fl_exist)
         if(fl_exist .and. iappend.eq.0)then
            open(iunit_ampl,file=filename_tdamps_out,status='REPLACE')
         elseif(fl_exist .and. iappend.eq.1)then
            open(iunit_ampl,file=filename_tdamps_out,status='UNKNOWN', &
                                                   position='APPEND')
         else
            open(iunit_ampl,file=filename_tdamps_out,status='NEW')
         endif
!
         print*,'The field amplitude(s) will be written at each step'
         print*,'to ',filename_tdamps_out
!
      endif
!
      end subroutine sub_openoutputfiles_td

!!!!!!!!!!!!!!!!!!!

      subroutine sub_calcrho_td
!
!  Integrate the optical Bloch equations for one or several fields, either
!  cw fields or fields with a time-dependent amplitude. First prepare the
!  calculation and check some of the parameters.
!
      implicit none
!
!  The tolerance parameters used by DOP853 are first passed to obe if
!  this solver is to be used:
      if(imethod.eq.4)call obe_set_tol_dop853(rtol,atol)
!
      if(inoncw.eq.0)then 
         call sub_calcrho_td_cw
      elseif(inoncw.eq.1)then 
         call sub_definetdfields
         call sub_calcrho_td_notcw
      else
         call printerror
         print*,'Illegal value of the control parameter inoncw. ',inoncw
         stop
      endif
 
      end subroutine sub_calcrho_td

!!!!!!!!!!!!!!!!!!!

      subroutine sub_calcrho_td_cw
!
!  Integrate the optical Bloch equations for one or several cw fields.
!
      implicit none
!
      if(iDoppler.eq.1)then
         if(iAorB.eq.1)then
            call obe_Doppler_av_td_A(irate,ti,tf,n_time_steps,  &
                           imethod,rhovec,iunit,popinit)
         elseif(iAorB.eq.2)then
            if(irate.eq.1)then
               call printerror
               print*,'The option irate = 1 is incompatible with using'
               print*,'obe_Doppler_av_td_B (iAorB = 2).'
               stop
            endif
            if(imethod.ne.-1 .and. imethod.ne.3)then
               call printerror
               print*,'Requesting any other method than that corresponding'
               print*,'to imethod = 3 is incompatible with using'
               print*,'obe_Doppler_av_td_B (iAorB = 2).'
               stop
            endif
            call obe_Doppler_av_td_B(ti,tf,n_time_steps,  &
                           rhovec,iunit,popinit)
         else
            call printerror
            print*,'An illegal value of iAorB is specified'
            print*,'in the controlparameters file, or'
            print*,'no value is specified. ',iAorB
         endif
      elseif(iDoppler.eq.0)then
            istore = 0
            call obe_tdint(irate,ti,tf,n_time_steps,imethod,  &
                           rhovec,istore,iunit,popinit)
      else
         call printerror
         print*,'The value of iDoppler specified in'
         print*,'the input file is illegal. ',iDoppler
         stop
      endif
!
      end subroutine sub_calcrho_td_cw

!!!!!!!!!!!!!!!!!!!

      subroutine sub_definetdfields
!
!  Prepare the details of the time-dependent field(s).
!
      implicit none
!
      if(nfields.gt.2)then
         call printerror
         print*,'A calculation for non-cw fields can be done only for'
         print*,'one or two fields.'
         stop
      endif
!
      if(itdfieldsAorB.eq.1)then
!
!  Analytical amplitude profile(s)
!
         if(iforce0(1).eq.0)then
            force0 = .false.
         elseif(iforce0(1).eq.1)then
            force0 = .true.
         else
            call printerror
            print*,'Illegal value of iforce0 for field 1. ',iforce0(1)
            stop
         endif
         if(force0 .and. iinterp.eq.0)then
            call printerror
            print*,'The value given to iinterp in the input file'
            print*,'is in conflict with the value given to iforce0'
            print*,'for field 1.'
            stop
         endif
         call mbe_set_envlp('probe',pulse_type(1),tw(1),t0(1),t1(1),force0)
         if(nfields.eq.2)then
            if(iforce0(2).eq.0)then
               force0 = .false.
            elseif(iforce0(2).eq.1)then
               force0 = .true.
            else
               call printerror
               print*,'Illegal value of iforce0 for field 2. ',iforce0(2)
               stop
            endif
            if(force0 .and. iinterp.eq.0)then
               call printerror
               print*,'The value given to iinterp in the input file'
               print*,'is in conflict with the value given to iforce0'
               print*,'for field 2.'
               stop
            endif
            call mbe_set_envlp('coupl',pulse_type(2),tw(2),t0(2),t1(2),force0)
         endif
!
         if(nfields.eq.1)then 
            call mbe_set_tdfields_A(ti,tf,n_time_steps,nfields,iinterp,    &
                                    cfields(1)%amplitude)
         else
            call mbe_set_tdfields_A(ti,tf,n_time_steps,nfields,iinterp,    &
               cfields(1)%amplitude,cfields(2)%amplitude)
         endif
!
      elseif(itdfieldsAorB.eq.2)then
!
!  The amplitude profile(s) are read from the tdamps_in file 
!
         if(filename_tdamps_in.eq.'undefined')then
            call printerror
            print*,'The name of the file containing the amplitude(s)'
            print*,'of the applied field(s) is not specified in'
            print*,'the input file.'
            stop
         endif
         if(n_time_steps.le.0)then
            call printerror
            print*,'No value or a wrong value n_time_steps is'
            print*,'specified in the controlparams file. ',n_time_steps
            stop
         endif
         allocate(t_mesh(1:n_time_steps+1),prbfield_ampl_arr(1:n_time_steps+1),&
                  stat=istat)
         if(istat.gt.0)then
            call printerror
            print*,'Error at the allocation of t_mesh and the prbfield array.'
            stop
         endif
         if(nfields.eq.2)then
            allocate(cplfield_ampl_arr(1:n_time_steps+1),stat=istat)
            if(istat.gt.0)then
               call printerror
               print*,'Error at the allocation of the cplfield array.'
               stop
            endif
         endif
         inquire(file=filename_tdamps_in,exist=fl_exist)
         if(.not.fl_exist)then
            call printerror
            print*,'The tdamps_in file ',filename_tdamps_in
            print*,'is not found.'
            stop
         endif
         open(1,file=filename_tdamps_in,status='OLD')
         do k=1,n_time_steps+1
            if(nfields.eq.1)then
               if(icmplxfld.eq.0)then
                  read(1,*,end=999,err=999)t_mesh(k),prbre
                  prbfield_ampl_arr(k) = dcmplx(prbre,0.d0)
               else
                  read(1,*,end=999,err=999)t_mesh(k),prbre,prbim
                  prbfield_ampl_arr(k) = dcmplx(prbre,prbim)
               endif
            else
               if(icmplxfld.eq.0)then
                  read(1,*,end=999,err=999)t_mesh(k),prbre,cplre
                  prbfield_ampl_arr(k) = dcmplx(prbre,0.d0)
                  cplfield_ampl_arr(k) = dcmplx(cplre,0.d0)
               else
                  read(1,*,end=999,err=999)t_mesh(k),prbre,prbim,cplre,cplim
                  prbfield_ampl_arr(k) = dcmplx(prbre,prbim)
                  cplfield_ampl_arr(k) = dcmplx(cplre,cplim)
               endif
            endif
         enddo
         close(1)
         print*,'Source file for the amplitude(s) of the field(s): ', &
                 filename_tdamps_in
!
         if(nfields.eq.1)then 
            call mbe_set_tdfields_B(n_time_steps,t_mesh,prbfield_ampl_arr)
         else
            call mbe_set_tdfields_B(n_time_steps,t_mesh,prbfield_ampl_arr,    &
                                    cplfield_ampl_arr)
         endif
!
      else
!
         call printerror
         print*,'itdfieldsAorB must be given a value of either 1'
         print*,'or 2. ',itdfieldsAorB
         stop
!
      endif
!
      return
!
 999  continue
!
      call printerror
      print*,'Error or premature end-of-file when reading the file'
      print*,'containing the time mesh and fields amplitude(s).'
      stop
!
      end subroutine sub_definetdfields

!!!!!!!!!!!!!!!!!!!

      subroutine sub_calcrho_td_notcw
!
!  Integrate the optical Bloch equations for one or several fields with
!  a time-dependent envelope.
!
      implicit none
!
      if(nsubsteps.eq.-1)then
         call printerror
         print*,'This program does not support calculations with'
         print*,'a variable number of intermediate steps between'
         print*,'each grid point (nsubsteps = -1).'
         stop
      endif
      if(istart.eq.-1)then
         call printerror
         print*,'This program does not support the option istart = -1.'
         stop
      endif
!
      if(nfields.eq.1)then 
         call mbe_tdint_1(istart,nsubsteps,imethod,iDoppler,rhovec,iunit, &
                          iunit_ampl,popinit)
      elseif(nfields.eq.2)then 
         call mbe_tdint_2(istart,nsubsteps,imethod,iDoppler,rhovec,iunit, &
                        iunit_ampl,popinit)
      else
         call printerror
         print*,'A calculation for non-cw fields can be done only for'
         print*,'one or two fields.'
         stop
      endif
!
      end subroutine sub_calcrho_td_notcw

!!!!!!!!!!!!!!!!!!!

      subroutine sub_outputrho_td
!
!  Write out the final density matrix at the end of an integration
!  in time, unless iprintrho = 0. The density matrix is written at each
!  time step by the relevant obe or ldbl routine when iunit = -1 or
!  iunit > 0.
!
      implicit none
!
      if(iprintrho.ne.0)then
         print 1000
         do j = nmin,nmax
            do i = nmin,j
               if(i.eq.j)then
                  call obe_pop_index(i,m)
                  print 1001,i,j,rhovec(m),0.d0
               else
                  call obe_coher_index(i,j,mre,mim)
                  print 1001,i,j,rhovec(mre),rhovec(mim)
               endif
            enddo
         enddo
      endif
!
 1000 format(/3x,'i',3x,'j',3x,'Re rho(i,j)',3x,'Im rho(i,j)'/)
 1001 format(1x,i3,1x,i3,2x,2(1pe12.5,2x))
!
      end subroutine sub_outputrho_td

!!!!!!!!!!!!!!!!!!!

      subroutine sub_closeoutputfiles_td
!
!  Close the files used to output specific elements of the
!  density matrix, if any were opened.
!
      implicit none
!
      if(any(filename_rho_out(:,:).ne.'undefined'))then
!
         do j = nmin,nmax
            do i = nmin,nmax
               if(iunits_arr(i,j).gt.0)close(iunits_arr(i,j))
            enddo
         enddo
      endif
!
!  Close the file used to output the whole density matrix
!  if this file exists.
!
      if(iunit.gt.0)close(iunit)
!
!  Close the file used to output the field(s) 
!  if this file exists.
!
      if(iunit_ampl.gt.0)close(iunit_ampl)
!
      end subroutine sub_closeoutputfiles_td

!!!!!!!!!!!!!!!!!!!

      subroutine sub_preparedetunings_steady
!
!  Prepare the list of detunings for which the density matrix needs to
!  be calculated (steady state calculations only).
!
      implicit none
!
      if(ivarydetuning.eq.1 .or. ivarydetuning.eq.2)then
         if(ivarydetuning.eq.2 .and. nfields.ne.1)then
            call printerror
            print*,'The option ivarydetuning = 2 cannot be used in calculations'
            print*,'involving more than one fields.'
            stop
         endif
         if(index_field.lt.1 .or. index_field.gt.nfields)then
            call printerror
            print*,'No value or a nonsensical value of index_field is'
            print*,'specified in the input file although ivarydetuning .ne. 0.'
            print*,index_field
            stop
         endif
         ndmax = (detuning_max-detuning_min)/detuning_step + 1
         if(ndmax.gt.10000)then
            call printerror
            print*,'The calculation is required for a suspiciously large'
            print*,'number of detuning values. Check the values of'
            print*,'detuning_min, detuning_max and detuning_step specified'
            print*,'in the input file.'
            print*,detuning_max,detuning_min,detuning_step
            stop
         endif
         if(ndmax.le.0)then
            call printerror
            print*,'The number of detunings for which the density matrix'
            print*,'must be calculated is not positive. Check the values of'
            print*,'detuning_min, detuning_max and detuning_step specified'
            print*,'in the input file.'
            print*,detuning_max,detuning_min,detuning_step
            stop
         endif
         allocate(detuningsarr(1:ndmax),stat=istat)
         if(istat.gt.0)then
            call printerror
            print*,'Error at the allocation of the detunings array.'
            stop
         endif
         do nd = 1,ndmax
            detuningsarr(nd) = detuning_min + dble(nd-1)*detuning_step
         enddo
         detuningsarr(ndmax) = detuning_max
      elseif(ivarydetuning.eq.0)then
         ndmax = 1
      else
         call printerror
         print*,'Illegal value of ivarydetuning. ',ivarydetuning
         stop
      endif
!
      end subroutine sub_preparedetunings_steady

!!!!!!!!!!!!!!!!!!!

      subroutine sub_preparechiout
!
!  Prepare the output file used to write out chi as necessary. This
!  file replaces any existing one of the same name if a value of 0
!  is specified for the control variable iappend in the input file.
!
      implicit none
!
      if(filename_chi_out.ne.'undefined' .and. isuscept.eq.1)then
!
         fl_no_chi_outfile = .false.
!
         if(iappend.eq.0)then
            inquire(file=filename_chi_out,exist=fl_exist)
            if(fl_exist)then
               open(1,file=filename_chi_out,status='REPLACE')
               close(1)
            endif
         elseif(iappend.ne.1)then
            call printerror
            print*,'The value of iappend specified in the input file'
            print*,'is illegal, or no value is specified, although'
            print*,'one or several output files are specified.'
            print*,iappend
            stop
         endif
!
         print*,'The complex susceptibility will be written to ', &
                filename_chi_out
!
      else
!
         fl_no_chi_outfile = .true.
!
      endif
!
      if(filename_chi_out.eq.'undefined' .and. isuscept.eq.1   &
         .and. ivarydetuning.ge.1)then
         call printerror
         print*,'Error: The options ivarydetuning = 1 or 2 and'
         print*,'isuscept.eq.1 have both been selected but no'
         print*,'output file has been specified for the susceptibility.'
         stop
      endif
!
      end subroutine sub_preparechiout

!!!!!!!!!!!!!!!!!!!

      subroutine sub_calcrho_steady
!
!  Calculate the steady state density matrix. Start by checking
!  there is a possibility for the program to output the results.
!
      implicit none
!
      if(fl_no_rho_outfile .and. fl_no_chi_outfile .and. ivarydetuning.ge.1)then
         call printerror
         print*,'Error: The option ivarydetuning = 1 or 2 has been'
         print*,'selected but no output file has been specified.'
         stop
      endif
!
      if(ivarydetuning.ge.1)then
         cfields(index_field)%detuning = detuningsarr(nd)
         call obe_reset_detuning(index_field,cfields(index_field))
      endif
!
      if(iladder_wkprb.gt.0 .and. iweakprb.ne.1)then
         call printerror
         print*,'Inconsistency: A value > 0 has been specified for'
         print*,'iladder_wkprb but iweakprb has not been set to 1.'
         stop
      endif
!
      if(iladder_wkprb.eq.1)then
!  The final populations are the same as the initial populations...
         call obe_steadystate_ladder(iDoppler,popinit,rhovec,urms)
      elseif(iladder_wkprb.eq.2)then
         call obe_steadystate_onefld_weakprb(iDoppler,popinit,rhovec,urms)
      elseif(iladder_wkprb.eq.3)then
         call obe_steadystate_onefld_powerbr(iDoppler,popinit,rhovec,urms)
      elseif(iladder_wkprb.eq.0)then
         if(iDoppler.eq.1)then
            if(iDoppler_numer_st.eq.1)then
               call obe_Doppler_av_st_numerical(rhovec,ioption,popinit)
            elseif(iDoppler_numer_st.eq.0)then
               if(ivarydetuning.eq.2)then
                  call obe_steadystate_onefld(iprep,1,rhovec,urms)
               else
                  call obe_Doppler_av_st(rhovec,urms)
               endif
            else
               call printerror
               print*,'The value of iDoppler_numer_st specified in'
               print*,'the input file is illegal. ',iDoppler_numer_st
               stop
            endif
         elseif(iDoppler.eq.0)then
            if(ivarydetuning.eq.2)then
               call obe_steadystate_onefld(iprep,0,rhovec)
            else
               call obe_steadystate(rhovec,ioption,popinit)
            endif
         else
            call printerror
            print*,'The value of iDoppler specified in'
            print*,'the input file is illegal. ',iDoppler
            stop
         endif
      else
         call printerror
         print*,'The value of iladder_wkprb specified in'
         print*,'the input file is illegal. ',iladder_wkprb
         stop
      endif
!
      end subroutine sub_calcrho_steady

!!!!!!!!!!!!!!!!!!!

      subroutine sub_outputrho_steady
!
!  Write out the density matrix for a steady state.
!
      implicit none
!
      if(ivarydetuning.eq.0 .and. iprintrho.eq.1)then
         print 1000
         do j = nmin,nmax
            do i = nmin,j
               if(i.eq.j)then
                  call obe_pop_index(i,m)
                  print 1001,i,j,rhovec(m),0.d0
               else
                  call obe_coher_index(i,j,mre,mim)
                  print 1001,i,j,rhovec(mre),rhovec(mim)
               endif
            enddo
         enddo
      endif
      if(any(filename_rho_out(:,:).ne.'undefined'))then
         do j = nmin,nmax
            do i = nmin,nmax
               if(filename_rho_out(i,j).ne.'undefined')then
                  open(1,file=filename_rho_out(i,j),status='UNKNOWN', &
                       position='APPEND')
                  if(i.eq.j)then
                     call obe_pop_index(i,m)
                     if(ivarydetuning.eq.0)then
                        write(1,1500)rhovec(m)
                     else
                        write(1,1501)detuningsarr(nd),rhovec(m)
                     endif
                  elseif(i.lt.j)then
                     call obe_coher_index(i,j,mre,mim)
                     if(ivarydetuning.eq.0)then
                        write(1,1501)rhovec(mre),rhovec(mim)
                     else
                        write(1,1502)detuningsarr(nd),rhovec(mre),rhovec(mim)
                     endif
                  else
                     call obe_coher_index(j,i,mre,mim)
                     if(ivarydetuning.eq.0)then
                        write(1,1501)rhovec(mre),-rhovec(mim)
                     else
                        write(1,1502)detuningsarr(nd),rhovec(mre),-rhovec(mim)
                     endif
                  endif
                  close(1)
               endif
            enddo
         enddo
!
      endif
!
 1000 format(/3x,'i',3x,'j',4x,'Re rho(i,j)',3x,'Im rho(i,j)'/)
 1001 format(1x,i3,1x,i3,2x,2(1pe12.5,2x))
 1500 format(1x,1(1pe12.5,2x))
 1501 format(1x,2(1pe12.5,2x))
 1502 format(1x,3(1pe12.5,2x))
!
      end subroutine sub_outputrho_steady

!!!!!!!!!!!!!!!!!!!

      subroutine sub_suscept
!
!  Calculate the susceptibility, refractive index and absorption
!  coefficient for the probe field (field 1) if required. To avoid
!  the risk of an inconsistency between the data used to generate the
!  density matrix and those required to calculate the susceptibility,
!  the calculation is done only if the density matrix has been
!  calculated with iRabi = 0.
!
      implicit none
!
      if(isuscept.eq.1)then
         if(iRabi.eq.1)then
            call printerror
            print*,'This program allows a calculation of the'
            print*,'susceptibility only for iRabi = 0.'
            stop
         endif
         call obe_susceptibility(cfields(1),rhovec,density,       &
                                            chi,refr_index,alpha)
         if(ivarydetuning.eq.0)then
            print 2000
            print 2001,dreal(chi),dimag(chi),refr_index,alpha
         endif
         if(filename_chi_out.ne.'undefined')then
            open(1,file=filename_chi_out,status='UNKNOWN',        &
                 position='APPEND')
            write(1,2501)detuningsarr(nd),dreal(chi),dimag(chi),  &
                         refr_index,alpha
            close(1)
         endif
      elseif(isuscept.ne.0)then
         call printerror
         print*,'The value of isuscept specified in'
         print*,'the input file is illegal. ',isuscept
         stop
      endif
!
 2000 format(/5x,'Re chi',8x,'Im chi',10x,'n',11x,'alpha'/)
 2001 format(1x,4(1pe12.5,2x))
 2501 format(1x,5(1pe12.5,2x))
!
      end subroutine sub_suscept

!!!!!!!!!!!!!!!!!!!

      subroutine sub_calcrho_mb
!
!  Propagate one or two fields by integrating the Maxwell-Bloch equations.
!
      implicit none
!
      if(iRabi.eq.1)then
         call printerror
         print*,'The field(s) must be defined in terms of field amplitudes'
         print*,'and dipole moments, not Rabi frequencies, for propagation'
         print*,'calculations.'
         stop
      endif
      if(iDoppler.eq.1 .and. iunit.ne.0)then
         call printerror
         print*,'The density matrix can be written out by this program'
         print*,'only for calculations without Doppler averaging.'
         stop
      endif
      if(iunit.eq.0 .and. iunit_ampl.le.0)then
         call printerror
         print*,'No output file has been defined.'
         stop
      endif
!
!  Define the applied field(s):
!
      call sub_definetdfields
!
!  The tolerance parameters used by DOP853 are first passed to obe if
!  this solver is to be used:
      if(imethod.eq.4)call obe_set_tol_dop853(rtol,atol)
!
      if(nfields.eq.1)then
!
         call mbe_propagate_1(istart,popinit,imethod,nsubsteps,izrule,  &
            zmax,n_z_steps,cfields(1),density,iDoppler,                 &
            nz_writeout,nt_writeout,extsub_outpr)
 
      elseif(nfields.eq.2)then
!
         call mbe_propagate_2(istart,popinit,imethod,nsubsteps,izrule,  &
            zmax,n_z_steps,cfields(1),cfields(2),density,iDoppler,  &
            nz_writeout,nt_writeout,extsub_outpr)
!
      else
!
         call printerror
         print*,'A propagation calculation can be done only for'
         print*,'one or two fields.'
         stop
!
      endif
!
      end subroutine sub_calcrho_mb

!!!!!!!!!!!!!!!!!!!

      subroutine printerror
!
      print*,' '
      print*,'Error...'
!
      end subroutine printerror

!!!!!!!!!!!!!!!!!!!

      end program driveall

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine extsub_outpr(i_in,k_in,iv,z,tp,rhovec,F_p_re,F_p_im, &
                              F_c_re,F_c_im)
!
!  To print the output from mbe_propagate_1 and mbe_propagate_2.
!
!  This subroutine does not print out Doppler averaged populations or
!  coherences.
!
!  Warning: The array iunits_arr is obtained and the values of iDoppler,
!  iunformatted, iunit, iunit_ampl and nfields are checked only at the first 
!  call to this subroutine. Subsequent changes to these variables will give
!  unpredictable results.
!
      use general_settings
      use obe
      use ldbl
      implicit none
      
      real(kd), dimension(nst*nst), intent(in) :: rhovec
      real(kd), intent(in) :: tp,z
      real(kd), intent(in) :: F_p_re,F_p_im
      real(kd), intent(in), optional :: F_c_re,F_c_im
      integer, intent(in) :: i_in,k_in,iv
      integer, allocatable, dimension(:,:), save :: iunits_arr
      integer :: i,iost,j,m,mre,mim
      logical, save :: fl_outputfiles
      logical :: fl
!
      character(len=11), save :: file_format
      character(len=24), save :: fmt
      logical, save :: fl_check = .true.
!
      integer :: iDoppler,iunformatted,iunit,iunit_ampl,nfields
      common /common_driveall/ iDoppler,iunformatted,iunit,iunit_ampl,nfields
!
!  Check that the values of iunit and iunit_ampl make sense (for efficiency,
!  this is done only at the first call to this subroutine).
!
      if(fl_check)then
!
         if(nfields.eq.2)then
            if(.not.present(F_c_re) .or. .not.present(F_c_im))then
               print*,''
               print*,'Error: The coupling field is not passed to' 
               print*,'extsub_outpr as expected for a 2-field calculation.'
               stop
            endif
         endif
!
         if(iunit.gt.0)then
!
            if(iDoppler.eq.1)then
               print*,''
               print*,'extsub_outpr: The density matrix written out'
               print*,'by this program is not Doppler averaged.'
               stop
            endif
            if(iunit.eq.5 .or. iunit.eq.6)then
               print*,''
               print*,'extsub_outpr: Wrong unit number. ',iunit
               stop
            endif
            inquire(unit=iunit,exist=fl)
            if(.not.fl)then
               print*,''
               print*,'extsub_outpr: the unit iunit does not exist ',iunit
               stop
            endif
            inquire(unit=iunit,opened=fl)
            if(.not.fl)then
               print*,''
               print*,'extsub_outpr: the unit iunit is not connected',iunit
               stop
            endif
            inquire(unit=iunit,form=file_format)
            if(iunformatted.eq.1)then
               if(file_format.ne.'UNFORMATTED')then
                  print*,''
                  print*,'extsub_outpr: Inconsistent value of iunformatted.'
                  print*,file_format,iunformatted
                  stop
               endif
            elseif(iunformatted.eq.0)then
               if(file_format.eq.'UNFORMATTED')then
                  print*,''
                  print*,'extsub_outpr: Inconsistent value of iunformatted.'
                  print*,file_format,iunformatted
                  stop
               endif
            else
               print*,''
               print*,'extsub_outpr: Illegal value of iunformatted.'
               print*,iunformatted
               stop
            endif
!
            if(iunformatted.eq.0)then
               if(nst*nst+2.le.999999)then
                  write(fmt,'("(1x,",i6,"(1pe18.11,1x))")')nst*nst+2
               else
                  print*,''
                  print*,'The number of states is larger than'
                  print*,'extsub_outpr can cope with in regards to'
                  print*,'preparing the correct format. The code'
                  print*,'should be updated in the unlikely event that'
                  print*,'nst*nst + 2 > 999,999.'
                  stop
               endif
            endif 
!
         elseif(iunit.eq.-1)then
!
            call obe_get_iunits(fl_outputfiles,iunits_arr)
            if(.not.fl_outputfiles)then
               print*,''
               print*,'extsub_outpr: The value of iunit is -1 but the unit'
               print*,'numbers of the output files have not been provided.'
               stop
            endif
            if(iDoppler.eq.1)then
               print*,''
               print*,'extsub_outpr: The density matrix written out'
               print*,'by this program is not Doppler averaged.'
               stop
            endif
!
         elseif(iunit.ne.0)then
!
            print*,''
            print*,'extsub_outpr: illegal value of iunit. ',iunit
            stop
!
         endif
!
         if(iunit_ampl.gt.0)then
            if(iunit_ampl.eq.5 .or. iunit_ampl.eq.6)then
               print*,''
               print*,'Wrong unit number in extsub_outpr for iunit_ampl.'
               print*,iunit_ampl
               stop
            endif
            inquire(unit=iunit_ampl,exist=fl)
            if(.not.fl)then
               print*,''
               print*,'extsub_outpr: the unit iunit_ampl does not exist ', &
                      iunit_ampl
               stop
            endif
            inquire(unit=iunit_ampl,opened=fl)
            if(.not.fl)then
               print*,''
               print*,'extsub_outpr: the unit iunit_ampl is not connected', &
                      iunit_ampl
               stop
            endif
            if(iunit.eq.iunit_ampl .and. iunit.ne.0 .and. iunit.ne.-1)then
               print*,''
               print*,'extsub_outpr: Conflict between iunit and iunit_ampl ', &
                      iunit,iunit_ampl
               stop
            endif
            if(fl_outputfiles)then
               do j = 1,nst
                  do i = 1,nst
                     if(iunits_arr(i,j).eq.iunit_ampl)then
                        print*,''
                        print*,'extsub_outpr: Conflict between iunit_ampl'
                        print*,'and iunits_arr for i,j = ',i,j
                        stop
                     endif
                  enddo
               enddo
            endif
         endif
!
         fl_check = .false.
!
      endif

      if(iunit.gt.0)then
!
         if(file_format.eq.'UNFORMATTED')then
            write(unit=iunit,iostat=iost)tp,z,(rhovec(j),j=1,nst*nst)
         else
            write(iunit,fmt,iostat=iost)tp,z,(rhovec(j),j=1,nst*nst)
         endif
         if(iost.ne.0)then
            print*,''
            print*,'extsub_outpr: write returns with iostat = ',iost
            stop
         endif
!
      elseif(iunit.eq.-1)then
!
         do j = 1,nst
            do i = 1,nst
               if(iunits_arr(i,j).gt.0)then
                  if(i.eq.j)then
                     call ldbl_pop_index(i,m)
                     write(iunits_arr(i,j),1551)tp,z,rhovec(m)
                  elseif(i.lt.j)then
                     call ldbl_coher_index(i,j,mre,mim)
                     write(iunits_arr(i,j),1552)tp,z,rhovec(mre), &
                                                     rhovec(mim)
                  else
                     call ldbl_coher_index(j,i,mre,mim)
                     write(iunits_arr(i,j),1552)tp,z,rhovec(mre), &
                                                    -rhovec(mim)
                  endif
               endif
            enddo
         enddo
      endif
!
      if(iunit_ampl.gt.0)then
!
         if(nfields.eq.1)then
            write(iunit_ampl,1651)tp,z,F_p_re,F_p_im
         elseif(nfields.eq.2)then
            write(iunit_ampl,1652)tp,z,F_p_re,F_p_im,F_c_re,F_c_im
         else
            print*,''
            print*,'Abnormal value of nfields in extsub_outpr.',nfields
            stop
         endif
!
      endif
!
 1551 format(1x,3(1pe12.5,1x))
 1552 format(1x,4(1pe12.5,1x))
 1651 format(1x,4(1pe12.5,1x))
 1652 format(1x,6(1pe12.5,1x))
!
      end subroutine extsub_outpr
