      module obe_constants
!
!  Physical and mathematical constants used by the obe and/or mbe
!  modules and/or associated program units.
!
      implicit none
!
!  The following statement makes the variables defined in obe_constants
!  "visible" to any program using this module.
!
      public   
!
!  CODATA 2014 recommended values:
!
!!!   double precision, parameter :: hbar = 1.054571800d-34
!!!   double precision, parameter :: epsilon0 = 8.854187817d-12
!!!   double precision, parameter :: Boltzk = 1.38064852d-23
!!!   double precision, parameter :: umass = 1.660539040d-27
!
!  CODATA 2018 recommended values:
!
      double precision, parameter :: hbar = 1.054571817d-34
      double precision, parameter :: epsilon0 = 8.8541878128d-12
      double precision, parameter :: Boltzk = 1.380649d-23
      double precision, parameter :: umass = 1.66053906660d-27
!
      end module obe_constants

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      module obe
!
!  The program units contained in this module integrate the optical Bloch
!  equations.
!
!  The parameters defining the system, such as decay rates, detunings, etc,
!  must be passed to obe through calls to the appropriate subroutines.
!
      use general_settings
      use obe_constants
      use ldbl
      implicit none
! 
!  The following statement makes the contents of obe "invisible" to other
!  parts of the programs, with the exceptions of the variables and subprograms
!  declared below as public.
!
      private   
!
!  The index labelling the states varies from nmn to nmx in the arrays
!  passed to obe by the user, while it varies from 1 to nst inside obe.
!  There is no difference if the parameter nmn is 1 as then nmx = nst.
!  The value of the parameter nmn is set in the general_settings module.
      integer, parameter :: nmx = nmn + nst - 1
      integer, save :: noffset = 1 - nmn
!
!------------------------------------------------------------------------------
!
!  The following is public and can be used or accessed from outside the module
!
!  The derived type obecfield, defined below, is used by obe to receive
!  the parameters of the laser fields from external program units. The detunings
!  must be frequencies (not angular frequencies) and must be expressed in MHz.
!  The wavelengths must be expressed in nm. The index idir defines whether
!  the field is considered to propagate towards the positive direction (+1)
!  or the negative direction (-1). This index and the wavelength are used only
!  to calculate Doppler shifts when averaging over a thermal velocity
!  distribution). 
!
!  The component detuning_fact(n) determines how the detuning of the field,
!  Delta, should appear in the (n,n) diagonal element of the Hamiltonian
!  (the energy of state n). In short, it is the factor multiplying Delta
!  (not including the 2 pi factor). Thus if detuning_fact(n) .eq. 0, this
!  diagonal element does not depend on the detuning, whereas if
!  detuning_fact(n) .eq. -1, -Delta is added to this diagonal element.
!
!  The amplitude component is the complex amplitude of the electric field
!  of the wave, expressed in V/m.
!
!  dip_mom: a non-zero value of dip_mom(np,n) signals that states n and np
!  are coupled by the field and state np is higher in energy than state
!  n (as this is how dip_mom needs to be configured). The corresponding
!  transition dipole moment, expressed in units of C m, is taken to be
!  the value of dip_mom.
!
!  Also defined is a derived type obefield, which is essentially the
!  same as obecfield apart that the field amplitude, the dipole moments
!  and the Rabi frequencies are double precision variables rather than
!  double complex variables. The type obefield is not used within obe or
!  mbe. It is defined for use outside these modules, should this be more
!  convenient than using the obecfield type (e.g., because the Rabi
!  frequencies and/or dipole moments and field amplitudes are real, not
!  complex). Input to obe and mbe is through variables of type
!  obecfield, not obefield. Variables can be converted from the obefield
!  type to the obecfield type by using the public subroutine
!  obe_fieldtocfield contained in this module.

      type obefield
         double precision :: detuning
         double precision :: amplitude
         double precision :: wavelength
         double precision, dimension(nmn:nmx) :: detuning_fact
         integer :: idir
         double precision, dimension(nmn:nmx,nmn:nmx) :: dip_mom
         double precision, dimension(nmn:nmx,nmn:nmx) :: Rabif
      end type
      public obefield

      type obecfield
         double precision :: detuning
         double complex :: amplitude
         double precision :: wavelength
         double precision, dimension(nmn:nmx) :: detuning_fact
         integer :: idir
         double complex, dimension(nmn:nmx,nmn:nmx) :: dip_mom
         double complex, dimension(nmn:nmx,nmn:nmx) :: Rabif
      end type
      public obecfield
!
!  The following subprograms contained in obe may be called from outside the
!  module:
!
      public obe_setcsts, obe_setfields, obe_reset_detuning, obe_pop_index,   &
             obe_coher_index, obe_get_campl, obe_reset_campl, obe_setsys,     &
             obe_init_rho, obe_steadystate, obe_tdint, obe_set_ommats,        &
             obe_set_Doppler, obe_Doppler_av_td_A, obe_Doppler_av_td_B,       &
             obe_Doppler_av_st_numerical, obe_susceptibility, obe_weakfield,  &
             obe_2state, obe_weakprb_3stladder, obe_weakprb_4stladder,        &
             obe_get_Dopplerpar, obe_set_tol_dop853, obe_get_tol_dop853,      &
             obe_find_Rabif, obe_find_campl, obe_Doppler_av_st,               &
             obe_steadystate_ladder, obe_fieldtocfield, obe_setoutputfiles,   &
             obe_get_iunits, obe_steadystate_onefld,                          &
             obe_steadystate_onefld_weakprb, obe_steadystate_onefld_powerbr
!
!  All the other subprograms contained in obe can be accessed only from obe.
!
!------------------------------------------------------------------------------
!
!  The following variables will be "known" by all the subprograms
!  contained in this module but are not directly accessible from the outside.
!  To change their values, call the subroutine obe_setcsts, obe_setfields,
!  obe_set_Doppler, obe_reset_detuning or obe_setoutputfiles (with the
!  exceptions of the logical flags, which are set within obe).
!
!  Delta: detunings.(*)
!  campl: complex field amplitudes defining the Rabi frequencies
!  through multiplication by the appropriate dipole moments.
!  wavel_nm: wavelength in nm.
!  idir_st: direction indicator (+1 or -1).
!  (*) These quantities are defined as frequencies (not angular frequencies)
!  and must be expressed in MHz.
      real(kd), dimension(:), allocatable, save :: Delta
      complex(kd), dimension(:), allocatable, save :: campl
      real(kd), dimension(:), allocatable, save :: wavel_nm
      integer, dimension(:), allocatable, save :: idir_st
      integer, dimension(:,:), allocatable, save :: iunits_arr
!
!  The number of fields for the current run (the value of nflds is set
!  by obe_setcsts and cannot be changed once set):
      integer, save :: nflds
!
!  The states which may be initially populated in a calculation within the
!  weak probe approximation (see obe_wkprbapp)
      integer, dimension(nst) :: n_init_wk
!
!  The following quantities are used by obe_setsys to set up the system of
!  equations
      real(kd), dimension (nst), save :: energ_f
      real(kd), dimension (:,:), allocatable, save :: detuning_fact
      complex(kd), dimension (:,:,:), allocatable, save :: dipmom,Rabifreq
      real(kd), dimension (:,:), allocatable, save :: add_dephas_f
      real(kd), dimension (:,:,:), allocatable, save :: Gamma_decay_f
      integer, save :: n_collapse_ops
!
!  Logical variables used to check or control the job flow
      logical, dimension(:), allocatable, save :: fl_fld_not_init
      logical, dimension(:), allocatable, save :: fl_cplg_not_init
      logical, save :: fl_par_not_init = .true.
      logical, save :: fl_initarrays = .true.
      logical, save :: fl_init_ommats = .true.
      logical, save :: fl_mltpls = .false.
      logical, save :: fl_wkprbapp = .false.
      logical, save :: fl_stDopp_init = .true.
      logical, save :: fl_tdDopp_init = .true.
      logical, save :: fl_add_dephas
      logical, save :: fl_Doppler_notinit = .true.
      logical, save :: fl_Doppl_td
      logical, save :: fl_Rabi
      logical, save :: fl_outputfiles=.false.
!
!  The following variables are initialized by obe_Doppler_av_td_A or B, or by
!  obe_Doppler_av_st, or by obe_set_Doppler, and used elsewhere in the module:
      real(kd), dimension(nst*nst), save :: rhovec_init_Dopp
      real(kd), dimension(:,:), allocatable, save :: rhs0v_Dopp, rhs1v_Dopp, &
         rhsv_Dopp 
      real(kd), dimension(:), allocatable, save :: Delta_store_Dopp, vmesh,  &
         fMvweight
      real(kd), save :: t1_Dopp, t2_Dopp, urms_saved  
      integer, save :: imethod_Dopp, ioption_Dopp, irate_Dopp, istore_Dopp,  &
         iunit_Dopp, n_time_steps_Dopp, n_v_values
!
!  Variables used to communicate with the dop853 ODE solver (part of
!  ldbl and mbe)
      double precision, save :: rtol_dop853,atol_dop853
      logical, save :: fl_set_tol_dop853 = .false.
!
!------------------------------------------------------------------------------

      contains
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
      subroutine obe_error
!
!  Control may be passed to this subprogram from elsewhere within obe 
!  following detection of an error.
!
!  At the moment obe_error simply stops the execution without saving files.
!  The user may want to replace this by a more sophisticated error handling
!  procedure.
!
      stop 'Error in the obe module...'
!
      end subroutine obe_error

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine obe_setoutputfiles(iunits_arr_in)
!
!  Receive and store the unit numbers of the files used by the
!  time-dependent routines for outputting the values of selected
!  elements of the density matrix.
!
      implicit none
      integer, dimension(nmn:nmx,nmn:nmx), intent(in) :: iunits_arr_in     
      integer :: i,j,k,istat
!
      if(allocated(iunits_arr))deallocate(iunits_arr)
      allocate(iunits_arr(1:nst,1:nst),stat=istat)
      if(istat.gt.0)then
         print*,'obe_setoutputfiles: Error at the allocation of'
         print*,'iunits_arr.' 
         call obe_error
      endif
!
      iunits_arr = -1
!
      do j = nmn,nmx
         do i = nmn,nmx
            if(iunits_arr_in(i,j).gt.0)then
               call obe_check_iunit(iunits_arr_in(i,j))
               iunits_arr(i+noffset,j+noffset) = iunits_arr_in(i,j)
            endif
         enddo
      enddo
!
      fl_outputfiles = .true.
!
      end subroutine obe_setoutputfiles

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine obe_get_iunits(fl_outputfiles_out,iunits_arr_out)
!
!  Return the units numbers of output files set by obe_setoutputfiles.
!  Warning: the array iunits_arr_out is allocated only if
!  fl_outputfiles_out is .true.
!  
      implicit none
      integer, dimension(:,:), allocatable, intent(out) :: iunits_arr_out
      logical, intent(out) :: fl_outputfiles_out
      integer :: istat
!
      if(fl_outputfiles)then
         if(allocated(iunits_arr_out))deallocate(iunits_arr_out)
         allocate(iunits_arr_out(1:nst,1:nst),stat=istat)
         if(istat.gt.0)then
            print*,'obe_get_iunits: Error at the allocation of iunits_arr_out.'
            call obe_error
         endif
         iunits_arr_out = iunits_arr
         fl_outputfiles_out = .true.
      else
         fl_outputfiles_out = .false.
      endif
!
      end subroutine obe_get_iunits

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine obe_fieldtocfield(field_in,cfield_out)
!
!  Convert the obefield object field_in into the obecfield object
!  cfield_out. The imaginary parts of the double complex components
!  of cfield_out are given a zero value in the conversion.
!
      implicit none
      type(obefield), intent(in) :: field_in
      type(obecfield), intent(out) :: cfield_out
!
      cfield_out%detuning = field_in%detuning
      cfield_out%amplitude = field_in%amplitude
      cfield_out%wavelength = field_in%wavelength
      cfield_out%detuning_fact = field_in%detuning_fact
      cfield_out%idir = field_in%idir
      cfield_out%dip_mom = field_in%dip_mom
      cfield_out%Rabif = field_in%Rabif
!
      end subroutine obe_fieldtocfield

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine obe_pop_index(i_in,m)
!
!  Return the index of the component corresponding to the population rho(i,i)
!  in the 1D array representing the density matrix in vector form.
!
!  obe_pop_index is basically an interface with ldbl_pop_index, but
!  with the difference that the state indexes start at nmn (set in
!  the general_settings module) rather than at 1.
!
      implicit none
      integer, intent(in) :: i_in
      integer, intent(out) :: m
      integer :: i_from1,istart
!
      istart = nmn
!
      i_from1 = i_in + 1 - istart
!
!  Check that the value of i_in makes sense and exits if it doesn't.
      if(i_from1.lt.1 .or. i_from1.gt.nst)then
         print*,'obe_pop_index: i out of bounds. ',i_in
         call obe_error
      endif
!
!  Get the result.
      call ldbl_pop_index(i_from1,m)
!
      end subroutine obe_pop_index

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine obe_coher_index(i_in,j_in,mrl,mim)
!
!  Return the indexes of the components corresponding to the real and the
!  imaginary parts of the coherence rho(i,j) (i < j) in the 1D array
!  representing the density matrix in vector form.
!
!  obe_coher_index is basically an interface with ldbl_coher_index, but
!  with the difference that the state indexes start at nmn (set in
!  the general_settings module) rather than at 1.
!
      implicit none
      integer, intent(in) :: i_in,j_in
      integer, intent(out) :: mrl,mim
      integer :: i_from1,j_from1,istart
!
      istart = nmn
!
      i_from1 = i_in + 1 - istart
      j_from1 = j_in + 1 - istart
!
!  Check that the values of i_in and j_in make sense and exits if they don't.
      if(i_from1.lt.1 .or. i_from1.gt.nst .or. j_from1.lt.1 .or.    &
         j_from1.gt.nst .or. j_from1.le.i_from1)then
         print*,'obe_coher_index was called with an illegal value of'
         print*,'i or of j: ',i_in,j_in
         call obe_error
      endif
!
!  Get the results.
      call ldbl_coher_index(i_from1,j_from1,mrl,mim)
!
      end subroutine obe_coher_index

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
      subroutine obe_setsys(h_mat,decrates_mat,t2dephas_mat,ncoll)
!
!  Define the system, using the parameters previously passed to obe
!  through obe_setcsts and obe_setfields.
!
!  See comments below about the setting of hmatdiag_f and rabi_f.
!
!  obe_setsys and the subroutine obe_getmats it contains also make use of
!  the array Gamma_decay_f, which is expected to have been populated by
!  obe_setcsts. Unless collapse operators are passed to obe_setcsts,
!  Gamma_decay_f(i,j,1) should be set to the rate Gamma for the
!  decay of state j to state i, expressed in MHz (here Gamma is a frequency).
!  The corresponding decay rates expressed as angular frequencies will be set
!  to 2 pi times Gamma_decay(i,j,:) when needed. Otherwise Gamma_decay_f(:,:,n)
!  should contain the n-th collapse operator in matrix form.
!
      implicit none
      complex(kd), dimension (nst,nst), intent(out) :: h_mat
      real(kd), dimension (:,:,:), allocatable, intent(out) :: decrates_mat
      complex(kd), dimension (nst,nst), intent(out) :: t2dephas_mat
      integer, intent(out) :: ncoll
      real(kd), dimension (nst):: hmatdiag_f
      complex(kd), dimension (nst,nst) :: rabi_f
      real(kd) :: twopi,conv
      integer :: i,j,k
!
!  Check that nst is not too small (parts of the program would be incorrect
!  if nst.eq.1)
      if(nst.lt.2)then
         print*,'obe_setsys: nst must be at least 2'
         call obe_error
      endif
!
!  Check that the arrays energ_f and Gamma_decay_f
!  have been initialized
      if(fl_par_not_init)then
         print*,'Attempt at using obe_setsys before a call to obe_setcsts'
         call obe_error
      endif
!
!  Check that Delta, campl, detuning_fact and dipmom (or Rabifreq) have been
!  initialized for each field in the problem.
      do k = 1,nflds
         if(fl_fld_not_init(k))then
            print*,'Attempt at using obe_setsys before field ',k, &
                   'was defined through a call to obe_setfields'
            call obe_error
         endif
         if(fl_cplg_not_init(k))then
            print*,'Attempt at using obe_setsys before detuning_fact'
            print*,'and dipmom or Rabifreq were defined for field ',k
            call obe_error
         endif
      enddo
!
!  hmatdiag_f(j) should be set to the diagonal element of the Hamiltonian
!  matrix for state j, expressed in MHz. The corresponding angular frequency
!  will be set to 2 pi times hmatdiag_f(j) by obe_getmats.
!
      do j=1,nst
!  First the energies:
         hmatdiag_f(j) = energ_f(j)
!  Add the detunings:
         do k=1,nflds 
            hmatdiag_f(j) = hmatdiag_f(j) + detuning_fact(j,k)*Delta(k)
         enddo
      enddo
!
!  rabi_f(i,j) should be set to the complex Rabi frequency Omega for the 
!  coupling of state j to state i expressed in MHz. The corresponding 
!  Rabi angular frequency will be set to 2 pi times rabi_f(i,j), which
!  is meant to correspond to the Omega multiplying the term in |i><j| in
!  the Hamiltonian. The Rabi frequency for the corresponding |j><i| term
!  is the complex conjugate of the Rabi frequency for the |i><j| term.
!  Only the latter needs to be specified here; the former
!  (involving the complex conjugate of campl) will be calculated as required
!  by the subroutine obe_getmats.
!
!  The hbar variable is defined in the module obe_constants.
!  The dipole moments passed to obe are assumed to be expressed in units of
!  C m. The reduced Planck constant hbar is expressed in units of J s.
!  If the Rabi frequencies to be used are those directly set through the
!  Rabif component of the field variables, those must be expressed in MHz.
!
      twopi = 8.0_kd*atan(1.0_kd)
      conv = 1.0e-6_kd/(twopi*hbar) ! The 1e-6 is for the conversion to MHz.
      rabi_f = (0.0_kd,0.0_kd)
      do k = 1,nflds
         do j = 1,nst
            do i = 1,nst
               if(fl_Rabi)then
                  rabi_f(i,j) = rabi_f(i,j) + Rabifreq(i,j,k)
               else
                  rabi_f(i,j) = rabi_f(i,j) + campl(k)*dipmom(i,j,k)*conv
               endif
            enddo
         enddo
      enddo
!
!  We now call obe_getmats, which builds the arrays h_mat, decrates_mat
!  and t2dephas_mat given the above information. The control returns to the
!  calling program once obe_getmats is finished.
      call obe_getmats
 
      contains

!!!!!!!!!!!!!!!!!!

      subroutine obe_getmats
!
!  Construct the matrices returned by obe_setsys (this includes converting
!  from frequencies to angular frequencies)
!
      implicit none
      integer :: i,istat,j
!
!  The twopi variable is defined in obe_setsys.
!
!  Construct the Hamiltonian matrix.
!  First diagonal elements:
      h_mat = (0.0_kd,0.0_kd)
      do i=1,nst
         h_mat(i,i) = cmplx(twopi*hmatdiag_f(i),0.0_kd,kd)
      enddo
!  Then lower and upper triangles (note the minus sign):
      do j=1,nst
         do i=1,nst
            if(i.eq.j)cycle   ! No dipole coupling of a state with itself
            if(fl_Rabi)then
               h_mat(i,j) = -twopi*rabi_f(i,j)/2.0_kd
            else
               if(rabi_f(i,j).ne.(0.0_kd,0.0_kd))then
                  h_mat(i,j) = -twopi*rabi_f(i,j)/2.0_kd
                  h_mat(j,i) = conjg(h_mat(i,j))
               endif
            endif
         enddo
      enddo
!
!  Construct the array containing the decay rates (expressed as angular
!  frequencies), for use by ldbl_set_superopmat.
      if(allocated(decrates_mat))deallocate(decrates_mat)
      if(n_collapse_ops .eq. -1)then
         allocate(decrates_mat(nst,nst,1),stat=istat)
      elseif(n_collapse_ops .ge. 1)then
         allocate(decrates_mat(nst,nst,n_collapse_ops),stat=istat)
      else
         print*,'obe_getmats: The value of n_collapse_ops is invalid.'
         print*,n_collapse_ops
         call obe_error
      endif
      if(istat.gt.0)then
         print*,'obe_getmats: Error at the allocation of decrates_mat'
         call obe_error
      endif
      ncoll = n_collapse_ops
      if(ncoll.eq.-1)then
         decrates_mat = twopi*Gamma_decay_f
      else
         decrates_mat = Gamma_decay_f  ! Note: There is no factor of 2 pi here
                                       ! as the entries of the user-provided
                                       ! collapse operator are meant to be the
                                       ! square roots of angular frequencies.
      endif
!
!  Define the matrix containing the additional dephasing rates expressed
!  as frequencies (these additional dephasing rates are meant to
!  be the rates of decay of the coherences, not including the contribution
!  of the decays of the populations - e.g., the rates of decay arising from
!  T2 processes such as collisional broadening). The decay rate of the (i,j)
!  coherence is taken to be the same as that of the (j,i) coherence. The
!  rate for either of these two coherences must be specified. If a rate
!  is specified in add_dephas for one of these two coherences, then the
!  element of add_dephas corresponding to the other one must be either
!  identical or be zero. The diagonal elements of add_dephas are ignored,
!  even if non-zero.
      t2dephas_mat = (0.0_kd,0.0_kd)
      if(fl_add_dephas)then
         do i=2,nst
            do j=1,i-1
               if(add_dephas_f(i,j).ne.0.0_kd)then
                  if(add_dephas_f(j,i).ne.0.0_kd .and.              &
                     add_dephas_f(j,i).ne.add_dephas_f(i,j) )then
                     print*,'obe_getmats: There seems to be an inconsistency'
                     print*,'in the T2 decays rates specified for the (i,j)'
                     print*,'and (j,i) coherences, with i, j = ',i,j
                     print*,add_dephas_f(j,i),add_dephas_f(i,j)
                     call obe_error
                  endif
                  t2dephas_mat(i,j) = cmplx(twopi*add_dephas_f(i,j),0.0_kd,kd)
                  t2dephas_mat(j,i) = t2dephas_mat(i,j)
               endif
            enddo
         enddo
      endif
!
      end subroutine obe_getmats

!!!!!!!!!!!!!!!!!!

      end subroutine obe_setsys
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine obe_setcsts(en_in_f,decay_in,nflds_in,iwkprbapp,   &
                 iRabi,imltpls,add_dephas,n_coll_in,collapse_in)
!
!  Interface with the user's own program. This subroutine merely copies
!  the content of the variables en_in_f and nflds_in (passed to
!  this subroutine through its arguments) into variables private to the obe
!  module, shifting the indexes as necessary if the starts indexing the 
!  states at a value of nmn not equal to 1.  These private variables
!  are called energ_f and nflds within this module.
!  Likewise, either decay_in or collapse_in is copied into the
!  array Gamma_decay_f, private to obe (see below).
!
!  At the first call, obe_setcsts also allocates a number of arrays whose
!  size depends on the value of nflds.
!
!  Finally, obe_setcsts resets the logical flags fl_par_not_init, fl_initarrays,
!  fl_stDopp_init and and fl_tdDopp_init (these variables are declared at
!  the start of the module).
!
!  iRabi is an optional argument. If present and .eq. 1, then the program
!  does not use the contents of the dip_mom and amplitude components
!  of the fields for calculating the Rabi frequencies, but instead
!  gets the latter directly from their Rabif components. If this argument is
!  present and not equal to 1, then it must be equal to 0 and the Rabif 
!  component is ignored. Setting iRabi = 1 makes it impossible to change
!  or query the complex field amplitudes and to use the subroutines of the mbe
!  module.
!
!  imltpls is an optional argument. If present and .eq. 1, then the program
!  calculates the matrix representing the r.h.s. of the Optical Bloch equations
!  in a way which is more efficient for calculations repeated for multiple
!  values of the complex field amplitudes (but which is slower otherwise).
!  If present and .eq. 0 or not present, this mode of calculation is disabled.
!  iRabi and imltpls cannot be both 1.
!
!  add_dephas is also an optional argument, which can be used to take into
!  account any dephasing not arising from the decay of the populations. This
!  argument is a matrix of dephasing frequencies (not angular frequencies).
!  It is present, the (i,j) element of add_dephas is taken to be the
!  additional decay rate of the (i,j) coherence. The program will multiply
!  this matrix by 2 pi for conversion to angular frequencies.
!
!  The arguments collapse_in and n_coll_in are also optional. If present, then
!  the program assumes that the collapse operators the ldbl module should use
!  are contained in collaps_in in matrix form, and that there are n_coll_in of
!  them (the argument n_coll_in must be provided if collapse_in is). If absent,
!  then the program assumes that the arrays decay_in contains all the necessary
!  information. (The size of decay_in does not matter if the arguments
!  collapse_in and n_coll_in are specified.)
!
      implicit none
      double precision, dimension(nmn:nmx), intent(in) :: en_in_f
      double precision, dimension(nmn:nmx,nmn:nmx), intent(in) :: decay_in
      double precision, dimension(nmn:nmx,nmn:nmx), intent(in), optional ::   &
                                                                  add_dephas
      double precision, dimension(:,:,:), intent(in), optional :: collapse_in
      integer, intent(in) :: iwkprbapp,nflds_in
      integer, intent(in), optional :: imltpls,iRabi,n_coll_in
      integer :: i,j,istat
      logical, save :: fl_first_call = .true.
      logical, save :: fl_Rabi_old = .false.
!
      do i = nmn,nmx
         energ_f(i+noffset) = en_in_f(i)
      enddo
!
!  Check that the array add_dephas, if present, is formatted in a meaningful
!  way, and store the corresponding information if required.
      if(present(add_dephas))then
         if(size(add_dephas,1).eq.nst .and. size(add_dephas,2).eq.nst)then
            if(allocated(add_dephas_f))deallocate(add_dephas_f)
            allocate(add_dephas_f(nst,nst),stat=istat)
            if(istat.gt.0)then
               print*,'obe_setcsts: Error at the allocation of add_dephas_f.'
               call obe_error
            endif
            add_dephas_f = 0.0_kd
            do j = nmn,nmx
               do i = nmn,nmx
                  if(add_dephas(i,j).ne.0.0_kd)then
                     if(add_dephas(j,i).ne.0.0_kd .and.              &
                        add_dephas(j,i).ne.add_dephas(i,j) )then
                        print*,'obe_setcsts: There seems to be an inconsistency'
                        print*,'in the T2 decays rates specified for the (i,j)'
                        print*,'and (j,i) coherences, with i, j = ',i,j
                        print*,add_dephas(j,i),add_dephas(i,j)
                        call obe_error
                     endif
                     add_dephas_f(i+noffset,j+noffset) = add_dephas(i,j)
                     add_dephas_f(j+noffset,i+noffset) = add_dephas(i,j)
                  endif
               enddo
            enddo
            fl_add_dephas = .true.
         else
            print*,'obe_setcsts: The subroutine is called with an add_dephas'
            print*,'argument but this argument is not dimensioned correctly.'
            print*,shape(add_dephas)
            call obe_error
         endif
      else
         fl_add_dephas = .false.
      endif
!
!  Check that the arrays decay_in and (if present) collapse_in are formatted
!  in a meaningful way, and store the corresponding information.
      if(allocated(Gamma_decay_f))deallocate(Gamma_decay_f)
      if(present(n_coll_in))then
         if(n_coll_in.lt.0)then
            print*,'obe_setcsts: Non-sensical value of n_coll_in.'
            print*,n_coll_in
            call obe_error
         elseif(.not.present(collapse_in))then
            print*,'obe_setcsts: An array of collapse operators must be'
            print*,'passed to this subroutine if the n_coll argument is'
            print*,'present.'
            call obe_error 
         endif
      endif
      if(present(n_coll_in))then
         if(n_coll_in.ne.0 .and. nmn.ne.1)then
            print*,'obe_setcsts: The current code can cater for the inclusion'
            print*,'of collapse operators (n_coll_in .ne. 0) only when' 
            print*,'the states are indexed from 1 upwards, i.e., only'
            print*,'for nmn = 1.'
            call obe_error
         endif
      endif
      if(present(collapse_in))then
         if(size(collapse_in,1).eq.nst .and. size(collapse_in,2).eq.nst)then
            if(present(n_coll_in))then
               if(size(collapse_in,3).lt.n_coll_in)then
                  print*,'obe_setcsts: the number of collapse operators'
                  print*,'specified by the input argument n_coll_in is'
                  print*,'inconsistent with the size of the collapse_in array.'
                  print*,size(collapse_in,3),n_coll_in
                  call obe_error
               elseif(n_coll_in.lt.0)then
                  print*,'obe_setcsts: Non-sensical value of n_coll_in.'
                  print*,n_coll_in
                  call obe_error
               else
                  n_collapse_ops = n_coll_in
               endif
            else
               print*,'obe_setcsts: The optional argument n_coll_in must'
               print*,'be specified if collapse_in is.'
               call obe_error
            endif
         else
            print*,'obe_setcsts: Error in the formatting of collapse_in...'
            print*,shape(collapse_in)
            call obe_error
         endif
         allocate(Gamma_decay_f(1:nst,1:nst,n_collapse_ops),stat=istat)
         if(istat.gt.0)then
            print*,'obe_setcsts: Error at the allocation of Gamma_decay_f.'
            call obe_error
         endif
         do j = 1,nst
            do i = 1,nst
               Gamma_decay_f(i,j,1:n_collapse_ops) =                     &
                                       collapse_in(i,j,1:n_collapse_ops)
            enddo
         enddo
      elseif(all(shape(decay_in).eq.(/nst,nst/)))then
         n_collapse_ops = -1
         allocate(Gamma_decay_f(nst,nst,1),stat=istat)
         if(istat.gt.0)then
            print*,'obe_setcsts: Error at the allocation of Gamma_decay_f.'
            call obe_error
         endif
         do j = nmn,nmx
            do i = nmn,nmx
               Gamma_decay_f(i+noffset,j+noffset,1) = decay_in(i,j)
            enddo
         enddo
      else
         print*,'obe_setcsts: No collapse matrices are provided but the'
         print*,'input variable decay_in is not a nst by nst array.'
         print*,shape(decay_in)
         call obe_error
      endif
!
!  Set fl_Rabi, to signal whether or not the Rabi frequencies are to be obtained
!  directly from the Rabif component of the fields. The default is .false.
      fl_Rabi = .false.
      if(present(iRabi))then
         if(iRabi.eq.1)then
            fl_Rabi = .true.
         elseif(iRabi.ne.0)then
            print*,'obs_setcsts: The iRabi argument is present but has'
            print*,'an illegal value. ',iRabi
            call obe_error
         endif
      endif
!
      if(fl_first_call)then
         if(nflds_in.lt.1)then
            print*,'The number of fields specified in a call to obe_setcsts'
            print*,'is wrong. ',nflds_in
            call obe_error
         endif
         nflds = nflds_in
!  Delta, campl, wavel_nm, idir_st, fl_fld_not_init are global
!  variables within obe.
         if(allocated(detuning_fact))deallocate(detuning_fact)
         if(allocated(Delta))deallocate(Delta)
         if(allocated(campl))deallocate(campl)
         if(allocated(wavel_nm))deallocate(wavel_nm)
         if(allocated(idir_st))deallocate(idir_st)
         if(allocated(fl_fld_not_init))deallocate(fl_fld_not_init)
         if(allocated(fl_cplg_not_init))deallocate(fl_cplg_not_init)
         allocate(detuning_fact(nst,nflds),                                &
            Delta(nflds), campl(nflds), wavel_nm(nflds), idir_st(nflds),  &
            fl_fld_not_init(nflds), fl_cplg_not_init(nflds),stat=istat)
         if(istat.gt.0)then
            print*,'obe_setcsts: Error at the allocation of detuning_fact etc.'
            call obe_error
         endif
         if(fl_Rabi)then
            allocate(Rabifreq(nst,nst,nflds), stat=istat)
         else
            allocate(dipmom(nst,nst,nflds), stat=istat)
         endif
         if(istat.gt.0)then
            print*,'obe_setcsts: Error at the allocation of dipmom or Rabifreq.'
            call obe_error
         endif
         fl_Rabi_old = fl_Rabi
         fl_fld_not_init = .true.
         fl_first_call = .false.
      else
         if(nflds_in.ne.nflds)then
            print*,'Error in obe_setcsts: '
            print*,'Attempt to update the number of fields.'
            print*,nflds,nflds_in
            call obe_error
         endif
         if(fl_Rabi.neqv.fl_Rabi_old)then
            print*,'Error in obe_setcsts: '
            print*,'Attempt to update the Rabi frequency option.'
            call obe_error
         endif
      endif
!
      call obe_wkprbapp
!
!  Set the flag fl_par_not_init to .false. to signal that the arrays
!  energ_f and Gamma_decay_f and the variable nflds are no longer
!  un-initialized.
      fl_par_not_init = .false.
!
!  Set the flag fl_mltpls to .true. if the calculation is meant to be
!  repeated for varying values of the field amplitude(s). This flag
!  is used by obe_calc_rhsmat.
      if(present(imltpls))then
         if(imltpls.eq.1)then
            if(fl_Rabi)then
               print*,'obe_setcsts: imltpls.eq.1 is incompatible with'
               print*,'iRabi.eq.1.'
               call obe_error
            endif
            fl_mltpls = .true.
         elseif(imltpls.eq.0)then
            fl_mltpls = .false.
         else
            print*,'obe_setcsts: Illegal value of imltpls. ',imltpls
            call obe_error
         endif
      else
         fl_mltpls = .false.
      endif
!
!  Sets the flags fl_initarrays, fl_init_ommats, fl_stDopp_init and
!  fl_tdDopp_init to .true., to force the re-initialization of the arrays used
!  by obe_calc_rhsmat, obe_Doppler_timedep, obe_steadystate_setarrays or
!  obe_Doppler_steadystate.
      fl_initarrays = .true.
      fl_init_ommats = .true.
      fl_stDopp_init = .true.
      fl_tdDopp_init = .true.
!
      contains
 
!!!!!!!!!!!!!!!!!!!

      subroutine obe_wkprbapp
!
!  To start or stop working within the weak probe approximation.
!
!  obe_wkrpbapp sets the logical variable fl_wkprbapp to .true. 
!  if iwkprbapp.eq.1, or to .false. if iwkprbapp.eq.0
!  (this logical variable is global within the module and is used by
!  obe_set_rhsmat). A message is printed on the standard output when the
!  value of fl_wkprbapp is changed (initially, its value is .false.).
!
      implicit none
!
      if(iwkprbapp.eq.0)then
         if(fl_wkprbapp)then
!!          print*,'The weak probe approximation is no longer assumed.'
            fl_wkprbapp = .false.
         endif
      elseif(iwkprbapp.eq.1)then
!!       print*,'The weak probe approximation is now assumed.'
         fl_wkprbapp = .true.
      else
         print*,'Illegal value of iwkprbapp at a call to obe_wkprbapp', &
                iwkprbapp
         call obe_error
      endif
!
      end subroutine obe_wkprbapp
 
!!!!!!!!!!!!!!!!!!!

      end subroutine obe_setcsts
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine obe_setfields(k,cfield)
!
!  Interface with the user's own program. This subroutine receives the
!  details of field k through the variable called cfield (of type obecfield)
!  and copies them into variables private to obe. The latter, as well as the
!  derived type obecfield, are defined in the declarations stated at the
!  beginning of the module obe.
!
!  In particular, this subroutine constructs the arrays dipmom(:,:,k) [or
!  Rabifreq(:,:,k)] and detuning_fact(:,k) according to the information passed
!  through the dip_mom and detuning_fact components of cfield.
!  dip_mom(np,n) must be non-zero if and only if the states n and np are coupled
!  by field k and the correct complex Rabi frequency is obtained by multiplying
!  cfield%dip_mom(np,n) by the complex field amplitude cfield%amplitude. Thus
!  cfield%dip_mom(n,np) must be zero if cfield%dip_mom(np,n) is non-zero.
!  Constructing dip_mom the other way round may lead to wrong results.
!
!  obe_setfields also resets the logical flags fl_fld_not_init and
!  fl_cplg_not_init for field k and the logical flags fl_initarrays,
!  fl_stDopp_init and fl_tdDopp_init.
!
!  This subroutine should not be called before the first call to the 
!  subroutine obe_setcsts, which sets the number of fields and allocate
!  the relevant arrays!
!
      implicit none
      integer, intent(in) :: k
      type(obecfield), intent(in) :: cfield
      integer :: n,np
!
!  Check that obe_setcsts has already been called; exit if hasn't.
      if(fl_par_not_init)then
         print*,'obe_setfields was called before the first call to obe_setcsts'
         call obe_error
      endif
!
!  Check that the index of the cfield is legal
      if(k.lt.1 .or. k.gt.nflds)then
         print*,'obe_setfields was called with an illegal value of k', k
         call obe_error
      endif
!
      Delta(k) = cfield%detuning
      campl(k) = cfield%amplitude
      wavel_nm(k) = cfield%wavelength
      idir_st(k) = cfield%idir
!
      do n = nmn,nmx
         detuning_fact(n+noffset,k) = real(cfield%detuning_fact(n),kd)
      enddo
!
      if(fl_Rabi)then
         Rabifreq(:,:,k) = (0.0_kd,0.0_kd)
      else
         dipmom(:,:,k) = (0.0_kd,0.0_kd)
      endif
      do n = nmn,nmx
         do np = nmn,nmx
            if(fl_Rabi)then
               if(np.eq.n .and. cfield%Rabif(np,n).ne.(0.d0,0.d0))then
                  print*,'obe_setfields has detected an error for field ',k
                  print*,'A non-zero Rabi frequency is specified for'
                  print*,'the coupling of state ',n,' with itself.'
                  call obe_error
               endif
               if(cfield%Rabif(np,n).ne.(0.d0,0.d0) .and.             &
                  cfield%Rabif(n,np).ne.(0.d0,0.d0) .and.             &
                  cfield%Rabif(np,n).ne.dconjg(cfield%Rabif(n,np)))then
                  print*,'obe_setfields has detected an error for field ',k
                  print*,'Inconsistent Rabi frequencies. ',n,np
                  call obe_error
               endif
               if(cfield%Rabif(np,n).ne.(0.d0,0.d0))then
                  Rabifreq(np+noffset,n+noffset,k) = cfield%Rabif(np,n)
                  Rabifreq(n+noffset,np+noffset,k) = dconjg(cfield%Rabif(np,n))
               endif
            else
               if(cfield%dip_mom(np,n).ne.(0.d0,0.d0))then
                  if(cfield%dip_mom(n,np).ne.(0.d0,0.d0))then
                     print*,'obe_setfields has detected an error for field ',k
                     print*,'the array dip_mom is incorrectly formatted.'
                     print*,n,np,cfield%dip_mom(n,np),cfield%dip_mom(np,n)
                     call obe_error
                  endif
                  dipmom(np+noffset,n+noffset,k) = cfield%dip_mom(np,n)
               endif
            endif
         enddo
      enddo
!
!  Sets the flag fl_fld_not_init(k) to .false. to signal that the k-components
!  of the arrays Delta and campl are no longer un-initialized.
      fl_fld_not_init(k) = .false.
!
!  Sets the flag fl_cplgs_not_init(k) to .false. to signal that the
!  k-blocks of the arrays delta_fact and dipmom
!  are no longer un-initialized.
      fl_cplg_not_init(k) = .false.
!
!  Sets the flags fl_initarrays, fl_init_ommats, fl_stDopp_init and
!  fl_tdDopp_init to .true., to force the re-initialization of the arrays used
!  by obe_calc_rhsmat, obe_Doppler_timedep, obe_steadystate_setarrays or
!  obe_Doppler_steadystate.
      fl_initarrays = .true.
      fl_init_ommats = .true.
      fl_stDopp_init = .true.
      fl_tdDopp_init = .true.
!
      end subroutine obe_setfields

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine obe_reset_detuning(k,cfield)
!
!  Interface with the user's own program. This is an abridged version of
!  obe_setfields, which only updates the detuning for field k. It does not
!  reset fl_fld_not_init, since only part of the details of the field k
!  are updated by this routine. It does not reset fl_initarrays either,
!  avoiding thereby that obe_steadystate_setarrays re-initializes the arrays.
!  (It resets fl_stDopp_init and fl_tdDopp_init, though, as the arrays
!  initialized by obe_Doppler_steadystate or obe_steadystate_timedep need
!  to be updated if the detunings change.)
!
      implicit none
      integer, intent(in) :: k
      type(obecfield), intent(in) :: cfield
!
!  Check that obe_setcsts has already been called; exit if hasn't.
      if(fl_par_not_init)then
         print*,'obe_reset_detuning was called before the first call to'
         print*,'obe_setcsts'
         call obe_error
      endif
!
!  Check that the index of the field is legal
      if(k.lt.1 .or. k.gt.nflds)then
         print*,'obe_reset_detuning was called with an illegal value of k', k
         call obe_error
      endif
!
!  Check that obe_setfields had already been called for this field (the
!  subroutine obe_reset_detuning must only by used to reset detunings,
!  not to set them for the first time).
      if(fl_fld_not_init(k))then
         print*,'obe_reset_detuning was called before the field was set ', k
         call obe_error
      endif
!
      Delta(k) = cfield%detuning
!
      fl_stDopp_init = .true.
      fl_tdDopp_init = .true.
!
      end subroutine obe_reset_detuning

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine obe_get_campl(nflds_out,campl_out)
!
!  Return the number of fields and their complex amplitudes.
!  
      implicit none
      integer, intent(out) :: nflds_out
      double complex, dimension(:), allocatable, intent(out) :: campl_out
      integer :: istat,n
!
      if(nflds.lt.1)then
         print*,'obe_get_campl: It appears that the variable nflds'
         print*,'has not been initialised correctly. ',nflds
         call obe_error
      endif
      if(fl_Rabi)then
         print*,'obe_get_campl cannot be called if the couplings are'
         print*,'defined by Rabi frequencies (iRabi.eq.1).'
         call obe_error
      endif
      if(allocated(campl_out))deallocate(campl_out)
      allocate(campl_out(nflds),stat=istat)
      if(istat.gt.0)then
         print*,'obe_get_campl: Error at the allocation of campl_out.'
         call obe_error
      endif
!
      nflds_out = nflds
      do n = 1,nflds
         if(fl_fld_not_init(n))then
            print*,'obe_get_campl: field ',n,' is not defined.'
            call obe_error
         else
            campl_out(n) = campl(n)
         endif
      enddo
!
      end subroutine obe_get_campl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine obe_reset_campl(k,campl_in)
!
!  Interface with the user's own program. This is an abridged version of
!  obe_setfields, which only updates the complex field amplitude for field k.
!  It does not reset fl_fld_not_init, since only part of the details of the
!  field k are updated by this routine. Contrary to obe_reset_detuning,
!  it resets fl_initarrays. Note that the flag fl_init_ommats is not changed.
!
!  Using this subroutine is incompatible with using the iRabi.eq.1 option
!  of obe_setcsts.
!
      implicit none
      integer, intent(in) :: k
      double complex, intent(in) :: campl_in
!
      if(fl_Rabi)then
         print*,'obe_reset_campl: using this subroutine is incompatible'
         print*,'with using the iRabi.eq.1 option of obe_setcsts.'
         call obe_error
      endif
!
!  Check that obe_setcsts has already been called; exit if hasn't.
      if(fl_par_not_init)then
         print*,'obe_reset_campl was called before the first call to'
         print*,'obe_setcsts'
         call obe_error
      endif
!
!  Check that the index of the field is legal
      if(k.lt.1 .or. k.gt.nflds)then
         print*,'obe_reset_campl was called with an illegal value of k', k
         call obe_error
      endif
!
!  Check that obe_setfields had already been called for this field (the
!  subroutine obe_reset_detuning must only by used to reset detunings,
!  not to set them for the first time).
      if(fl_fld_not_init(k))then
         print*,'obe_reset_campl was called before the field was set ', k
         call obe_error
      endif
!
      campl(k) = campl_in
!
      fl_initarrays = .true.
      fl_stDopp_init = .true.
      fl_tdDopp_init = .true.
!
      end subroutine obe_reset_campl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
      subroutine obe_set_rhsmat(rhs_mat,v)
!
!  Returns the matrix representing the right-hand sides of the optical Bloch
!  equations through the arrays rhs_mat.
!
!  obe_set_rhsmat gets the right-hand sides for the full system from
!  ldbl_set_rhsmat. In case of a calculation assuming the weak probe
!  approximation, these right-hand sides are changed as required. If the weak
!  probe approximation is not assumed, obe_set_rhsmat simply passes on the
!  rhs_mat array received from ldbl_set_rhsmat, without any change.
!
!  In this implementation, the probe field is assumed to be field 1 and
!  the initial state is meant to be state 1...
!
!  v is an optional argument. If present, the detunings are Doppler shifted 
!  for a velocity v before rhs_mat is calculated (the detunings are then reset
!  to their original values). If absent, rhs_mat is calculated without Doppler
!  shift.
!
      implicit none
      integer, parameter :: neqs = nst*nst
      real(kd), intent(in), optional :: v
      real(kd), dimension (neqs,neqs), intent(out) :: rhs_mat
      real(kd), dimension(:), allocatable :: Delta_store
      logical, save :: fl_first_call = .true.
!
      if(present(v) .and. fl_first_call)then
!  Check that obe_setcsts has already been called; exit if hasn't.
         if(fl_par_not_init)then
            print*,'obe_set_rhsmat was called before the first'
            print*,'call to obe_setcsts, for a calculation with Doppler shift'
            call obe_error
         endif
         allocate(Delta_store(nflds))
         fl_first_call = .false.
      endif
!
      if(present(v))then
         Delta_store = Delta
         call obe_Doppler_shift(v)
      endif
!
      call ldbl_set_rhsmat(rhs_mat)
!
!  Transform rhs_mat if the weak probe approximation is assumed:
      if(fl_wkprbapp)call obe_set_wkprb_rhsmat
!
      if(present(v))then
         Delta = Delta_store
      endif
!
      contains

!!!!!!!!!!!!!!!!!!

      subroutine obe_set_wkprb_rhsmat
!
      real(kd), dimension (:,:), allocatable :: rhs_mat0
      complex(kd), dimension (:,:), allocatable :: Rabifreq_store
      integer, dimension (:), allocatable :: k
      complex(kd) :: campl1_store
      integer :: i,istat,j,m
!
!  Get the matrix representing the right-hand sides for a zero probe field.
      allocate(rhs_mat0(neqs,neqs),stat=istat)
      if(istat.gt.0)then
         print*,'obe_set_wkprb_rhsmat: Error at the allocation of rhs_mat0.'
         call obe_error
      endif
      if(fl_Rabi)then
         allocate(Rabifreq_store(nst,nst),stat=istat)
         if(istat.gt.0)then
            print*,'obe_set_wkprb_rhsmat: Error at the allocation'
            print*,'of the array Rabifreq_store.'
            call obe_error
         endif
         Rabifreq_store = Rabifreq(:,:,1)
         Rabifreq(:,:,1) = (0.0_kd,0.0_kd)
         call ldbl_set_rhsmat(rhs_mat0)
         Rabifreq(:,:,1) = Rabifreq_store
         deallocate(Rabifreq_store)
      else
         campl1_store = campl(1)
         campl(1) = (0.0_kd,0.0_kd)
         call ldbl_set_rhsmat(rhs_mat0)
         campl(1) = campl1_store
      endif
!
!  Find out which of the components of the density matrix form group A.
      allocate(k(neqs),stat=istat)
      if(istat.gt.0)then
         print*,'obe_set_wkprb_rhsmat: Error at the allocation of k.'
         call obe_error
      endif
      k = 0
!  First, include the populations in group A:
      do i = 1,nst
         call ldbl_pop_index(i,m)
         k(m) = 1
      enddo
!  Then scan the other states repeatedly, until the entire group A is found
!  (k(i).eq.1 or 2 means that i is in group A):
      i = 1
      do
         if(k(i).eq.1)then
            do j = 1,neqs
               if(k(j).eq.0)then
                  if(rhs_mat0(i,j).ne.0 .or. rhs_mat0(j,i).ne.0)k(j) = 1
               endif
            enddo
            k(i) = 2
            i = 1
         else
            if(i.ge.neqs)exit
            i = i + 1
         endif
      enddo
!  Make sure that the elements of rhs_mat0 forming the AA block do
!  not have a dependence in the probe field:
      do i=1,neqs
         if(k(i).eq.2)then
            do j=1,neqs
               if(k(j).eq.2)then
                  if(rhs_mat0(i,j).ne.rhs_mat(i,j))then
                     print*,'obe_set_wkprb_rhsmat: The A block has'
                     print*,'a dependence in the probe field... ',i,j
                     call obe_error
                  endif
               endif
            enddo
         endif
      enddo
!
!  And, in rhs_mat, replace the columns corresponding to components in the 
!  B group by their zero-field_probe counterpart.
     do i =1,neqs
        if(k(i).ne.2)rhs_mat(:,i) = rhs_mat0(:,i)
     enddo
!
!  That's all... Deallocate the scratch space and return.
      deallocate(k,rhs_mat0)
!
      end subroutine obe_set_wkprb_rhsmat

!!!!!!!!!!!!!!!!!!
 
      end subroutine obe_set_rhsmat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
      subroutine obe_set_ommats(omeg0_mat,omegrl_mat,omegim_mat,omg0_opt,optch)
!
!  Calculate the arrays omeg0_mat, omegrl_mat and omegim_mat such that
!  the matrix representing the right-hand sides, rhs_mat, is given by
!    rhs_mat(:,:) = omeg0_mat(:,:) +
!      sum_j [ Re(omega_j) omegrl_mat(:,:,j) + Im(omega_j) omegim_mat(:,:,j)].
!
!  Also calculate the array omg0_opt if the optional arguments omg0_opt and
!  optch are present. If optch.eq.'v', then rhs_mat(:,:) at a velocity v
!  (for a calculation with Doppler shift) would be
!    rhs_mat(:,:) = omeg0_mat(:,:) + v omg0_opt(:,:,1) +
!      sum_j [ Re(omega_j) omegrl_mat(:,:,j) + Im(omega_j) omegim_mat(:,:,j)],
!  where v is in m/s. If optch.eq.'d', then the detunings are are varied:
!    rhs_mat(:,:) = omeg0_mat(:,:) + sum_j [Delta_j omg0_opt(:,:,j)] +
!      sum_j [ Re(omega_j) omegrl_mat(:,:,j) + Im(omega_j) omegim_mat(:,:,j)],
!
!  These arrays depend on the detunings, on the dephasing rates and possibly
!  also on other parameters of the field (apart from the complex field
!  amplitudes) and must be re-calculated if the value of any of these is
!  changed.
!
!  This subroutine is not used by the other program units contained in obe.
!  It can be called from outside obe, e.g., in the context of a propagation
!  calculation.
!
!  The arrays omegrl_mat and omegim_mat must given the dimension
!  (neqs,neqs,nflds) in the calling program.
!
!  Using this subroutine is incompatible with using the iRabi.eq.1 option
!  of obe_setcsts.
!
!  Programming note: As the dummy arguments of this subprogram include
!  assumed-shape arrays and an optional argument, it will require an interface
!  block in each calling program if taken out of the module.
!
      implicit none
      integer, parameter :: neqs = nst*nst
      real(kd), dimension (neqs,neqs), intent(out) :: omeg0_mat
      real(kd), dimension (:,:,:), intent(out), optional :: omg0_opt
      real(kd), dimension (:,:,:), intent(out) :: omegrl_mat
      real(kd), dimension (:,:,:), intent(out) :: omegim_mat
      character*1, intent(in), optional :: optch
      real(kd), dimension(:), allocatable :: Delta_store
      complex(kd), dimension(:), allocatable :: campl_store
      integer :: istat,k
      logical :: fl_Deltas
!
      if(fl_Rabi)then
         print*,'obe_set_ommats: using this subroutine is incompatible'
         print*,'with using the iRabi.eq.1 option of obe_setcsts.'
         call obe_error
      endif
!
!  Check that obe_setcsts has already been called; exit if hasn't.
      if(fl_par_not_init)then
         print*,'obe_set_ommats was called before the first call to'
         print*,'obe_setcsts'
         call obe_error
      endif
!
!  Check that the arrays omegrl_mat and omegim_mat have the right size
!  in all three dimensions; exit if they haven't.
      if(size(omegrl_mat,1).ne.neqs .or. &
         size(omegrl_mat,2).ne.neqs .or. &
         size(omegrl_mat,3).ne.nflds)then
         print*,'obe_set_ommats: Wrong size for omegrl_mat'
         print*,size(omegrl_mat,1),neqs
         print*,size(omegrl_mat,2),neqs
         print*,size(omegrl_mat,3),nflds
         call obe_error
      elseif(size(omegim_mat,1).ne.neqs .or. &
             size(omegim_mat,2).ne.neqs .or. &
             size(omegim_mat,3).ne.nflds)then
         print*,'obe_set_ommats: Wrong size for omegim_mat'
         print*,size(omegim_mat,1),neqs
         print*,size(omegim_mat,2),neqs
         print*,size(omegim_mat,3),nflds
         call obe_error
      endif
!
!  Check that the array omg0_opt is suitable if present
      if(present(omg0_opt))then
         if(present(optch))then
            if(size(omg0_opt,1).ne.neqs .or. &
               size(omg0_opt,2).ne.neqs)then
               print*,'obe_set_ommats: Wrong size for omg0_opt'
               print*,size(omg0_opt,1),neqs
               print*,size(omg0_opt,2),neqs
               call obe_error
            endif
            if(optch.eq.'v')then
               if(size(omg0_opt,3).ne.1)then
                  print*,'obe_set_ommats: Wrong third size for omg0_opt'
                  print*,size(omg0_opt,3)
                  call obe_error
               endif
            elseif(optch.eq.'d')then
               if(size(omg0_opt,3).ne.nflds)then
                  print*,'obe_set_ommats: Wrong third size for omg0_opt'
                  print*,size(omg0_opt,3)
                  call obe_error
               endif
            else
               print*,'obe_set_ommats: Illegal calculation option. ',optch
               call obe_error
            endif
         else
            print*,'obe_set_ommats: The calculation option is not specified.'
            call obe_error
         endif
      endif
!
!  Store the current values of the complex field amplitudes
!  (recall that campl is a global variable):
      if(.not.allocated(campl_store))allocate(campl_store(nflds),stat=istat)
      if(istat.gt.0)then
         print*,'obe_set_ommats: Error at the allocation of campl_store'
         call obe_error
      endif
      campl_store = campl
!
!  Calculate arrays with zero field amplitudes
      campl = (0.0_kd,0.0_kd)
      if(present(omg0_opt))then
         if(optch.eq.'v')then
            call obe_set_rhsmat(omeg0_mat)
            call obe_set_rhsmat(omg0_opt(:,:,1),1.0_kd)
            omg0_opt(:,:,1) = omg0_opt(:,:,1) - omeg0_mat
            fl_Deltas = .false.
         elseif(optch.eq.'d')then
            if(.not.allocated(Delta_store)) &
               allocate(Delta_store(nflds),stat=istat)
            if(istat.gt.0)then
               print*,'obe_set_ommats: Error at the allocation of Delta_store'
               call obe_error
            endif
            Delta_store = Delta
            Delta = 0.0_kd
            call obe_set_rhsmat(omeg0_mat)
            do k=1,nflds
               Delta = 0.0_kd
               Delta(k) = 100.0_kd
               call obe_set_rhsmat(omg0_opt(:,:,k))
               omg0_opt(:,:,k) = (omg0_opt(:,:,k) - omeg0_mat)/100.0_kd
            enddo
            Delta = 0.d0
            fl_Deltas = .true.
         endif               
      else
         call obe_set_rhsmat(omeg0_mat)
         fl_Deltas = .false.
      endif

!  Calculate the first order matrices, changing the fields one by one
      do k=1,nflds
         campl = (0.0_kd,0.0_kd)
         campl(k) = (100.0_kd,0.0_kd)
         call obe_set_rhsmat(omegrl_mat(:,:,k))
         omegrl_mat(:,:,k) = (omegrl_mat(:,:,k) - omeg0_mat)/100.0_kd
         campl(k) = (0.0_kd,100.0_kd)
         call obe_set_rhsmat(omegim_mat(:,:,k))
         omegim_mat(:,:,k) = (omegim_mat(:,:,k) - omeg0_mat)/100.0_kd
      enddo
!
!  When done, reset the field amplitudes (and, if necessary, the
!  detunings), to their original values
      campl = campl_store
      if(fl_Deltas)Delta = Delta_store
!
      deallocate(campl_store)
      if(fl_Deltas)deallocate(Delta_store)
!
      end subroutine obe_set_ommats

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine obe_calc_rhsmat(rhs_mat)
!
!  Calculate the matrix representing the right-hand sides of the optical
!  Bloch equations, rhs_mat, using
!    rhs_mat(:,:) = omeg0_mat(:,:) +
!      sum_j [ Re(omega_j) omegrl_mat(:,:,j) + Im(omega_j) omegim_mat(:,:,j)] +
!      sum_j [ Delta_j omeg0_d_mat(:,:,j)].
!
!  The arrays appearing in this equation are calculated by obe_set_ommats.
!  A value of .true. for the flag fl_init_ommats is the signal to invoke
!  obe_set_ommats to recalculate these arrays. (This flag is reset to .false.
!  by this subroutine.)
!
!  Using this subroutine is incompatible with using the iRabi.eq.1 option
!  of obe_setcsts.
!
      implicit none
      integer, parameter :: neqs = nst*nst
      real(kd), dimension (:,:), allocatable, save  :: omeg0_mat
      real(kd), dimension (:,:,:), allocatable, save  :: omeg0_d_mat
      real(kd), dimension (:,:,:), allocatable, save  :: omegrl_mat
      real(kd), dimension (:,:,:), allocatable, save  :: omegim_mat
      real(kd), dimension (:,:), allocatable, intent(inout) :: rhs_mat
      integer :: istat,k
!
      if(fl_Rabi)then
         print*,'obe_calc_rhsmat: using this subroutine is incompatible'
         print*,'with using the iRabi.eq.1 option of obe_setcsts.'
         call obe_error
      endif
!
      if(fl_init_ommats)then
!  Get the necessary arrays and reset fl_init_ommats.
         if(.not.allocated(omeg0_mat))    &
            allocate(omeg0_mat(neqs,neqs),stat=istat)
         if(.not.allocated(omeg0_d_mat))  &
            allocate(omeg0_d_mat(neqs,neqs,nflds),stat=istat)
         if(.not.allocated(omegrl_mat))   &
            allocate(omegrl_mat(neqs,neqs,nflds),stat=istat)
         if(.not.allocated(omegim_mat))   &
            allocate(omegim_mat(neqs,neqs,nflds),stat=istat)
         if(istat.gt.0)then
            print*,'obe_calc_rhsmat: Error at the allocation of the arrays.'
            call obe_error
         endif
         call obe_set_ommats(omeg0_mat,omegrl_mat,omegim_mat,omeg0_d_mat,'d')
         fl_init_ommats = .false.
      endif
!
!  Check that rhs_mat has been allocated prior to entry and has the
!  correct size.
      if(allocated(rhs_mat))then
         if(size(rhs_mat,1).ne.neqs .or. size(rhs_mat,2).ne.neqs)then
            print*,'obe_calc_rhsmat was called with a wrongly sized argument'
            print*,size(rhs_mat,1),size(rhs_mat,2)
            call obe_error
         endif
      else
         print*,'obe_calc_rhsmat: rhs_mat is not allocated.'
         call obe_error
      endif
!
!  Calculate the first order matrices, changing the fields one by one
      rhs_mat = omeg0_mat
      do k=1,nflds
         rhs_mat = rhs_mat + Delta(k)*omeg0_d_mat(:,:,k) &
                           + real(campl(k),kd)*omegrl_mat(:,:,k) &
                           + dimag(1.d0*campl(k))*omegim_mat(:,:,k) 
      enddo
!
      end subroutine obe_calc_rhsmat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine obe_init_rho(popinit,rho_vec)
!
!  This subroutine can be used to fill the real, 1D array rho_vec with the
!  density matrix defined by the elements of popinit, assuming all the
!  coherences are zero.
!
!  The initial populations must be passed to obe_init_rho through the array
!  popinit. At exit, the density matrix stored in rho_vec will have these
!  populations as its diagonal elements [i.e., rho_{jj} will be popinit(j)] and
!  all its off-diagonal elements will be zero [rho(i,j) = 0 for i .ne. j].
!
      implicit none
      integer, parameter :: neqs = nst*nst
      double precision, dimension(nmn:nmx), intent(in) :: popinit
      real(kd), dimension(neqs), intent(out) :: rho_vec
      integer :: i,m
!
      rho_vec = 0.0_kd
      do i=nmn,nmx
         call obe_pop_index(i,m)
         rho_vec(m) = real(popinit(i),kd)
      enddo
!
      end subroutine obe_init_rho

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine obe_set_tol_dop853(rtol_in,atol_in)
!
!  This subroutine is used to communicate the values of the tolerance
!  parameters rtol and atol to ldbl and to mbe, for use by the
!  dop853 ode solver. It is a simple mailbox.
!
      implicit none
      double precision, intent(in) :: rtol_in,atol_in
!
!  rtol_dop853, atol_dop853 and fl_set_tol_dop853 are global variables.
      rtol_dop853 = rtol_in
      atol_dop853 = atol_in
      fl_set_tol_dop853 = .true.
!
!  Communicate these values to ldbl, in case of need. mbe will get them
!  from obe, if required, through a call to obe_get_tol_dop853.
      call ldbl_set_tol_dop853(rtol_in,atol_in,fl_set_tol_dop853)
!
      end subroutine obe_set_tol_dop853

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine obe_get_tol_dop853(rtol_out,atol_out,flag_out)
!
!  This subroutine is used by mbe to obtain the values of the
!  tolerance parameters rtol_dop853 and atol_dop853 as well as the
!  logical variable fl_set_tol_dop853. It is a simple mailbox.
!
      implicit none
      double precision, intent(out) :: rtol_out,atol_out
      logical, intent(out) :: flag_out
!
!  rtol_dop853, atol_dop853 and fl_set_tol_dop853 are global variables.
      rtol_out = rtol_dop853
      atol_out = atol_dop853
      flag_out = fl_set_tol_dop853
!
      end subroutine obe_get_tol_dop853

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      subroutine obe_tdint(irate,t1_in,t2_in,n_time_steps,imethod,rho_vec, &
         istore,iunit,popinit)
!
!  A simple interface with the subroutine ldbl_tdint of the ldbl module,
!  for integration of the time-dependent system.
!
!  obe_tdint first obtains the matrix representing the right-hand sides from
!  obe_set_rhsmat and, optionally, (re)computes the initial density matrix.
!  obe_tdint then passes the control to ldbl_tdint for the
!  time-dependent calculation, before returning to the calling program
!  with the density matrix at the final time. Passing through obe_tdint
!  rather than calling ldbl_tdint directly makes it possible to solve the
!  system within the weak probe limit (as set up by obe_set_rhsmat) and
!  to easily define the density matrix at t = t_1.
!
!  The values of istore are restricted to 0 or 1 here. ldbl_tdint also admits
!  istore = -1, for special operations such as the Doppler averaging as done
!  by obe_av_td_A. This option is disabled here, as it does not seem relevant
!  for the intended use of obe_tdint.
!
!  popinit is an optional argument. If absent, the density matrix at t = t_1
!  is taken to be that defined by rho_vec at entry. If present, the initial
!  coherences are set to zero and the initial populations to the values
!  prescribed by popinit.
! 
      implicit none
      integer, parameter :: neqs = nst*nst
      real(kd), dimension (neqs), intent(inout) :: rho_vec
      real(kd) :: t1,t2
      double precision, dimension (nmn:nmx), optional, intent(in) :: popinit
      double precision, intent(in) :: t1_in,t2_in
      integer, intent(in) :: imethod,irate,istore,iunit,n_time_steps
      real(kd), dimension (neqs,neqs) :: rhs_mat
!
      if(istore.lt.0 .or. istore.gt.1)then
         print*,'obe_tdint: The values of istore are limited to 0 and 1 here.'
         print*,'See comments in the Fortran code about istore.eq.-1.'
         call obe_error
      endif
      if(n_time_steps.lt.1)then
         print*,'obe_tdint has been called with a non-sensical value'
         print*,'of n_time_steps. ',n_time_steps
         call obe_error
      endif
!
      call obe_set_rhsmat(rhs_mat)
!
!  Redefines the initial density matrix if popinit is present:
      if(present(popinit))call obe_init_rho(popinit,rho_vec)
!
      t1 = real(t1_in,kd)
      t2 = real(t2_in,kd)
      if(iunit.eq.-1)then
         if(fl_outputfiles)then
            call ldbl_tdint(irate,t1,t2,n_time_steps,imethod,rho_vec,istore, &
               iunit,rhs_mat,iunits_arr,isteady=0)
         else
            print*,'obe_tdint: The value of iunit is -1 but the unit numbers'
            print*,'of the output files have not been provided.'
            call obe_error
         endif
      else
         call ldbl_tdint(irate,t1,t2,n_time_steps,imethod,rho_vec,istore, &
            iunit,rhs_mat,isteady=0)
      endif
!
      end subroutine obe_tdint

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine obe_steadystate(rhovec_st,ioption,popinit)
!
!  Calculate the steady-state density matrix for the current values of
!  the detunings. The actual calculation is done by subprograms of the
!  ldbl module. Using obe_steadystate makes it possible to work in the
!  weak probe approximation and also to work in a mode of operation better
!  adapted to calculations for a range of detunings.
!
!  ioption: If 0, 1 or 2, the calculation is done by ldbl_steadystate using
!           this value of ioption. If -1, the calculation is done by Method 1
!           of ldbl_steadystate but here with the a_mat and b arrays 
!           calculated in a way well suited to repeated calculations for
!           varying detunings.
!
!  The content of rhovec_st at entry is not used if ioption is 1 or -1.
!  It is used unless the optional argument popinit is specified if ioption is 2:
!  if popinit is not present, then the initial content of rhovec_st is passed
!  to and used by ldbl_steadystate, while if popinit is present then this
!  content if first replaced by the density matrix derived from the populations
!  specified in popinit. If ioption is 0, then rhovec_st / popinit may or may
!  not be used by ldbl_steadystate.
!
      implicit none
      integer, parameter :: neqs = nst*nst
      real(kd), dimension (neqs), intent(inout) :: rhovec_st
      double precision, dimension(nmn:nmx), optional, intent(in) :: popinit
      integer, intent(in) :: ioption
      real(kd), dimension (neqs,neqs) :: rhs_mat
      real(kd), dimension (neqs-1,neqs-1) :: a_mat
      real(kd), dimension (neqs-1) :: b
      integer :: ifail
!
      if(ioption.eq.0 .or. ioption.eq.1 .or. ioption.eq.2)then
!
         call obe_set_rhsmat(rhs_mat)
         if(ioption.eq.0 .or. ioption.eq.2) then
!  Overwrite the content of rhovec_st passed to obe_steadystate by
!  the density matrix derived from the array popinit, if popinit was passed to
!  obe_steadystate (it is an optional argument). Otherwise, don't change
!  rhovec_st before passing it to ldbl_steadystate.
            if(present(popinit)) call obe_init_rho(popinit,rhovec_st)
         endif
         call ldbl_steadystate(rhovec_st,ioption,rhs_mat)
!
      elseif(ioption.eq.-1)then
!
!  Calculate the arrays a_mat and b for the current values of the detunings:
         call obe_steadystate_setarrays(a_mat,b)
!  Solve the system and rearrange the result:
         call ldbl_steadystate_calc(a_mat,b,rhovec_st,ifail)
         if(ifail.ne.0)then
            print*,'obe_steadystate: Unsuccessful calculation.'
            call obe_error
         endif
!
      else
!
         print*,'obe_steadystate: Illegal value of ioption ',ioption
         call obe_error
!
      endif
!
      end subroutine obe_steadystate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
      subroutine obe_steadystate_setarrays(a_mat,b)
!
!  If fl_initarrays .eqv .false., calculate the arrays a_mat and b used by the  
!  subroutine steadystate.
!  If fl_initarrays .eqv .true, initialize this calculation before 
!  calculating a_mat and b.
!
!  This subroutine allocate a number of arrays used in the calculation at
!  its first call.
!
      implicit none
      integer, parameter :: neqs = nst*nst
      real(kd), dimension (neqs-1,neqs-1), intent(out) :: a_mat
      real(kd), dimension (neqs-1), intent(out) :: b
      real(kd), dimension (:,:), allocatable, save :: a_mat0
      real(kd), dimension (:), allocatable, save :: b0
      real(kd), dimension (:,:,:), allocatable, save :: a_mat1
      real(kd), dimension (:,:), allocatable, save :: b1
      real(kd), dimension (:,:), allocatable :: rhs_mat
      real(kd), dimension (:), allocatable :: Delta_store
      integer istat,k
      logical, save :: fl_first_call = .true.
!
      if(fl_first_call)then
!  Check that obe_setcsts has already been called; exit if it hasn't.
         if(fl_par_not_init)then
            print*,'obe_steadydstate_setarrays was called before the first'
            print*,'call to obe_setcsts'
            call obe_error
         endif
         allocate(a_mat0(neqs-1,neqs-1),b0(neqs-1),  &
                  a_mat1(neqs-1,neqs-1,nflds),b1(neqs-1,nflds),stat=istat)
         if(istat.gt.0)then
            print*,'Error at the allocation of the a or b arrays in'
            print*,'obe_steadystate_setarrays.'
            call obe_error
         endif
         fl_first_call = .false.
      endif
!
      if(fl_initarrays)then
!
!  Initialize...
!
         allocate(rhs_mat(neqs,neqs),Delta_store(nflds),stat=istat)
         if(istat.gt.0)then
            print*,'Error at the allocation of rhs_mat or Delta_store in'
            print*,'obe_steadystate_setarrays.'
            call obe_error
         endif
!
!  Store the current values of the detunings
         Delta_store = Delta
!
!  Calculate arrays with zero detunings
         Delta = 0.0_kd
         if(fl_mltpls)then
            call obe_calc_rhsmat(rhs_mat)
         else
            call obe_set_rhsmat(rhs_mat)
         endif
         call ldbl_steadystate_prep(rhs_mat,a_mat0,b0)
!
!  Calculate first order matrices, changing the detunings one by one
         do k=1,nflds
            Delta = 0.0_kd
            Delta(k) = 1.0_kd
            if(fl_mltpls)then
               call obe_calc_rhsmat(rhs_mat)
            else
               call obe_set_rhsmat(rhs_mat)
            endif
            call ldbl_steadystate_prep(rhs_mat,a_mat1(:,:,k),b1(:,k))
            a_mat1(:,:,k) = a_mat1(:,:,k) - a_mat0
            b1(:,k) = b1(:,k) - b0
         enddo
!
!  When done, reset the detunings to their original values, reset the
!  flag fl_initarrays, deallocate the arrays which are not needed anymore
!  and proceed with the calculation:
         Delta = Delta_store
         fl_initarrays = .false.
         deallocate(rhs_mat,Delta_store)
!
      endif
!
!  Calculate the arrays a_mat and b for the current values of the
!  detunings and return to the calling program:
!
      a_mat = a_mat0
      b = b0
      do k=1,nflds
         a_mat = a_mat + a_mat1(:,:,k)*Delta(k) 
         b = b + b1(:,k)*Delta(k)
      enddo
!
      end subroutine obe_steadystate_setarrays
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
      subroutine obe_check_iunit(iunit)
!
!  Check that the value of iunit makes sense for the unit number(s) to be
!  used by obe_Doppler_av_td_A or obe_Doppler_av_td_B.
!
      implicit none
      integer, intent(in) :: iunit
      logical :: fl
!
      if(iunit.gt.0)then
         if(iunit.eq.5 .or. iunit.eq.6)then
            print*,'Wrong unit number in obe_check_iunit ',iunit
            print*,'The unit numbers 5 and 6 are reserved for the standard'
            print*,'input and standard output files.'
            call obe_error
         endif
         inquire(unit=iunit,exist=fl)
         if(.not.fl)then
            print*,'obe_check_iunit: the unit iunit does not exist ',iunit
            call obe_error
         endif
         inquire(unit=iunit,opened=fl)
         if(.not.fl)then
            print*,'obe_check_iunit: the unit iunit is not connected ',iunit
            call obe_error
         endif
      else
         print*,'Wrong unit number in obe_check_iunit ',iunit
         call obe_error
      endif
!
      end subroutine obe_check_iunit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
      subroutine obe_set_Doppler(urms,npoints,irule,           &
                                 vmax,vmesh_in,vweight_in)
!
!  Set variables used by other parts of the module for Doppler-averaging:
!
!  urms is the r.m.s. speed in the light propagation direction, in m/s.
!  npoints is the number of abscissas in the integration over velocity
!  classes, vmesh_in are the abscissas and vweight_in are the corresponding
!  weights. The latter are multiplied by the Maxwell distribution
!  to form fMvweight.
!
      implicit none
      double precision, dimension (:), intent(in) :: vmesh_in,vweight_in
      double precision, intent(in) :: urms,vmax
      integer, intent(in) :: npoints,irule
      real(kd) :: fMax,sqrtpi
      double precision, allocatable :: x(:),w(:)
      integer :: imin,imax,istat,nv,nvsize
      character*8 :: rule
      optional :: vmesh_in,vweight_in
!
!  Preliminary check
      if(npoints.le.0)then
         print*,'obe_set_Doppler: the value of npoints is nonsensical. '
         print*,npoints
         call obe_error
      endif
!
      if(irule.eq.0)then
!
         if(.not.present(vmesh_in) .or. .not.present(vweight_in))then
            print*,'obe_set_Doppler: One or both of the optional arguments'
            print*,'is absent, although irule = 0.'
            call obe_error
         endif
         nvsize = size(vmesh_in)
         if(nvsize.ne.size(vweight_in))then
            print*,'obe_set_Doppler: vmesh_in and vweight_in have inconsistent'
            print*,'sizes. ',nvsize,size(vweight_in)
            call obe_error
         elseif(npoints.gt.nvsize)then
            print*,'obe_set_Doppler: the value of npoints seems incorrect.'
            print*,npoints,nvsize
            call obe_error
         elseif(npoints.lt.nvsize)then
            print*,'Warning: obe_set_Doppler is called with an oversized'
            print*,'vmesh_in array. ',npoints,nvsize
         endif
         imin = lbound(vmesh_in,1)
         imax = ubound(vmesh_in,1)
         if(lbound(vweight_in,1).ne.imin .or. ubound(vweight_in,1).ne.imax)then
            print*,'obe_set_Doppler detected an incorrect array bound.'
            call obe_error
         endif
!
      elseif(irule.le.4)then
!
         if(vmax.le.0.d0)then
            print*,'obe_set_Doppler: the value of vmax is non-sensical.'
            print*,vmax
         endif
         if(irule.eq.1)rule = 'Clenshaw'
         if(irule.eq.2)rule = 'Simpsons'
         if(irule.eq.3)rule = 'trapezdl'
         if(irule.eq.4)rule = 'Gaussian'
!  The arrays x and w are allocated by obe_quad_rule and deallocated
!  below.
         call obe_quad_rule(rule,npoints,x,w)
!  We integrate over v from -vmax to vmax, not from -1 to 1:
         x = vmax*x
         w = vmax*w
!
      else
!
         print*,'obe_set_Doppler: the value of irule is illegal. '
         print*,irule
         call obe_error
      endif
!
!  Copy this data into the relevant variables global within obe.
      urms_saved = real(urms,kd)
      n_v_values = npoints
      if(allocated(vmesh))deallocate(vmesh)
      if(allocated(fMvweight))deallocate(fMvweight)
      allocate(vmesh(n_v_values),fMvweight(n_v_values),stat=istat)
      if(istat.gt.0)then
         print*,'obe_set_Doppler: Error at the allocation of vmesh or '
         print*,'fMvweight.'
         call obe_error
      endif
!
!  We transcribe vmesh_in and vweight_in into vmesh and fMvweight,
!  multiplying the integration weights by the Maxwell distribution
      sqrtpi = sqrt(4.d0*atan(1.0_kd))
      if(irule.eq.0)then
         do nv = 1,n_v_values
            vmesh(nv) = real(vmesh_in(imin+nv-1),kd)
            fMax = exp(-(vmesh(nv)/urms_saved)**2)/(urms_saved*sqrtpi)
            fMvweight(nv) = fMax*real(vweight_in(imin+nv-1),kd)
         enddo
      else
         do nv = 1,n_v_values
            vmesh(nv) = real(x(nv),kd)
            fMax = exp(-(vmesh(nv)/urms_saved)**2)/(urms_saved*sqrtpi)
            fMvweight(nv) = fMax*real(w(nv),kd)
         enddo
         deallocate(x,w)
      endif
!
!  Set fl_Doppler_notinit to .false. to signal that urms_saved, n_v_values,
!  vmesh and fMvweight have been initialised.
      fl_Doppler_notinit = .false.
!
      end subroutine obe_set_Doppler

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
      subroutine obe_Doppler_av_td_A(irate,t1_in,t2_in,n_time_steps,imethod,   &
                                             rhovec_td_av,iunit,popinit)
!
!  Gets the Doppler-average density matrix over values of t ranging from
!  t1 to t2. Compared to a calculation using obe_Doppler_av_td_B, a calculation
!  using obe_Doppler_av_td_A is normally faster but more demanding in memory.
!
!  At entry, t1 should be the initial time, t2 the final (largest) time,
!  n_time_steps the number of time steps at which the density matrix needs
!  to be calculated, imethod the index of the method to be used to solve
!  the ODEs, and iunit the unit nunber of the file used to store the
!  Doppler-average density matrix at the end of the calculation (nothing
!  is written on file if iunit.eq.0, and special output files are used
!  if iunit.eq.-1). It is assumed that the density matrix
!  at t1 is rhovec_td_av for all the velocity classes, unless the optional
!  argument popinit is specified. In the latter case, it is assumed that the
!  coherences are all initially zero and the populations are those specified
!  by the content of popinit.
!
!  On return, rhovec_td_av contains the density matrix at t2.
!
!  This subroutine differs from obe_Doppler_av_td_B primarily in that 
!  the integration over velocity classes is done for all time steps at once
!  rather than time step by by time step. Moreover, obe_Doppler_av_td_A creates
!  a potentially very large array in the ldblstore module, whilst
!  obe_Doppler_av_td_B does not.
!
      use ldblstore
      implicit none
      integer, parameter :: neqs=nst*nst
      double precision, intent(in) :: t1_in,t2_in
      double precision, dimension(nmn:nmx), optional, intent(in) :: popinit
      integer, intent(in) :: imethod,irate,iunit,n_time_steps
      real(kd), dimension(nst*nst), intent(inout) :: rhovec_td_av
      double precision :: t
      integer :: iost,istat,nt,i,j,m,mre,mim
      logical :: fl_store
      character(len=11) :: file_format
!
!  Check that the necessary information about the Doppler averaging
!  has been fed to obe through a prior call to obe_set_Doppler.
      if(fl_Doppler_notinit)then
         print*,'obe_Doppler_av_td_A: obe_set_Doppler has not been'
         print*,'called prior to a call to this subroutine'
         call obe_error
      endif
!
!  Check that the value of iunit makes sense, and find out whether the
!  corresponding file (if any) is meant to be formatted or unformatted.
!  Check the value of n_time_steps, too.
!
      if(iunit.gt.0)then
         call obe_check_iunit(iunit)
         inquire(unit=iunit,form=file_format)
      elseif(iunit.lt.-1)then
         print*,'obe_Doppler_av_td_A: illegal value of iunit. ',iunit
         call obe_error
      endif
      if(n_time_steps.lt.1)then
         print*,'obe_Doppler_av_td_A has been called with a non-sensical'
         print*,'value of n_time_steps. ',n_time_steps
         call obe_error
      endif
!
      irate_Dopp = irate
      t1_Dopp = real(t1_in,kd)
      t2_Dopp = real(t2_in,kd)
      n_time_steps_Dopp = n_time_steps
      imethod_Dopp = imethod
      if(present(popinit))then
         call obe_init_rho(popinit,rhovec_init_Dopp)
      else
         rhovec_init_Dopp = rhovec_td_av
      endif
      istore_Dopp = -1 ! Use the istore.eq.-1 option of ldbl_tdin.
      iunit_Dopp = -1  ! No file writing within ldbl_tdint here
!
!  fl_Doppl_td .eqv. .true. is the signal for obe_Doppler_integrand to call
!  obe_Doppler_timedep rather than obe_Doppler_steadystate to calculate the
!  integrand of the integral over the velocity classes. fl_Doppl_td is a
!  global variable within the obe module.
      fl_Doppl_td = .true.
!
!  The arrays rho_allt belongs to the module ldblstore.
      if(allocated(rho_allt))deallocate(rho_allt)
      allocate(rho_allt(0:neqs,0:n_time_steps),stat=istat)
      if(istat.gt.0)then
         print*,'obe_Doppler_av_td_A: Error at the allocation of rho_allt.'
         call obe_error
      endif
!
!  Integrate over the velocity classes (rho_allt is calculated within the
!  ldble module, as driven indirectly by obe_Doppler_vint) 
      rho_allt = 0.0_kd
      fl_store = .true.
      call obe_Doppler_vint(fl_store,rhovec_td_av)
!
!  Write the Doppler-average density matrix on the unit iunit if iunit > 0
!  or on the unit(s) specified in a prior call to obe_setoutputfiles if
!  iunit = 0. The first element of each column of rho_allt is the
!  corresponding time, t.
      if(iunit.gt.0)then
         rewind iunit
         do nt = 0,n_time_steps
            if(file_format.eq.'UNFORMATTED')then
               write(iunit,iostat=iost)rho_allt(0:neqs,nt)
            else
               write(iunit,*,iostat=iost)rho_allt(0:neqs,nt)
            endif
            if(iost.ne.0)then
               print*,'obe_Doppler_av_td_A: write returns with iostat = ',iost
               call obe_error
            endif
         enddo
      elseif(iunit.eq.-1)then
         if(fl_outputfiles)then
            do nt = 0,n_time_steps
               t = rho_allt(0,nt)
               do j = 1,nst
                  do i = 1,nst
                     if(iunits_arr(i,j).gt.0)then
                        if(i.eq.j)then
                           call ldbl_pop_index(i,m)
                           write(iunits_arr(i,j),1501)t,rho_allt(m,nt)
                        elseif(i.lt.j)then
                           call ldbl_coher_index(i,j,mre,mim)
                           write(iunits_arr(i,j),1502)t,rho_allt(mre,nt), &
                                                        rho_allt(mim,nt)
                        else
                           call ldbl_coher_index(j,i,mre,mim)
                           write(iunits_arr(i,j),1502)t,rho_allt(mre,nt), &
                                                       -rho_allt(mim,nt)
                        endif
                     endif
                  enddo
               enddo
            enddo
         else
            print*,'obe_Doppler_av_td_A: The value of iunit is -1 but the unit'
            print*,'numbers of the output files have not been provided.'
            call obe_error
         endif
      endif
      deallocate(rho_allt)
!
 1501 format(1x,2(1pe12.5,2x))
 1502 format(1x,3(1pe12.5,2x))
!
!  and return...
!
      end subroutine obe_Doppler_av_td_A
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine obe_Doppler_vint(fl_store,rhovec_av)
!
!  Integrate the density matrix over the velocity classes.
!
!  fl_store: a logical variable. The variable factor, defined
!  within the module ldblstore, is given a value by obe_Doppler_vint
!  if fl_store is .true. Otherwise this variable is ignored.
!  (fl_store is normally .true. only if obe_Doppler_vint is
!  called from obe_Doppler_av_td_A.)
!
!  The integration is done within this subroutine for rhovec_av
!  and, if required, within the ldbl module for rho_allt.
!
!  The arrays vmesh and fMvweights, containing, respectively, the
!  abscissas of the numerical quadrature and the corresponding weights
!  multiplied by the Maxwell distribution, are global variables within obe.
!
      use ldblstore
      implicit none
      integer, parameter :: neqs=nst*nst
      real(kd), dimension(nst*nst), intent(out) :: rhovec_av
      logical, intent(in) :: fl_store
      real(kd), dimension(neqs) :: rhovec
      integer :: nv
!
      rhovec_av = 0.0_kd
      do nv=1,n_v_values
         if(fl_store)factor = fMvweight(nv)
         call obe_Doppler_integrand(vmesh(nv),rhovec)
         rhovec_av = rhovec_av + fMvweight(nv)*rhovec
      enddo
!
      end subroutine obe_Doppler_vint

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
      subroutine obe_Doppler_av_td_B(t1_in,t2_in,n_time_steps,            &
                                       rhovec_td_av,iunit,popinit)
!
!  Average the density matrix at time t2 over the velocity distribution.
!
!  The time-integration itself is actually performed by ldbl_tdint, which
!  is accessed through obe_Doppler_timedep. The dummy variables
!  t1, t2, n_time_steps, imethod and iunit have the same meanings as in the
!  subroutine ldbl_tdint. Here, the interval [t1,t2] is divided into
!  n_time_steps of duration tstep, and Doppl_int is called to obtain the
!  Doppler-averaged density matrix first at t1 + tstep, then at t1 + 2 tstep,
!  then at t1 + 3 tstep, etc., each time starting at t1. The values of the
!  parameters the subroutine obe_Doppler_timedep must use are stored in the
!  variables t1_Dopp, t2_Dopp, n_time_steps_Dopp, imethod_Dopp and iunit_Dopp,
!  which are global within the obe module. The subroutine obe_Doppler_timedep
!  uses obe_tdint to integrate the optical Bloch equations using the
!  integration method defined by the value of the variable imethod.
!  As the Runge-Kutta method is not suitable for the large integration
!  intervals normally involved in here, imethod must be either 2 or 3.
!
!  rhovec_td_av, popinit:  As for obe_Doppler_av_td_A.
!
      implicit none
      integer, parameter :: neqs=nst*nst
      double precision, intent(in) :: t1_in,t2_in
      integer, intent(in) :: iunit,n_time_steps
      integer :: imethod
      real(kd), dimension(neqs), intent(inout) :: rhovec_td_av
      double precision, dimension(nmn:nmx), optional, intent(in) :: popinit
      real(kd) :: t,t1,t2,tstep
      integer :: i,iost,j,m,mre,mim,n
      logical :: fl_store
      character(len=11) :: file_format
!
!  The only method option possible in the current state of development of
!  the module is imethod = 3...
      imethod = 3
!
!  Check that the necessary information about the Doppler averaging
!  has been fed to obe through a previous call to obe_set_Doppler.
      if(fl_Doppler_notinit)then
         print*,'obe_Doppler_av_td_B: obe_set_Doppler has not been'
         print*,'called prior to a call to this subroutine'
         call obe_error
      endif
!
!  Calculations within the rate equation approximation are not possible
!  with this subroutine. Hence irate_Dopp is necessarily not equal to 1.
      irate_Dopp = 0
!
!  Check that the value of iunit makes sense, and find out whether the
!  corresponding file (if any) is meant to be formatted or unformatted.
!
      if(iunit.gt.0)then
         call obe_check_iunit(iunit)
         inquire(unit=iunit,form=file_format)
      elseif(iunit.lt.-1)then
         print*,'obe_Doppler_av_td_B: illegal value of iunit. ',iunit
         call obe_error
      endif
!
!  Check that the value of imethod is appropriate. Check the value of
!  n_time_steps, too.
!
      if(imethod.ne.2 .and. imethod.ne.3)then
         print*,'obe_Doppler_av_td_B: imethod must be 2 or 3 ',imethod
         call obe_error
      endif
      if(n_time_steps.lt.1)then
         print*,'obe_Doppler_av_td_B has been called with a non-sensical'
         print*,'value of n_time_steps. ',n_time_steps
         call obe_error
      endif
!
!  fl_Doppl_td .eqv. .true. is the signal for obe_Doppler_integrand to call
!  obe_Doppler_timedep rather than obe_Doppler_steadystate to calculate the
!  integrand of the integral over the velocity classes. fl_Doppl_td is a
!  global variable within the obe module.
      fl_Doppl_td = .true.
!
!  Convert t1_in and t2_in from double precision to real, if necessary:
      t1 = real(t1_in,kd)
      t2 = real(t2_in,kd)
!
!  Calculate the time step for the required number of steps:
      tstep = (t2-t1)/real(n_time_steps,kd)
!
!  Initial values (at t1)
      n = 0
      t = t1
      if(present(popinit))then
         call obe_init_rho(popinit,rhovec_init_Dopp)
      else
         rhovec_init_Dopp = rhovec_td_av
      endif
!
!  Loop over time, starting at t1.
      do n=0,n_time_steps
         t = t1 + n*tstep
         t1_Dopp = t1
         t2_Dopp = t
         n_time_steps_Dopp =  1 ! ldbl_tdint should go from t1 to t in one step
         imethod_Dopp = imethod
         istore_Dopp = 0  ! no use of ldbstore here...
         iunit_Dopp = -1  ! no writing on file within ldbl_tdint...
!
         fl_store = .false.
         call obe_Doppler_vint(fl_store,rhovec_td_av)
!
!  ... and write the result on file(s)
         if(iunit.gt.0)then
            if(file_format.eq.'UNFORMATTED')then
               write(iunit,iostat=iost)t,(rhovec_td_av(j),j=1,neqs)
            else
               write(iunit,*,iostat=iost)t,(rhovec_td_av(j),j=1,neqs)
            endif
            if(iost.ne.0)then
               print*,'obe_Doppler_av_td_B: write returns with iostat = ', &
                  iost,n
               call obe_error
            endif
         elseif(iunit.eq.-1)then
            if(fl_outputfiles)then
               do j = 1,nst
                  do i = 1,nst
                     if(iunits_arr(i,j).gt.0)then
                        if(i.eq.j)then
                           call ldbl_pop_index(i,m)
                           write(iunits_arr(i,j),1501)t,rhovec_td_av(m)
                        elseif(i.lt.j)then
                           call ldbl_coher_index(i,j,mre,mim)
                           write(iunits_arr(i,j),1502)t,rhovec_td_av(mre), &
                                                        rhovec_td_av(mim)
                        else
                           call ldbl_coher_index(j,i,mre,mim)
                           write(iunits_arr(i,j),1502)t,rhovec_td_av(mre), &
                                                       -rhovec_td_av(mim)
                        endif
                     endif
                  enddo
               enddo
            else
               print*,'obe_Doppler_av_td_B: The value of iunit is -1 but the unit'
               print*,'numbers of the output files have not been provided.'
               call obe_error
            endif
         endif
      enddo
!
 1501 format(1x,2(1pe12.5,2x))
 1502 format(1x,3(1pe12.5,2x))
!
      end subroutine obe_Doppler_av_td_B

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine obe_Doppler_timedep(v,rho_td_vec)
!
!  Calculate the time-dependent density matrix for Doppler-shifted detunings.
!
!  See ldbl_tdint for the method.
!
      implicit none
      integer, parameter :: neqs=nst*nst
      real(kd), intent(in) :: v
      real(kd), dimension(neqs), intent(out) :: rho_td_vec
      real(kd), dimension(neqs,neqs) :: rhs_mat
      real(kd), dimension(neqs,neqs), save :: rhs_mat0,rhs_mat1
      real(kd) :: t1,t2
      integer :: irate
!
!  Initialize the arrays if required
      if(fl_tdDopp_init)then
!  Calculation for zero velocity:
         call obe_set_rhsmat(rhs_mat0,0.0_kd)
!  Calculation for a velocity of 1 m/s
         call obe_set_rhsmat(rhs_mat1,1.0_kd)
!  Put the contribution to rhs_mat proportional to v into rhs_mat1:
         rhs_mat1 = rhs_mat1 - rhs_mat0
!  Reset fl_tdDopp_init as the arrays are now initialized
         fl_tdDopp_init = .false.
      endif
!
!  Calculate rhs_mat for the velocity v:
      rhs_mat = rhs_mat0 + v*rhs_mat1
!
!  Solve the system, using the density matrix saved in rhovec_init_Dopp
!  by obe_Doppler_av_td_A or B as initial value.
      irate = irate_Dopp
      rho_td_vec = rhovec_init_Dopp
      t1 = t1_Dopp
      t2 = t2_Dopp
      call ldbl_tdint(irate,t1,t2,n_time_steps_Dopp,imethod_Dopp,rho_td_vec,  &
                        istore_Dopp,iunit_Dopp,rhs_mat,isteady=0)
!
      end subroutine obe_Doppler_timedep
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine obe_Doppler_av_st_numerical(rhovec_st_av,ioption,popinit)
!
!  Average the steady-state density matrix over the velocity distribution.
!  The calculation is done by numerical integration.
!
!  ioption has the same meaning as for the subroutine obe_steadystate; however,
!  here the possible values of this variable are restricted to 1 or 2.
!
      implicit none
      integer, parameter :: neqs = nst*nst
      real(kd), dimension (nst*nst), intent(inout) :: rhovec_st_av
      double precision, dimension (nmn:nmx), optional, intent(in) :: popinit
      integer, intent(in) :: ioption
      integer :: istat
      logical :: fl_store
!
!  Check that the value of ioption is suitable. If it is, save it for later
!  use by obe_Doppler_steadystate.
      if(ioption.ne.1 .and. ioption.ne.2)then
         print*,'obe_Doppler_av_st_numerical: illegal value of ioption.',ioption
         call obe_error
      else
         ioption_Dopp = ioption
      endif
!
!  Allocate the arrays used by obe_Doppler_steadystate and prepare them
!  as necessary:
      if(ioption.eq.1)then
         if(allocated(Delta_store_Dopp))deallocate(Delta_store_Dopp)
         allocate(Delta_store_Dopp(nflds),stat=istat)
         if(istat.gt.0)then
            print*,'obe_Doppler_av_st_numerical: error at the allocation of'
            print*,'Delta_store_Dopp.'
            call obe_error
         endif
      else
         if(allocated(rhs0v_Dopp))deallocate(rhs0v_Dopp)
         if(allocated(rhs1v_Dopp))deallocate(rhs1v_Dopp)
         if(allocated(rhsv_Dopp))deallocate(rhsv_Dopp)
         allocate(rhs0v_Dopp(neqs,neqs),rhs1v_Dopp(neqs,neqs), &
                  rhsv_Dopp(neqs,neqs),stat=istat)
         if(istat.gt.0)then
            print*,'obe_Doppler_av_st_numerical: error at the allocation of'
            print*,'rhs0v_Dopp, rhs1v_Dopp and rhsv_Dopp.'
            call obe_error
         endif
         call obe_set_rhsmat(rhs0v_Dopp)  
         call obe_set_rhsmat(rhs1v_Dopp,1.0_kd)  
         rhs1v_Dopp = rhs1v_Dopp - rhs0v_Dopp
         if(present(popinit))then
            call obe_init_rho(popinit,rhovec_init_Dopp)
         else
            rhovec_init_Dopp = rhovec_st_av
         endif
      endif
!
!  fl_Doppl_td .eqv. .false. is the signal for obe_Doppler_integrand to call
!  obe_Doppler_steadystate to calculate the integrand of the integral over
!  the velocity classes, not obe_Doppler_timedep. fl_Doppl_td is a global
!  variable within the obe module.
      fl_Doppl_td = .false.
!
!  Integrate...
      fl_store = .false.
      call obe_Doppler_vint(fl_store,rhovec_st_av)
!
!  deallocate the arrays not needed anymore...
      if(ioption.eq.1)then
         deallocate(Delta_store_Dopp)
      else
         deallocate(rhs0v_Dopp,rhs1v_Dopp,rhsv_Dopp)
      endif
!
!  ... and return.
!
      end subroutine obe_Doppler_av_st_numerical

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine obe_Doppler_steadystate(v,rho_st_vec)
!
!  Calculate the steady-state density matrix for Doppler-shifted detunings.
!
!  See obe_steadystate for the method.
!
      implicit none
      integer, parameter :: neqs = nst*nst
      real(kd), dimension (neqs-1,neqs-1) :: a_mat
      real(kd), dimension (neqs-1,neqs-1), save :: a_mat0,a_mat1
      real(kd), dimension (neqs-1) :: b
      real(kd), dimension (neqs-1), save :: b0,b1
      real(kd), intent(in) :: v
      real(kd), dimension (neqs), intent(out) :: rho_st_vec
      real(kd), dimension(:), allocatable :: Delta_store
      integer :: ifail
      logical, save :: fl_first_call = .true.
!
      if(fl_first_call)then
!  Check that obe_setcsts has already been called; exit if hasn't.
         if(fl_par_not_init)then
            print*,'obe_Doppler_steadystate was called before the first'
            print*,'call to obe_setcsts'
            call obe_error
         endif
         fl_first_call = .false.
      endif
!
!  Initialize the arrays if required
      if(fl_stDopp_init .and. ioption_Dopp.eq.1)then
!
!  Calculation without Doppler shift
         call obe_steadystate_setarrays(a_mat0,b0)
!  Calculation for a Doppler shift of 1 m/s (obe_Doppler_shift shifts the
!  detunings; the original values are saved into Delta_store_Dopp during this
!  operation)
         Delta_store_Dopp = Delta
         call obe_Doppler_shift(1.0_kd) 
         call obe_steadystate_setarrays(a_mat1,b1)
         Delta = Delta_store_Dopp
!  Put the contributions to a_mat and b proportional to v into a_mat1 and b1:
         a_mat1 = a_mat1 - a_mat0
         b1 = b1 - b0
!  Reset fl_stDopp_init as the arrays are now initialized
         fl_stDopp_init = .false.
      endif
!
      if(ioption_Dopp.eq.1)then
!  Calculate the arrays a_mat and b for the velocity v:
         a_mat = a_mat0 + v*a_mat1
         b = b0 + v*b1
!  Solve the system. In case where the calculation of the steady state
!  would be unsucessful, keep going but issue an error message
!  and set the elements of the density matrix to nonsensical values.
         call ldbl_steadystate_calc(a_mat,b,rho_st_vec,ifail)
         if(ifail.ne.0)then
            print*,'obe_Doppler_steadystate: calculation unsuccessful! ', v
            rho_st_vec = -999.0_kd
         endif
      else
!  Calculate the matrix representing the right-hand-side
         rhsv_Dopp = rhs0v_Dopp + v*rhs1v_Dopp
!  Solve the system
         rho_st_vec = rhovec_init_Dopp
         call ldbl_steadystate(rho_st_vec,ioption_Dopp,rhsv_Dopp)
      endif
!
      end subroutine obe_Doppler_steadystate
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
      subroutine obe_Doppler_integrand(v,rho_vec)
!
!  The integrand of the integral over velocities, for Doppler averaging in 1D,
!  without the Maxwell distribution.
!
      implicit none
      real(kd), intent (in) :: v
      real(kd), dimension (nst*nst), intent(out) :: rho_vec
!
      if(fl_Doppl_td)then
         call obe_Doppler_timedep(v,rho_vec)
      else
         call obe_Doppler_steadystate(v,rho_vec)
      endif
!
      end subroutine obe_Doppler_integrand
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine obe_Doppler_shift(v)
!
!  Shift the detunings by the Doppler shift for a 1D velocity v.
!  The atomic mass and the wavelenghts of the fields are meant to be passed
!  to obe through obe_setcsts and obe_setfields. The atomic mass must be
!  specified in unified atomic mass units (u), the wavelenghts in nm, and
!  v in m/s. As the detunings are frequencies, not angular frequencies, the
!  shifts are plus or minus kv/(2 pi), where k is the wave number.
!  Hence, these frequency shifts are plus or minus v/wavelength.
!
!  Given that the detunings are in MHz units, the wavelengths in nm
!  and v in m/s, the ratios v/wavelength must be multiplied by a factor
!  of 1d9/1d6.
!
      implicit none
      real(kd), parameter :: conversion_factor = 1.e9_kd/1.e6_kd
      real(kd), intent(in) :: v
      integer :: k
!
!  Add the shift to the detuning if the field propagates in the negative
!  direction and subtract it from the detuning if the field propagates in
!  the positive direction.
      do k = 1,nflds
         if(idir_st(k).ne.1 .and. idir_st(k).ne.-1)then
            print*,'obe_Doppler_shift: Wrong idir attribute for field ',k
            call obe_error
         endif
         if(wavel_nm(k).eq.0.0_kd)then
            print*,'obe_Doppler_shift found a zero wavelength for field ',k
            call obe_error
         endif
         Delta(k) = Delta(k) - real(idir_st(k),kd)*conversion_factor*     &
                                                          v/wavel_nm(k)
      enddo
!
      end subroutine obe_Doppler_shift

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine obe_susceptibility(field,rhovec,density,chi,refr_index,alpha)
!
!  Calculates the complex susceptibility and, optionally, the corresponding
!  refractive index (refr_index) and absorption coefficient (alpha), given the
!  coherences, for the field defined by the variable field and the atomic
!  density defined by the variable density.
!
!  refr_index and alpha are optional arguments.
!
!  Warning: It is essential, in calculations using the iRabi.eq.1 option of
!  obe_setcsts, that the dip_mom and amplitude attributes of the
!  variable field are consistent with the Rabi frequencies used to calculate
!  the density matrix.
!
      implicit none
      integer, parameter :: neqs = nst*nst
!
      type(obecfield), intent(in) :: field
      real(kd), dimension(neqs), intent(in) :: rhovec
      double precision, intent(in) :: density
      double complex, intent(out) :: chi
      double precision, intent(out), optional :: alpha,refr_index
      double complex :: sumrho
      double precision, save :: pi
      double precision :: ak, beta
      integer :: i,j,idxre,idxim,n
      logical, save :: fl_init = .true.
!
      if(fl_init)then
         pi = 4.0d0 * datan(1.0d0)
         fl_init = .false.
      endif
!
      sumrho = (0.0_kd,0.0_kd)
      do i = nmn,nmx
         do j = nmn,nmx
            if(field%dip_mom(i,j).ne.(0.d0,0.d0))then
               if(field%dip_mom(j,i).ne.(0.d0,0.d0))then
                  print*,'obe_susceptibility has detected an error:'
                  print*,'the dip_mom component of field is incorrectly'
                  print*,'formatted.'
                  print*,i,j,field%dip_mom(i,j),field%dip_mom(j,i)
                  call obe_error
               endif
               if(j.gt.i)then
                  call obe_coher_index(i,j,idxre,idxim)
                  sumrho = sumrho + dconjg(field%dip_mom(i,j)) *               &
                                   cmplx(rhovec(idxre),rhovec(idxim),kd)
               elseif(j.lt.i)then
                  call obe_coher_index(j,i,idxre,idxim)
                  sumrho = sumrho + dconjg(field%dip_mom(i,j)) *               &
                                   cmplx(rhovec(idxre),-rhovec(idxim),kd)
               else
                  print*,'obe_susceptibility: error in field%dip_mom. ',i
                  call obe_error
               endif
            endif
         enddo
      enddo
!
!  The following piece of code mixes double precision and real*4 precision
!  if kd refers to real*4 precision rather than real*8 precision, but
!  the necessary conversion to double precision will be done automatically.
      chi = 2.0d0 * density * sumrho / (epsilon0 * field%amplitude)
      if(present(refr_index))then
         refr_index = dreal(sqrt((1.0d0,0.0d0)+chi))
      endif
      if(present(alpha))then
         beta = dimag(sqrt((1.0d0,0.0d0)+chi))  ! "extinction coefficient"
         ak = 2.0d0 * pi / (field%wavelength*1.0d-9) ! vacuum wave number
         alpha = 2.0d0 * ak * beta
      endif
!
      end subroutine obe_susceptibility

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine obe_weakfield(energ_f,gammatot_f,field,popinit,        &
                       detunings,density,refr_index,alpha,iDoppler,urms)
!
!  energ_f, gammatot_f, detunings: as frequencies in MHz, not
!                                  angular frequencies!
!  iDoppler.eq.1 : with Doppler averaging. Otherwise iDoppler must be 0.
!
!  gammatot_f refers to what is denoted gamma^tot_{ij} in Appendix E of the
!  notes. I.e., the total coherence decay rates, including spontaneous
!  decay, laser linewidth and other dephasing rates. E.g., without the
!  laser linewidth, gammatot_f = Gamma/2 + Gamma_coll/2, where Gamma is the
!  natural width of the upper state (assuming the lower state is stable) and
!  Gamma_coll is, typically, the collisional width (FWHM).
!
!  It is assumed that field%dip_mom(i,j) can be non-zero only if state j is
!  lower in energy than state i, which is the normal way of configuring
!  this component. Accordingly, it is assumed that state i can decay
!  to state j if dip_mom(i,j) is non-zero, and therefore that the
!  population of state i is zero in the weak probe limit. 
!
!  popinit: initial populations (typically, the population in the upper state
!  of each of the transitions must be 0).
!
!  density: must be expressed in units of 1/m^3.
!
      implicit none
      type(obecfield), intent(in) :: field
      double precision, dimension(nmn:nmx), intent(in) :: energ_f, popinit
      double precision, dimension(nmn:nmx,nmn:nmx), intent(in) :: gammatot_f
      double precision, dimension(:), intent(in) :: detunings
      double precision, dimension(:), intent(out) :: alpha,refr_index
      double precision, intent(in) :: density,urms
      integer, intent(in) :: iDoppler
      optional :: urms
!
      double complex :: chi,coher,w
      double precision, dimension(:), allocatable :: fact_tr,gamma_tr,   &
                                                     delta_offset
      double precision :: x,y,u,v
      double precision :: ak,beta,delta,delta_trans,pi
      integer :: i,id,imin,imax,istat,j,n,ndet,ntrans
      logical :: fl_Doppler
!
      pi = 4.d0 * datan(1.d0)
!
!  Check the sizes of the arrays provided to obe_weakfields
      imin = lbound(detunings,1)
      imax = ubound(detunings,1)
      if(lbound(refr_index,1).ne.imin .or. ubound(refr_index,1).ne.imax .or.       &
         lbound(alpha,     1).ne.imin .or. ubound(alpha,     1).ne.imax)then
         print*,'obe_weakfield: Inconsistent array bounds.'
         call obe_error
      endif
      ndet = imax - imin + 1
!
!  Check that the value of iDoppler makes sense and that urms 
!  is specified if the calculation requires Doppler-averaging
      if(iDoppler.eq.1)then
         fl_Doppler = .true.
         if(.not.present(urms))then
            print*,'obe_weakfield: urms must be specified'
            print*,'for a calculation with Doppler averaging.'
            call obe_error
         endif
      elseif(iDoppler.eq.0)then
         fl_Doppler = .false.
      else
         print*,'obe_weakfield: illegal value of iDoppler ',iDoppler
         call obe_error
      endif
!
!  Find out how many pairs of states need to be taken into account.
!  We skips pairs for which the population in the lower state is zero.
      ntrans = 0
      do j = nmn,nmx
         do i = nmn,nmx
            if(field%dip_mom(i,j).ne.(0.d0,0.d0))then
               if(field%dip_mom(j,i).ne.(0.d0,0.d0))then
                  print*,'obe_weakfield has detected an error:'
                  print*,'the dip_mom component of field is incorrectly'
                  print*,'formatted.'
                  print*,i,j,field%dip_mom(i,j),field%dip_mom(j,i)
                  call obe_error
               endif
               if(popinit(j).gt.0.d0)ntrans = ntrans + 1
            endif
         enddo
      enddo
      if(ntrans.eq.0)then
         print*,'obe_weakfield: Error, the field does not couple any of'
         print*,'the states...'
         call obe_error
      endif
!
!  Store the details of the transitions.
      if(allocated(fact_tr))deallocate(fact_tr)
      if(allocated(gamma_tr))deallocate(gamma_tr)
      if(allocated(delta_offset))deallocate(delta_offset)
      allocate(fact_tr(ntrans),gamma_tr(ntrans),delta_offset(ntrans),     &
                                                             stat=istat)
      if(istat.gt.0)then
         print*,'obe_weakfield: error at the allocation of the arrays.'
         call obe_error
      endif
      n = 0
      do j = nmn,nmx
         do i = nmn,nmx
            if(field%dip_mom(i,j).ne.(0.d0,0.d0) .and. popinit(j).gt.0.d0)then
               n = n + 1
               if(n.gt.ntrans)then
                  print*,'obe_weakfield: the program is broken (1).'
                  call obe_error
               endif
               if(popinit(i).gt.0.d0)then
                  print*,'obe_weakfield: the upper states populations'
                  print*,'must be zero.', i,popinit(i)
                  call obe_error
               endif
!  Homogeneous broadening for this pair of states:
               gamma_tr(n) = 2.d0*pi*gammatot_f(i,j)*1.d6
               delta_offset(n) = 2.d0*pi*(energ_f(i) - energ_f(j) )*1.d6
!  Factor multiplying the coherence for the j to i transition in the
!  susceptibility. The variables hbar and epsilon are defined in the
!  module obe_constants.
               fact_tr(n) = cdabs(field%dip_mom(i,j))**2 *         &
                          (popinit(j)*density) / (hbar*epsilon0)
            endif
          enddo
      enddo
      if(n.ne.ntrans)then
         print*,'obe_weakfield: the program is broken (2).'
         call obe_error
      endif
!
      ak = 2.0d0 * pi / (field%wavelength*1.0d-9) ! vacuum wave number
!
!  Loop over the detunings specified in the array detunings:
!
      do id = 1,ndet 
         delta = 2.d0*pi*detunings(imin+id-1)*1.d6
!  Sum the susceptibilities
         chi = (0.d0,0.d0)
         do n = 1,ntrans
            delta_trans = delta - delta_offset(n)
!
            if(fl_Doppler)then
!  Coherence with Doppler averaging:
               x = delta_trans/(ak*urms)
               y = gamma_tr(n)/(ak*urms)
               call obe_Faddeeva(x,y,u,v)
               w = dcmplx(u,v)
               coher = (0.d0,1.d0)*sqrt(pi)*w/(ak*urms)
            else
!  Coherence without Doppler averaging:
               coher = (0.d0,1.d0)/(gamma_tr(n)-(0.d0,1.d0)*delta_trans)
            endif
!
            chi = chi + fact_tr(n) * coher 
         enddo
!
         refr_index(id) = dreal(sqrt((1.d0,0.d0)+chi))
         beta = dimag(sqrt((1.d0,0.d0)+chi))  ! "extinction coefficient"
         alpha(id) = 2.d0 * ak * beta
!
      enddo
!
      deallocate(fact_tr,gamma_tr,delta_offset)
!
      end subroutine obe_weakfield

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine obe_2state(Omega_p_0,Delta_p,Gamma_AB,frqwidth,            &
                            other_dephas_AB,iweakprb,iDoppler,rhoAA,rhoAB,  &
                            wvl_p,urms)
!
!  Two-state ladder system, single field.
!  Weak probe limit with or without Doppler averaging.
!
      implicit none
      double complex, intent(in) :: Omega_p_0
      double precision, intent(in) :: Delta_p,Gamma_AB,other_dephas_AB
      double precision, intent(in) :: frqwidth,wvl_p,urms
      integer, intent(in) :: iDoppler,iweakprb
      double complex, intent(out) :: rhoAB
      double precision, intent(out) :: rhoAA
      optional :: wvl_p,urms
!
!  Check whether the calculation is with or without Doppler averaging
!  and proceed accordingly.
!
      if(iDoppler.eq.1 .and. iweakprb.eq.1)then
         call obe_2state_weak_withav(rhoAA,rhoAB)
      elseif(iDoppler.eq.1 .and. iweakprb.eq.0)then
         call obe_2state_withav(rhoAA,rhoAB)
      elseif(iDoppler.eq.0)then
         call obe_2state_noav(rhoAA,rhoAB)
      else
         print*,'obe_2state: illegal value of iDoppler or iweakprb '
         print*,iDoppler,iweakprb
         call obe_error
      endif
!
      contains
  
!!!!!!!!!!!!!!!!!!

      subroutine obe_2state_noav(rhoAA,rhoAB)
!
!  Calculation without Doppler averaging
!
      implicit none
      double complex, intent(out) :: rhoAB
      double precision, intent(out) :: rhoAA
      double complex :: rhoBA,a
      double precision :: sqr,gamtot
      double complex :: oneoveretc
!
      if(iweakprb.eq.1)then
         a = Gamma_AB/2.d0 + other_dephas_AB + frqwidth/2.d0 -        &
                                                (0.d0,1.d0)*Delta_p
         rhoAA = 1.d0
         rhoBA = ((0.d0,1.d0)*Omega_p_0/2.d0)/a
         rhoAB = dconjg(rhoBA)
      elseif(iweakprb.eq.0)then
         gamtot = Gamma_AB/2.d0 + other_dephas_AB + frqwidth/2.d0
         sqr = sqrt(gamtot**2+cdabs(Omega_p_0)**2*gamtot/Gamma_AB)
         oneoveretc = (0.d0,1.d0)/dcmplx(Delta_p,sqr)
         rhoAA = 1.d0 - cdabs(Omega_p_0)**2*gamtot*                            &
                    dreal(oneoveretc)/(2.d0*sqr*Gamma_AB)
         rhoAB = -(0.d0,1.d0)*dconjg(Omega_p_0)/2.d0* &
                    dcmplx((gamtot/sqr)*dreal(oneoveretc),-dimag(oneoveretc))
      else
         print*,'obe_2state_noav: illegal value of iweakprb '
         print*,iweakprb
         call obe_error
      endif
!
      end subroutine obe_2state_noav
  
!!!!!!!!!!!!!!!!!!

      subroutine obe_2state_weak_withav(rhoAA,rhoAB)
!
!  Calculation in the weak field limit with Doppler averaging
!
      implicit none
      double complex, intent(out) :: rhoAB
      double precision, intent(out) :: rhoAA
      double complex :: rhoBA
      double complex :: Omega_p,eta,d,w,rp
      double complex :: complexint
      double precision :: gammaBtot_f,gammaCtot_f
      double precision :: x,y,u,v,pi,akp
!
!  Consistency check:
      if(.not.present(wvl_p) .or. .not.present(urms))then
         print*,'obe_2state_weak_withav: One or both of the optional'
         print*,'arguments required for a calculation with Doppler'
         print*,'averaging is not specified...'
         call obe_error
      endif
!
!  Calculate the wavenumber
      pi = 4.d0 * datan(1.d0)
      akp = 2.0d0 * pi / (wvl_p*1.0d-9)
!
      Omega_p = Omega_p_0*2.d0*pi*1.d6
!
      gammaBtot_f = Gamma_AB/2.d0 + other_dephas_AB + frqwidth/2.d0
!
      rp = dcmplx(gammaBtot_f,-Delta_p)*2.d0*pi*1.d6/akp
      eta = (0.d0,1.d0)*rp/urms
!
      if(dimag(eta).gt.0.d0)then
         x = dreal(eta)
         y = dimag(eta)
         call obe_Faddeeva(x,y,u,v)
         w = dcmplx(u,v)
         complexint = (0.d0,1.d0)*pi*w
      elseif(dimag(eta).lt.0.d0)then
         x = dreal(eta)
         y = -dimag(eta)
         call obe_Faddeeva(x,y,u,v)
         w = dcmplx(u,v)
         complexint = dconjg((0.d0,1.d0)*pi*w)
      else
         print*,'obe_2state_weak_withav: Zero imaginary part for eta'
         call obe_error
      endif
      d = -(0.d0,1.d0)/urms*rp/eta
      rhoBA = (d-1.d0)*complexint - (d+1.d0)*(0.d0,1.d0)*pi
      rhoBA = rhoBA*(0.d0,1.d0)/(2.d0*urms*sqrt(pi)*akp)
      rhoBA = ((0.d0,1.d0)*Omega_p/2.d0) * rhoBA
!
      rhoAB = dconjg(rhoBA)
      rhoAA = 1.d0
!
      end subroutine obe_2state_weak_withav
  
!!!!!!!!!!!!!!!!!!
 
      subroutine obe_2state_withav(rhoAA,rhoAB)
!
!  Calculation beyond the weak field approximation with Doppler averaging
!
      implicit none
      double complex, intent(out) :: rhoAB
      double precision, intent(out) :: rhoAA
      double complex :: rhoBA
      double complex :: Omega_p,eta1,eta2,d1,d2,w,rp,rm,c1,c2
      double complex :: complexint1,complexint2
      double precision :: gammaBtot_f,sq,sq_f,fact1,fact2,rhoBB
      double precision :: x,y,u,v,pi,akp
!
!  Consistency check:
      if(.not.present(wvl_p) .or. .not.present(urms))then
         print*,'obe_2state_withav: One or both of the optional'
         print*,'arguments required for a calculation with Doppler'
         print*,'averaging is not specified...'
         call obe_error
      endif
!
!  Calculate the wavenumber
      pi = 4.d0 * datan(1.d0)
      akp = 2.0d0 * pi / (wvl_p*1.0d-9)
!
      gammaBtot_f = Gamma_AB/2.d0 + other_dephas_AB + frqwidth/2.d0
      sq_f = sqrt(gammaBtot_f**2+cdabs(Omega_p_0)**2*gammaBtot_f/Gamma_AB)
      fact1 = (gammaBtot_f + sq_f)/(2.d0*sq_f)
      fact2 = (gammaBtot_f - sq_f)/(2.d0*sq_f)
!
      Omega_p = Omega_p_0*2.d0*pi*1.d6
      sq = sq_f*2.d0*pi*1.d6
!
      rp = dcmplx(sq_f,-Delta_p)*2.d0*pi*1.d6/akp
      rm = dcmplx(sq_f, Delta_p)*2.d0*pi*1.d6/akp
      eta1 = (0.d0,1.d0)*rp/urms
      eta2 = (0.d0,1.d0)*rm/urms
!
      if(dimag(eta1).gt.0.d0)then
         x = dreal(eta1)
         y = dimag(eta1)
         call obe_Faddeeva(x,y,u,v)
         w = dcmplx(u,v)
         complexint1 = (0.d0,1.d0)*pi*w
      elseif(dimag(eta1).lt.0.d0)then
         x = dreal(eta1)
         y = -dimag(eta1)
         call obe_Faddeeva(x,y,u,v)
         w = dcmplx(u,v)
         complexint1 = dconjg((0.d0,1.d0)*pi*w)
      else
         print*,'obe_2state_withav: Zero imaginary part for eta1'
         call obe_error
      endif
      d1 = -(0.d0,1.d0)/urms*rp/eta1
!
      if(dimag(eta2).gt.0.d0)then
         x = dreal(eta2)
         y = dimag(eta2)
         call obe_Faddeeva(x,y,u,v)
         w = dcmplx(u,v)
         complexint2 = (0.d0,1.d0)*pi*w
      elseif(dimag(eta2).lt.0.d0)then
         x = dreal(eta2)
         y = -dimag(eta2)
         call obe_Faddeeva(x,y,u,v)
         w = dcmplx(u,v)
         complexint2 = dconjg((0.d0,1.d0)*pi*w)
      else
         print*,'obe_2state_withav: Zero imaginary part for eta2'
         call obe_error
      endif
      d2 = -(0.d0,1.d0)/urms*rm/eta2
!
      c1 = (d1-1.d0)*complexint1 - (d1+1.d0)*(0.d0,1.d0)*pi
      c2 = (d2-1.d0)*complexint2 - (d2+1.d0)*(0.d0,1.d0)*pi
      rhoBB = -cdabs(Omega_p)**2*dimag(c1 + c2)*gammaBtot_f/    &
                  (4.d0*sq*Gamma_AB)/(2.d0*urms*sqrt(pi)*akp)
      rhoAA = 1.d0 - rhoBB
      rhoBA = fact1*c1 + fact2*c2
      rhoBA = rhoBA*(0.d0,1.d0)/(2.d0*urms*sqrt(pi)*akp)
      rhoBA = ((0.d0,1.d0)*Omega_p/2.d0) * rhoBA
!
      rhoAB = dconjg(rhoBA)
!
      end subroutine obe_2state_withav
  
!!!!!!!!!!!!!!!!!!
 
      end subroutine obe_2state
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine obe_weakprb_3stladder                                       &
                        (Omega_p_0,Omega_c_0,Delta_p,Delta_c,                &
                         Gamma_AB,Gamma_BC,frqwidth_p,frqwidth_c,            &
                         other_dephas_AB,other_dephas_AC,iDoppler,rhoAB,     &
                         wvl_p,wvl_c,idir_c,urms)
!
!  Three-state ladder system with two co-linear fields.
!  Weak probe limit with or without Doppler averaging.
!
      implicit none
      double complex, intent(in) :: Omega_p_0
      double precision, intent(in) :: Omega_c_0,Delta_p,Delta_c
      double precision, intent(in) :: Gamma_AB,Gamma_BC
      double precision, intent(in) :: other_dephas_AB,other_dephas_AC
      double precision, intent(in) :: frqwidth_p,frqwidth_c
      double precision, intent(in) :: wvl_p,wvl_c,urms
      integer, intent(in) :: iDoppler,idir_c
      double complex, intent(out) :: rhoAB
      optional :: wvl_p,wvl_c,idir_c,urms
!
!  Check whether the calculation is with or without Doppler averaging
!  and proceed accordingly.
!
      if(iDoppler.eq.1)then
         call obe_weakprb_3stladder_withav(rhoAB)
      elseif(iDoppler.eq.0)then
         call obe_weakprb_3stladder_noav(rhoAB)
      else
         print*,'obe_weakprb_3stladder: illegal value of iDoppler ',iDoppler
         call obe_error
      endif
!
      contains
  
!!!!!!!!!!!!!!!!!!

      subroutine obe_weakprb_3stladder_noav(rhoAB)
!
!  Calculation without Doppler averaging
!
      implicit none
      double complex, intent(out) :: rhoAB
      double complex :: rhoBA
      double precision :: pi,Omega_c
      double complex :: Omega_p,a,b
!
      pi = 4.d0 * datan(1.d0)
!
      Omega_p = Omega_p_0*2d0*pi*1.d6
      Omega_c = Omega_c_0*2d0*pi*1.d6
!
      a = Gamma_AB/2.d0 + other_dephas_AB + frqwidth_p/2.d0 -     &
                                            (0.d0,1.d0)*Delta_p
      b = Gamma_BC/2.d0 + other_dephas_AC + frqwidth_p/2.d0 +     &
                 frqwidth_c/2.d0 - (0.d0,1.d0)*(Delta_p+Delta_c)
!  Conversion to angular frequencies expressed in 1/s:
      a = 2.d0*pi*a*1.d6
      b = 2.d0*pi*b*1.d6
!
      rhoBA = 1.d0/(a + Omega_c**2/4.d0/b)
      rhoBA = ((0.d0,1.d0)*Omega_p/2.d0) * rhoBA
!
      rhoAB = dconjg(rhoBA)
!
      end subroutine obe_weakprb_3stladder_noav
  
!!!!!!!!!!!!!!!!!!

      subroutine obe_weakprb_3stladder_withav(rhoAB)
!
!  Calculation with Doppler averaging
!
      implicit none
      double complex, intent(out) :: rhoAB
      double complex :: rhoBA
      double complex :: Omega_p,eta1,eta2,d,w,discr,rp,rc
      double complex :: complexint1,complexint2
      double precision :: gammaBtot_f,gammaCtot_f,Omega_c
      double precision :: Delta_R,x,y,u,v,pi,akp,akc,akdiff
!
!  Consistency check:
      if(.not.present(wvl_p) .or. .not.present(wvl_c)                 &
         .or. .not.present(idir_c) .or. .not.present(urms))then
         print*,'obe_weakprb_3stladder_withav: One or several of the optional'
         print*,'arguments required for a calculation with Doppler'
         print*,'averaging is not specified...'
         call obe_error
      endif
      if(iabs(idir_c).ne.1)then
         print*,'obe_weakprb_3stladder_withav: illegal value of idir_c'
         print*,idir_c
         call obe_error
      endif
!
!  Calculate the wavenumbers
      pi = 4.d0 * datan(1.d0)
      akp = 2.0d0 * pi / (wvl_p*1.0d-9)
      akc = 2.0d0 * pi / (wvl_c*1.0d-9)
      akdiff = akp + dble(idir_c)*akc
!
      Omega_p = Omega_p_0*2.d0*pi*1.d6
      Omega_c = Omega_c_0*2.d0*pi*1.d6
      Delta_R = Delta_p + Delta_c
!
      gammaBtot_f = Gamma_AB/2.d0 + other_dephas_AB + frqwidth_p/2.d0
      gammaCtot_f = Gamma_BC/2.d0 + other_dephas_AC + frqwidth_p/2.d0 +  &
                                                      frqwidth_c/2.d0    
!
      rp = dcmplx(gammaBtot_f,-Delta_p)*2.d0*pi*1.d6/akp
      rc = dcmplx(gammaCtot_f,-Delta_R)*2.d0*pi*1.d6/akdiff
      discr = sqrt( (rp-rc)**2-                                    &
        (1.d0,0.d0)*Omega_c**2/(akp*akdiff) )
      eta1 = (0.d0,1.d0)/(2.d0*urms)*(rp+rc+discr)
      eta2 = (0.d0,1.d0)/(2.d0*urms)*(rp+rc-discr)
!
      d = -(0.d0,1.d0)/urms*(rp-rc)/(eta1-eta2)
      if(dimag(eta1).gt.0.d0)then
         x = dreal(eta1)
         y = dimag(eta1)
         call obe_Faddeeva(x,y,u,v)
         w = dcmplx(u,v)
         complexint1 = (0.d0,1.d0)*pi*w
      elseif(dimag(eta1).lt.0.d0)then
         x = dreal(eta1)
         y = -dimag(eta1)
         call obe_Faddeeva(x,y,u,v)
         w = dcmplx(u,v)
         complexint1 = dconjg((0.d0,1.d0)*pi*w)
      else
         print*,'obe_weakprb_3stladder_withav: Zero imaginary part for eta1'
         call obe_error
      endif
      if(dimag(eta2).gt.0.d0)then
         x = dreal(eta2)
         y = dimag(eta2)
         call obe_Faddeeva(x,y,u,v)
         w = dcmplx(u,v)
         complexint2 = (0.d0,1.d0)*pi*w
      elseif(dimag(eta2).lt.0.d0)then
         x = dreal(eta2)
         y = -dimag(eta2)
         call obe_Faddeeva(x,y,u,v)
         w = dcmplx(u,v)
         complexint2 = dconjg((0.d0,1.d0)*pi*w)
      else
         print*,'obe_weakprb_3stladder_withav: Zero imaginary part for eta2'
         call obe_error
      endif
      rhoBA = (d-1.d0)*complexint1 - (d+1.d0)*complexint2
      rhoBA = rhoBA*(0.d0,1.d0)/(2.d0*urms*sqrt(pi)*akp)
      rhoBA = ((0.d0,1.d0)*Omega_p/2.d0) * rhoBA
!
      rhoAB = dconjg(rhoBA)
!
      end subroutine obe_weakprb_3stladder_withav
  
!!!!!!!!!!!!!!!!!!
 
      end subroutine obe_weakprb_3stladder
 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

      subroutine obe_weakprb_4stladder  &
             (Omega_p_0,Omega_c1_0,Omega_c2_0,Delta_p,Delta_c1,Delta_c2,  &
              Gamma_AB,Gamma_BC,frqwidth_p,frqwidth_c1,frqwidth_c2,       &   
              other_dephas_AB,other_dephas_AC,other_dephas_AD,            &
              iDoppler,rhoAB,wvl_p,wvl_c1,wvl_c2,idir_c1,idir_c2,urms)
!
!  Four-state ladder system with three co-linear fields.
!  Weak probe limit with or without Doppler averaging.
!
      implicit none
      double complex, intent(in) :: Omega_p_0
      double precision, intent(in) :: Omega_c1_0,Omega_c2_0,               &
                                      Delta_p,Delta_c1,Delta_c2,           &
                                      Gamma_AB,Gamma_BC,other_dephas_AB,   &
                                      other_dephas_AC,other_dephas_AD,     &
                                      frqwidth_p,frqwidth_c1,frqwidth_c2,  &
                                      wvl_p,wvl_c1,wvl_c2,urms
      integer, intent(in) :: iDoppler,idir_c1,idir_c2
      double complex, intent(out) :: rhoAB
      double complex :: Omega_p,a,b,c
      double precision :: pi,Omega_c1,Omega_c2
      optional :: wvl_p,wvl_c1,wvl_c2,idir_c1,idir_c2,urms
!
      pi = 4.d0 * datan(1.d0)
!
      Omega_p = Omega_p_0*2d0*pi*1.d6
      Omega_c1 = Omega_c1_0*2d0*pi*1.d6
      Omega_c2 = Omega_c2_0*2d0*pi*1.d6
!
      a = Gamma_AB/2.d0 + other_dephas_AB + frqwidth_p/2.d0 -      &
                                         (0.d0,1.d0)*Delta_p
      b = Gamma_BC/2.d0 + other_dephas_AC + frqwidth_p/2.d0 +      &
                 frqwidth_c1/2.d0 - (0.d0,1.d0)*(Delta_p+Delta_c1)
      c = other_dephas_AD + frqwidth_p/2.d0 +                      &
                      frqwidth_c1/2.d0 + frqwidth_c2/2.d0 -        &
                      (0.d0,1.d0)*(Delta_p+Delta_c1+Delta_c2)
!  Conversion to angular frequencies expressed in 1/s:
      a = 2.d0*pi*a*1.d6
      b = 2.d0*pi*b*1.d6
      c = 2.d0*pi*c*1.d6
!
!  Check whether the calculation is with or without Doppler averaging
!  and proceed accordingly
!
      if(iDoppler.eq.1)then
         call obe_weakprb_4stladder_withav(rhoAB)
      else
         call obe_weakprb_4stladder_noav(rhoAB)
      endif
!
      contains

!!!!!!!!!!!!!!!!!!

      subroutine obe_weakprb_4stladder_noav(rhoAB)
!
!  Calculation without Doppler averaging
!
      implicit none
      double complex, intent(out) :: rhoAB
      double complex :: rhoBA
!
      if(c.eq.(0.d0,0.d0))then
         print*,'obe_weakprb_4stladder_withav: The variable c is zero'
         print*,'for the choice of detunings and dephasings, which'
         print*,'the computation impossible as programmed.'
         call obe_error
      endif
!
      rhoBA = 1.d0/(a + Omega_c1**2/4.d0/(b+Omega_c2**2/4.d0/c))
      rhoBA = ((0.d0,1.d0)*Omega_p/2.d0) * rhoBA
!
      rhoAB = dconjg(rhoBA)
!
      end subroutine obe_weakprb_4stladder_noav

!!!!!!!!!!!!!!!!!!

      subroutine obe_weakprb_4stladder_withav(rhoAB)
!
!  Calculation with Doppler averaging
!
      implicit none
      double complex, intent(out) :: rhoAB
      double complex :: rhoBA
      double complex :: p,q,clam,cmu,ctheta,uu,vv,z1,z2,z3
      double complex :: w,complexint1,complexint2,complexint3
      double precision :: akp,akc1,akc2,dk1,dk2,u,v,x,y
!
!  Consistency check:
      if(.not.present(wvl_p) .or. .not.present(wvl_c1) .or.           & 
            .not.present(wvl_c2) .or.                                 &
            .not.present(idir_c1) .or. .not.present(idir_c2) .or.     &
            .not.present(urms))then
         print*,'obe_weakprb_4stladder_withav: One or several of the optional'
         print*,'arguments required for a calculation with Doppler'
         print*,'averaging is not specified...'
         call obe_error
      endif
      if(iabs(idir_c1).ne.1 .or. iabs(idir_c2).ne.1)then
         print*,'obe_weakprb_4stladder_withav: illegal value of '
         print*,'idir_c1 and/or idir_c2 ',idir_c1,idir_c2
         call obe_error
      endif
!
!  Calculate the wavenumbers
      akp = 2.0d0 * pi / (wvl_p*1.0d-9)
      akc1 = 2.0d0 * pi / (wvl_c1*1.0d-9)
      akc2 = 2.0d0 * pi / (wvl_c2*1.0d-9)
      dk1 = akp + dble(idir_c1)*akc1
      dk2 = akp + dble(idir_c1)*akc1 + dble(idir_c2)*akc2
!
      p = -(0.d0,1.d0)/urms*(b/dk1+c/dk2)
      q = -1.d0/urms**2*(4.d0*b*c+Omega_c2**2)/(4.d0*dk1*dk2)
      clam = -(0.d0,1.d0)/urms*(a/akp+b/dk1+c/dk2)
      cmu = -1.d0/urms**2*(a*b/(akp*dk1)+a*c/(akp*dk2)+b*c/(dk1*dk2) + &
                Omega_c2**2/(4.d0*dk1*dk2)+                      &
                Omega_c1**2/(4.d0*akp*dk1))
      ctheta = (0.d0,1.d0)/urms**3*(4.d0*a*b*c+                        &
                a*Omega_c2**2+c*Omega_c1**2)/              &
                (4.d0*akp*dk1*dk2)
      uu = 81.d0*ctheta**2+12.d0*ctheta*clam**3-54.d0*ctheta*clam*cmu- &
                3.d0*clam**2*cmu**2+12.d0*cmu**3
      uu = -27.d0*ctheta-2.d0*clam**3+9.d0*clam*cmu+3.d0*sqrt(uu)
      uu = (4.d0*uu)**(1.d0/3.d0)
      vv = 4.d0*(clam**2-3.d0*cmu)/uu
      z1 = (-2.d0*clam+uu+vv)/6.d0
      z2 = (-4.d0*clam-dcmplx(1.d0,sqrt(3.d0))*uu-                     &
                       dcmplx(1.d0,-sqrt(3.d0))*vv)/12.d0
      z3 = (-4.d0*clam-dcmplx(1.d0,-sqrt(3.d0))*uu-                    &
                       dcmplx(1.d0,sqrt(3.d0))*vv)/12.d0
!
      if(dimag(z1).gt.0.d0)then
         x = dreal(z1)
         y = dimag(z1)
         call obe_Faddeeva(x,y,u,v)
         w = dcmplx(u,v)
         complexint1 = (0.d0,1.d0)*pi*w
      elseif(dimag(z1).lt.0.d0)then
         x = dreal(z1)
         y = -dimag(z1)
         call obe_Faddeeva(x,y,u,v)
         w = dcmplx(u,v)
         complexint1 = dconjg((0.d0,1.d0)*pi*w)
      else
         print*,'obe_weakprb_4stladder_withav: Zero imaginary part for z1'
         print*,'This problem is likely to disappear if the wavelength'
         print*,'of one of the three fields is given a very slightly'
         print*,'different value.'
         call obe_error
      endif
      if(dimag(z2).gt.0.d0)then
         x = dreal(z2)
         y = dimag(z2)
         call obe_Faddeeva(x,y,u,v)
         w = dcmplx(u,v)
         complexint2 = (0.d0,1.d0)*pi*w
      elseif(dimag(z2).lt.0.d0)then
         x = dreal(z2)
         y = -dimag(z2)
         call obe_Faddeeva(x,y,u,v)
         w = dcmplx(u,v)
         complexint2 = dconjg((0.d0,1.d0)*pi*w)
      else
         print*,'obe_weakprb_4stladder_withav: Zero imaginary part for z2'
         call obe_error
      endif
      if(dimag(z3).gt.0.d0)then
         x = dreal(z3)
         y = dimag(z3)
         call obe_Faddeeva(x,y,u,v)
         w = dcmplx(u,v)
         complexint3 = (0.d0,1.d0)*pi*w
      elseif(dimag(z3).lt.0.d0)then
         x = dreal(z3)
         y = -dimag(z3)
         call obe_Faddeeva(x,y,u,v)
         w = dcmplx(u,v)
         complexint3 = dconjg((0.d0,1.d0)*pi*w)
      else
         print*,'obe_weakprb_4stladder_withav: Zero imaginary part for z3'
         call obe_error
      endif
!
      rhoBA =         (q+p*z1+z1**2)/((z1-z2)*(z1-z3))*complexint1
      rhoBA = rhoBA + (q+p*z2+z2**2)/((z2-z1)*(z2-z3))*complexint2
      rhoBA = rhoBA + (q+p*z3+z3**2)/((z3-z1)*(z3-z2))*complexint3
!
      rhoBA = -2.d0*rhoBA*(0.d0,1.d0)/(2.d0*urms*sqrt(pi)*akp)
!
      rhoBA = ((0.d0,1.d0)*Omega_p/2.d0) * rhoBA
!
      rhoAB = dconjg(rhoBA)
!
      end subroutine obe_weakprb_4stladder_withav
  
!!!!!!!!!!!!!!!!!!
 
      end subroutine obe_weakprb_4stladder
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine obe_Faddeeva(x,y,u,v)
!
!  Calculates the Faddeeva function of argument x + iy. The calculation
!  is entirely done by the internal subroutine obe_wofz.
!
!  At return, u and v are, respectively, the real part and imaginary
!  part of the Faddeeva function.
!
      implicit none
      double precision, intent(in) :: x,y
      double precision, intent(out) :: u,v
      logical :: flag
!
      call obe_wofz(x,y,u,v,flag)
!
      if(flag)then
         print*,'obe_Faddeeva: something went wrong in obe_wofz.'
         call obe_error
      endif
!
      contains
  
!!!!!!!!!!!!!!!!!!

      SUBROUTINE obe_WOFZ (XI, YI, U, V, FLAG)
!
! This is a Fortran-90 implementation of the Fortran-77 program
! published as Algorithm 680 of the Collected Algorithms from ACM (see
! reference below), without any substantive change to the original code.
! In particular, this version of WOFZ does not include the changes 
! proposed by M R Zaghloul which aim at improving the accuracy of the
! computation near the real axis. 
!
! While not the original software, this code contains material copyrighted
! by ACM and distributed within the licence provisions spelled out at
! https://www.acm.org/publications/policies/software-copyright-notice 
! In particular, please note that its usage is limited solely for
! academic, research and other similar noncommercial uses, and also that
! this code is being supplied "as is", without any support services from
! ACM. Neither ACM nor the author makes any representations or
! warranties, express or implied, including, without limitation, any
! representations or warranties of the merchantability or fitness for
! any particular purpose, or that the application of the software, will
! not infringe on any patents or other proprietary rights of others.
!
!  ALGORITHM 680, COLLECTED ALGORITHMS FROM ACM.
!  THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
!  VOL. 16, NO. 1, PP. 47.
!
!  GIVEN A COMPLEX NUMBER Z = (XI,YI), THIS SUBROUTINE COMPUTES
!  THE VALUE OF THE FADDEEVA-FUNCTION W(Z) = EXP(-Z**2)*ERFC(-I*Z),
!  WHERE ERFC IS THE COMPLEX COMPLEMENTARY ERROR-FUNCTION AND I
!  MEANS SQRT(-1).
!  THE ACCURACY OF THE ALGORITHM FOR Z IN THE 1ST AND 2ND QUADRANT
!  IS 14 SIGNIFICANT DIGITS; IN THE 3RD AND 4TH IT IS 13 SIGNIFICANT
!  DIGITS OUTSIDE A CIRCULAR REGION WITH RADIUS 0.126 AROUND A ZERO
!  OF THE FUNCTION.
!  ALL REAL VARIABLES IN THE PROGRAM ARE DOUBLE PRECISION.
!
!
!  THE CODE CONTAINS A FEW COMPILER-DEPENDENT PARAMETERS :
!     RMAXREAL = THE MAXIMUM VALUE OF RMAXREAL EQUALS THE ROOT OF
!                RMAX = THE LARGEST NUMBER WHICH CAN STILL BE
!                IMPLEMENTED ON THE COMPUTER IN DOUBLE PRECISION
!                FLOATING-POINT ARITHMETIC
!     RMAXEXP  = LN(RMAX) - LN(2)
!     RMAXGONI = THE LARGEST POSSIBLE ARGUMENT OF A DOUBLE PRECISION
!                GONIOMETRIC FUNCTION (DCOS, DSIN, ...)
!  THE REASON WHY THESE PARAMETERS ARE NEEDED AS THEY ARE DEFINED WILL
!  BE EXPLAINED IN THE CODE BY MEANS OF COMMENTS
!
!
!  PARAMETER LIST
!     XI     = REAL      PART OF Z
!     YI     = IMAGINARY PART OF Z
!     U      = REAL      PART OF W(Z)
!     V      = IMAGINARY PART OF W(Z)
!     FLAG   = AN ERROR FLAG INDICATING WHETHER OVERFLOW WILL
!              OCCUR OR NOT; TYPE LOGICAL;
!              THE VALUES OF THIS VARIABLE HAVE THE FOLLOWING
!              MEANING :
!              FLAG=.FALSE. : NO ERROR CONDITION
!              FLAG=.TRUE.  : OVERFLOW WILL OCCUR, THE ROUTINE
!                             BECOMES INACTIVE
!  XI, YI      ARE THE INPUT-PARAMETERS
!  U, V, FLAG  ARE THE OUTPUT-PARAMETERS
!
!  FURTHERMORE THE PARAMETER FACTOR EQUALS 2/SQRT(PI)
!
!  THE ROUTINE IS NOT UNDERFLOW-PROTECTED BUT ANY VARIABLE CAN BE
!  PUT TO 0 UPON UNDERFLOW;
!
!  REFERENCE - GPM POPPE, CMJ WIJERS; MORE EFFICIENT COMPUTATION OF
!  THE COMPLEX ERROR-FUNCTION, ACM TRANS. MATH. SOFTWARE.
!
!
      IMPLICIT NONE
!
      DOUBLE PRECISION, INTENT(IN) :: XI,YI
      DOUBLE PRECISION, INTENT(OUT) :: U,V
      LOGICAL, INTENT(OUT) :: FLAG
      DOUBLE PRECISION, PARAMETER :: FACTOR = 1.12837916709551257388D0,   &
                                     RMAXREAL = 0.5D+154,                 &
                                     RMAXEXP  = 708.503061461606D0,       &
                                     RMAXGONI = 3.53711887601422D+15
      DOUBLE PRECISION :: C,DAUX,H,H2,QLAMBDA,QRHO,       &
         RX,RY,SX,SY,TX,TY,U1,U2,V1,V2,W1,                &
         X,XABS,XABSQ,XAUX,XQUAD,XSUM,Y,YABS,YQUAD,YSUM
      INTEGER :: I,J,KAPN,N,NP1,NU
      LOGICAL :: A, B
!
      FLAG = .FALSE.
!
      XABS = DABS(XI)
      YABS = DABS(YI)
      X    = XABS/6.3
      Y    = YABS/4.4
!
!
!     THE FOLLOWING IF-STATEMENT PROTECTS
!     QRHO = (X**2 + Y**2) AGAINST OVERFLOW
!
      IF ((XABS.GT.RMAXREAL).OR.(YABS.GT.RMAXREAL)) GOTO 100
!
      QRHO = X**2 + Y**2
!
      XABSQ = XABS**2
      XQUAD = XABSQ - YABS**2
      YQUAD = 2*XABS*YABS
!
      A     = QRHO.LT.0.085264D0
!
      IF (A) THEN
!
!  IF (QRHO.LT.0.085264D0) THEN THE FADDEEVA-FUNCTION IS EVALUATED
!  USING A POWER-SERIES (ABRAMOWITZ/STEGUN, EQUATION (7.1.5), P.297)
!  N IS THE MINIMUM NUMBER OF TERMS NEEDED TO OBTAIN THE REQUIRED
!  ACCURACY
!
        QRHO  = (1-0.85D0*Y)*DSQRT(QRHO)
        N     = IDNINT(6 + 72*QRHO)
        J     = 2*N+1
        XSUM  = 1.0D0/J
        YSUM  = 0.0D0
        DO I=N, 1, -1
          J    = J - 2
          XAUX = (XSUM*XQUAD - YSUM*YQUAD)/I
          YSUM = (XSUM*YQUAD + YSUM*XQUAD)/I
          XSUM = XAUX + 1.0D0/J
        ENDDO
        U1   = -FACTOR*(XSUM*YABS + YSUM*XABS) + 1.0
        V1   =  FACTOR*(XSUM*XABS - YSUM*YABS)
        DAUX =  DEXP(-XQUAD)
        U2   =  DAUX*DCOS(YQUAD)
        V2   = -DAUX*DSIN(YQUAD)
!
        U    = U1*U2 - V1*V2
        V    = U1*V2 + V1*U2
!
      ELSE
!
!  IF (QRHO.GT.1.O) THEN W(Z) IS EVALUATED USING THE LAPLACE
!  CONTINUED FRACTION
!  NU IS THE MINIMUM NUMBER OF TERMS NEEDED TO OBTAIN THE REQUIRED
!  ACCURACY
!
!  IF ((QRHO.GT.0.085264D0).AND.(QRHO.LT.1.0)) THEN W(Z) IS EVALUATED
!  BY A TRUNCATED TAYLOR EXPANSION, WHERE THE LAPLACE CONTINUED FRACTION
!  IS USED TO CALCULATE THE DERIVATIVES OF W(Z)
!  KAPN IS THE MINIMUM NUMBER OF TERMS IN THE TAYLOR EXPANSION NEEDED
!  TO OBTAIN THE REQUIRED ACCURACY
!  NU IS THE MINIMUM NUMBER OF TERMS OF THE CONTINUED FRACTION NEEDED
!  TO CALCULATE THE DERIVATIVES WITH THE REQUIRED ACCURACY
!
        IF (QRHO.GT.1.0) THEN
          H    = 0.0D0
          KAPN = 0
          QRHO = DSQRT(QRHO)
          NU   = IDINT(3 + (1442/(26*QRHO+77)))
        ELSE
          QRHO = (1-Y)*DSQRT(1-QRHO)
          H    = 1.88*QRHO
          H2   = 2*H
          KAPN = IDNINT(7  + 34*QRHO)
          NU   = IDNINT(16 + 26*QRHO)
        ENDIF
!
        B = (H.GT.0.0)
!
        IF (B) QLAMBDA = H2**KAPN
!
        RX = 0.0
        RY = 0.0
        SX = 0.0
        SY = 0.0
!
        DO N=NU, 0, -1
          NP1 = N + 1
          TX  = YABS + H + NP1*RX
          TY  = XABS - NP1*RY
          C   = 0.5/(TX**2 + TY**2)
          RX  = C*TX
          RY  = C*TY
          IF ((B).AND.(N.LE.KAPN)) THEN
            TX = QLAMBDA + SX
            SX = RX*TX - RY*SY
            SY = RY*TX + RX*SY
            QLAMBDA = QLAMBDA/H2
          ENDIF
        ENDDO
!
        IF (H.EQ.0.0) THEN
          U = FACTOR*RX
          V = FACTOR*RY
        ELSE
          U = FACTOR*SX
          V = FACTOR*SY
        END IF
!
        IF (YABS.EQ.0.0) U = DEXP(-XABS**2)
!
      END IF
!
!
!
!  EVALUATION OF W(Z) IN THE OTHER QUADRANTS
!
!
      IF (YI.LT.0.0) THEN
!
        IF (A) THEN
          U2    = 2*U2
          V2    = 2*V2
        ELSE
          XQUAD =  -XQUAD
!
!
!         THE FOLLOWING IF-STATEMENT PROTECTS 2*EXP(-Z**2)
!         AGAINST OVERFLOW
!
          IF ((YQUAD.GT.RMAXGONI).OR.(XQUAD.GT.RMAXEXP)) GOTO 100
!
          W1 =  2*DEXP(XQUAD)
          U2  =  W1*DCOS(YQUAD)
          V2  = -W1*DSIN(YQUAD)
        END IF
!
        U = U2 - U
        V = V2 - V
        IF (XI.GT.0.0) V = -V
      ELSE
        IF (XI.LT.0.0) V = -V
      END IF
!
      RETURN
!
  100 FLAG = .TRUE.
      RETURN
!
      END SUBROUTINE obe_WOFZ

!!!!!!!!!!!!!!!!!!
 
      end subroutine obe_Faddeeva
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine obe_quad_rule(rule,npoints,x,w)
!
!  Returns the abscissas and weights for a Gauss-Legendre quadrature
!  (rule.eq.'Gaussian') or a trapezoidal quadrature (rule.eq.'trapezdl')
!  or an extended Simpson's quadrature (rule.eq.'Simpsons') or a Clenshaw-Curtis
!  quadrature (rule.eq.'Clenshaw'), of order npoints.
!  This subroutine allocates and returns the arrays x and w, of size
!  1:npoints, which contain, respectively, the abscissas and the weights
!  for an integration over the range [-1,1].
!
!  For a Gaussian quadrature, npoints must be 9, 12 or 16. npoints must be odd
!  for an extended Simpson's quadrature. The value of npoints is arbitrary for
!  a trapezoidal or a Clenshaw-Curtis quadrature, as long as it is positive.
!
!  The abscissas and weights for the Gaussian quadratures are copied from
!  Abramowitz and Stegun. Those for the Clenshaw-Curtis quadrature are
!  calculated by the subroutine clenshaw_curtis_compute , which is contained
!  in quad_rule.
!
      double precision, dimension(:), allocatable, intent(out) :: w,x
      integer, intent(in) :: npoints
      character*8, intent(in) :: rule
      integer :: istat
      double precision, dimension(9), parameter :: xGauss9 =   &
         (/ 0.d0,                                              &
            0.324253423403809d0, -0.324253423403809d0,         &
            0.613371432700590d0, -0.613371432700590d0,         &
            0.836031107326636d0, -0.836031107326636d0,         &
            0.968160239507626d0, -0.968160239507626d0 /)
      double precision, dimension(9), parameter :: wGauss9 =   &
         (/ 0.330239355001260d0,                               &
            0.312347077040003d0,  0.312347077040003d0,         &
            0.260610696402935d0,  0.260610696402935d0,         &
            0.180648160694857d0,  0.180648160694857d0,         &
            0.081274388361574d0,  0.081274388361574d0 /)
      double precision, dimension(12), parameter :: xGauss12 = &
         (/ 0.125233408511469d0, -0.125233408511469d0,         &
            0.367831498998180d0, -0.367831498998180d0,         &
            0.587317954286617d0, -0.587317954286617d0,         &
            0.769902674194305d0, -0.769902674194305d0,         &
            0.904117256370475d0, -0.904117256370475d0,         &
            0.981560634246719d0, -0.981560634246719d0 /)
      double precision, dimension(12), parameter :: wGauss12 = &
         (/ 0.249147045813403d0,  0.249147045813403d0,         &
            0.233492536538355d0,  0.233492536538355d0,         &
            0.203167426723066d0,  0.203167426723066d0,         &
            0.160078328543346d0,  0.160078328543346d0,         &
            0.106939325995318d0,  0.106939325995318d0,         &
            0.047175336386512d0,  0.047175336386512d0 /)
      double precision, dimension(16), parameter :: xGauss16 =       &
         (/ 0.095012509837637440185d0, -0.095012509837637440185d0,   &
            0.281603550779258913230d0, -0.281603550779258913230d0,   &
            0.458016777657227386342d0, -0.458016777657227386342d0,   &
            0.617876244402643748447d0, -0.617876244402643748447d0,   &
            0.755404408355003033895d0, -0.755404408355003033895d0,   &
            0.865631202387831743880d0, -0.865631202387831743880d0,   &
            0.944575023073232576078d0, -0.944575023073232576078d0,   &
            0.989400934991649932596d0, -0.989400934991649932596d0 /)
      double precision, dimension(16), parameter :: wGauss16 =       &
         (/ 0.189450610455068496285d0, 0.189450610455068496285d0,    &
            0.182603415044923588867d0, 0.182603415044923588867d0,    &
            0.169156519395002538189d0, 0.169156519395002538189d0,    &
            0.149595988816576732081d0, 0.149595988816576732081d0,    &
            0.124628971255533872052d0, 0.124628971255533872052d0,    &
            0.095158511682492784810d0, 0.095158511682492784810d0,    &
            0.062253523938647892863d0, 0.062253523938647892863d0,    &
            0.027152459411754094852d0, 0.027152459411754094852d0 /)
!
      if(npoints.gt.0)then
         if(allocated(x))deallocate(x)
         if(allocated(w))deallocate(w)
         allocate(x(npoints),w(npoints),stat=istat)
         if(istat.gt.0)then
            print*,'Error at the allocation of the arrays in obe_quad_rule'
            call obe_error
         endif
      else
         print*,'obe_quad_rule was called with npoints .leq. 0 ',npoints
         call obe_error
      endif
!
      if(rule.eq.'Gaussian')then
         if(npoints.eq.9)then
            x = xGauss9
            w = wGauss9
         elseif(npoints.eq.12)then
            x = xGauss12
            w = wGauss12
         elseif(npoints.eq.16)then
            x = xGauss16
            w = wGauss16
         else
            print*,'obe_quad_rule: The value of npoints is unsuitable for a'
            print*,'Gaussian quadrature.'
            call obe_error
         endif
      elseif(rule.eq.'trapezdl')then
         call obe_trapezoidal_compute
      elseif(rule.eq.'Simpsons')then
         call obe_Simpson_compute
      elseif(rule.eq.'Clenshaw')then
         call obe_clenshaw_curtis_compute(npoints,x,w)
      else
         print*,'obe_quad_rule: Unknown rule: ',rule
         call obe_error
      endif
!
      contains

!!!!!!!!!!!!!!!!!!

      subroutine obe_trapezoidal_compute
!
!  Abscissas and weights for a npoints-points trapezoidal integration
!  between -1 and 1.
!
      implicit none
      double precision :: step
      integer :: i
!
      step = 2.d0/dble(npoints-1)
!
      do i = 1,npoints
         x(i) = -1.d0 + dble(i-1)*step
         if(i.eq.1 .or. i.eq.npoints)then
            w(i) = 0.5d0
         else
            w(i) = 1.d0
         endif
      enddo
      w = w * step
!
      end subroutine obe_trapezoidal_compute
 
!!!!!!!!!!!!!!!!!!
 
      subroutine obe_Simpson_compute
!
!  Abscissas and weights for a npoints-points Simpson's integration
!  between -1 and 1.
!
      implicit none
      double precision :: s_simps,step
      integer :: i
!
      if(mod(npoints,2).eq.0)then
         print*,'obe_Simpson_compute: the number of points must be'
         print*,'odd for a Simpson''s integration.'
         call obe_error
      endif
!
      step = 2.d0/dble(npoints-1)
!
!  s_simps is used to generate the quadrature weights 4/3, 2/3, 4/3, etc.
! 
      x(1) = -1.d0
      w(1) = 1.d0/3.d0
      s_simps = -1.d0
      do i = 2,npoints-1
         x(i) = -1.d0 + dble(i-1)*step
         s_simps = -s_simps
         w(i) = 1.d0 + s_simps/3.d0
      enddo
      x(npoints) = 1.d0
      w(npoints) = 1.d0/3.d0
      w = w * step
!
      end subroutine obe_Simpson_compute
 
!!!!!!!!!!!!!!!!!!

      subroutine obe_clenshaw_curtis_compute ( order, x, w )
!
!  The original code, with the following changes:
!
!  The integer(kind = 4) declarations have been changed to
!      integer(kind = kind(1)) declarations.
!  The real(kind = 8) declarations have been changed to
!      real(kind = kind(1.d0)) declarations.
!  The D+00 have been changed into d0.
!  
!*****************************************************************************80
!
!! CLENSHAW_CURTIS_COMPUTE computes a Clenshaw Curtis quadrature rule.
!
!  Discussion:
!
!    Our convention is that the abscissas are numbered from left to right.
!
!    The rule is defined on [-1,1].
!
!    The integral to approximate:
!
!      Integral ( -1 <= X <= 1 ) F(X) dX
!
!    The quadrature rule:
!
!      Sum ( 1 <= I <= ORDER ) W(I) * F ( X(I) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 February 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER, the order of the rule.
!    1 <= ORDER.
!
!    Output, real ( kind = 8 ) X(ORDER), the abscissas.
!
!    Output, real ( kind = 8 ) W(ORDER), the weights.
!
  implicit none

  integer ( kind = kind(1) ) order

  real ( kind = kind(1.d0) ) b
  integer ( kind = kind(1) ) i
  integer ( kind = kind(1) ) j
  real ( kind = kind(1.d0) ), parameter :: r8_pi = 3.141592653589793d0
  real ( kind = kind(1.d0) ) theta
  real ( kind = kind(1.d0) ) w(order)
  real ( kind = kind(1.d0) ) x(order)

  if ( order < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CLENSHAW_CURTIS_COMPUTE - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of ORDER = ', order
!!  stop
    call obe_error
  end if

  if ( order == 1 ) then
    x(1) = 0.0d0
    w(1) = 2.0d0
    return
  end if

  do i = 1, order
    x(i) = cos ( real ( order - i, kind = kind(1.d0) ) * r8_pi &
               / real ( order - 1, kind = kind(1.d0) ) )
  end do

  x(1) = -1.0d0
  if ( mod ( order, 2 ) == 1 ) then
    x((order+1)/2) = 0.0d0
  end if
  x(order) = +1.0d0

  do i = 1, order

    theta = real ( i - 1, kind = kind(1.d0) ) * r8_pi &
          / real ( order - 1, kind = kind(1.d0) )

    w(i) = 1.0d0

    do j = 1, ( order - 1 ) / 2

      if ( 2 * j == ( order - 1 ) ) then
        b = 1.0d0
      else
        b = 2.0d0
      end if

      w(i) = w(i) - b * cos ( 2.0d0 * real ( j, kind = kind(1.d0) ) * theta )&
           / real ( 4 * j * j - 1, kind = kind(1.d0) )

    end do

  end do

  w(1)         =           w(1)         / real ( order - 1, kind = kind(1.d0) )
  w(2:order-1) = 2.0d0   * w(2:order-1) / real ( order - 1, kind = kind(1.d0) )
  w(order)     =           w(order)     / real ( order - 1, kind = kind(1.d0) )

  return

      end subroutine obe_clenshaw_curtis_compute

!!!!!!!!!!!!!!!!!!

      end subroutine obe_quad_rule


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine obe_get_Dopplerpar(fl_notinit_out,urms_out,n_v_values_out, &
                                    vmesh_out,fMvweight_out)
!
!  Returns the contents of the global variable fl_Doppler_notinit
!  through the argument fl_notinit_out. If fl_Doppler_notinit is false,
!  which is the case if and only if the parameters of the integration
!  over the velocity classes have been set through a prior call to
!  obe_set_Doppler, the allocatable arrays vmesh_out and fMvweight_out
!  are allocated and the contents of the global variables urms_saved,
!  n_v_values, vmesh and fMvweight are copied into the corresponding
!  variables of the list of arguments
!
      implicit none
      real(kd), allocatable, intent(out) :: vmesh_out(:),fMvweight_out(:)
      real(kd), intent(out) :: urms_out
      integer, intent(out) :: n_v_values_out
      logical, intent(out) :: fl_notinit_out
      integer :: ndim,istat
!
      fl_notinit_out = fl_Doppler_notinit
      if(.not.fl_Doppler_notinit)then
         ndim = size(vmesh)
         if(ndim.le.0 .or. size(fMvweight).ne.ndim)then
            print*,'obe_get_Dopplerpar: Logic error, ndim is not positive'
            print*,'or vmesh and fMvweight have inconsistent lengths.'
            print*,ndim,size(fMvweight)
            call obe_error
         endif
         if(allocated(vmesh_out))deallocate(vmesh_out)
         if(allocated(fMvweight_out))deallocate(fMvweight_out)
         allocate(vmesh_out(ndim),fMvweight_out(ndim),stat=istat)
         if(istat.gt.0)then
            print*,'obe_getDopplerpar: error at the allocation of'
            print*,'the arrays.'
            call obe_error
         endif
         urms_out = urms_saved
         n_v_values_out = n_v_values
         vmesh_out = vmesh
         fMvweight_out = fMvweight
      endif

      end subroutine obe_get_Dopplerpar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine obe_find_Rabif(ihigher,dipmom,campl,rabif)
!
!  Given a complex dipole moment expressed in C m and a
!  complex electric field amplitude expressed in V/m, 
!  returns the corresponding Rabi frequency expressed
!  as a frequency in MHz. Specifically, dipmom is
!  expected to be the matrix element <i | epsilon dot D | j>
!  if state i is higher in energy than state j and
!  <i | epsilon^* dot D | j> if state i is lower in energy
!  than state j, where D is the dipole operator and epsilon is
!  the complex polarisation vector of the field. At return,
!  rabif is Omega_{ij}. Whether state i is higher or lower
!  in energy than state j is identified by the value of
!  ihigher: ihigher = 1 means that i is higher in energy,
!  and ihigher = 0 means that i is lower in energy.
!
      double complex, intent(in) :: dipmom,campl
      double complex, intent(out) :: rabif
      integer, intent(in) :: ihigher
      double precision, parameter :: pi = 4.d0*datan(1.d0)
!
      if(ihigher.eq.1)then
         rabif = campl*dipmom*1.d-6/(2.d0*pi*hbar)
      elseif(ihigher.eq.0)then
         rabif = dconjg(campl)*dipmom*1.d-6/(2.d0*pi*hbar)
      else
         print*,'obe_find_Rabif: the value of ihigher is illegal.'
         print*,ihigher
         call obe_error
      endif
!
      end subroutine obe_find_Rabif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine obe_find_campl(ihigher,dipmom,rabif,campl)
!
!  Given a complex dipole moment expressed in C m and a
!  complex Rabi frequency expressed as a frequency in MHz,
!  returns the corresponding complex electric field amplitude
!  expressed in V/m.
!  
!  The variables dipmom and rabif are as in obe_find_Rabif,
!  and all necessary explanations can be found in that subroutine.
!
!
      double complex, intent(in) :: dipmom,rabif
      double complex, intent(out) :: campl
      integer, intent(in) :: ihigher
      double precision, parameter :: pi = 4.d0*datan(1.d0)
!
      campl = rabif*1.d6*2.d0*pi*hbar/dipmom
      if(ihigher.eq.0)then
         campl = dconjg(campl)
      elseif(ihigher.ne.1)then
         print*,'obe_find_campl: the value of ihigher is illegal.'
         print*,ihigher
         call obe_error
      endif
!
      end subroutine obe_find_campl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine obe_transform_rhs_mat(nindx,rhs_mat,c2)
!
!  Transform rhs_mat into a double complex matrix acting on a column vector
!  of real populations and complex coherences, returned through the array c2.
!  nindx is meant to contain the indexes of the different populations
!  and coherences, as set up in the calling program.
!
      integer, parameter :: neqs = nst*nst
      real(kd), dimension(neqs,neqs), intent(in) :: rhs_mat
      integer, dimension(nst,nst), intent(in) :: nindx
      complex(kd), dimension(neqs,neqs), intent(out) :: c2
      complex(kd), dimension(1,neqs) :: cmplxcoh1,cmplxcoh2
      integer :: m,n,k,l,i,j,mre,mim
!
      c2 = (0.0_kd,0.0_kd)
      n = 0
      do l = 1,nst
         do k = 1,l
            n = n + 1
            do j = 1,nst
               do i = 1,j
                  if(i.eq.j)then
                     call ldbl_pop_index(i,m)
                     c2(nindx(k,l),nindx(i,j)) = c2(nindx(k,l),nindx(i,j)) +  &
                        (1.0_kd,0.0_kd)*rhs_mat(n,m)
                  else
                     call ldbl_coher_index(i,j,mre,mim)
                     c2(nindx(k,l),nindx(i,j)) = c2(nindx(k,l),nindx(i,j)) +  &
                        (0.5_kd,0.0_kd)*rhs_mat(n,mre)
                     c2(nindx(k,l),nindx(j,i)) = c2(nindx(k,l),nindx(j,i)) +  &
                        (0.5_kd,0.0_kd)*rhs_mat(n,mre)
                     c2(nindx(k,l),nindx(i,j)) = c2(nindx(k,l),nindx(i,j)) +  &
                        (0.0_kd,-0.5_kd)*rhs_mat(n,mim)
                     c2(nindx(k,l),nindx(j,i)) = c2(nindx(k,l),nindx(j,i)) +  &
                        (0.0_kd,0.5_kd)*rhs_mat(n,mim)
                  endif
               enddo
            enddo
            if(k.ne.l)then
               n = n + 1
               do j = 1,nst
                  do i = 1,j
                     if(i.eq.j)then
                        call ldbl_pop_index(i,m)
                        c2(nindx(l,k),nindx(i,j)) = c2(nindx(l,k),nindx(i,j)) +&
                           (1.0_kd,0.0_kd)*rhs_mat(n,m)
                     else
                        call ldbl_coher_index(i,j,mre,mim)
                        c2(nindx(l,k),nindx(i,j)) = c2(nindx(l,k),nindx(i,j)) +&
                           (0.5_kd,0.0_kd)*rhs_mat(n,mre)
                        c2(nindx(l,k),nindx(j,i)) = c2(nindx(l,k),nindx(j,i)) +&
                           (0.5_kd,0.0_kd)*rhs_mat(n,mre)
                        c2(nindx(l,k),nindx(i,j)) = c2(nindx(l,k),nindx(i,j)) +&
                           (0.0_kd,-0.5_kd)*rhs_mat(n,mim)
                        c2(nindx(l,k),nindx(j,i)) = c2(nindx(l,k),nindx(j,i)) +&
                           (0.0_kd,0.5_kd)*rhs_mat(n,mim)
                     endif
                  enddo
               enddo
            endif
         enddo
      enddo
!
      do l = 1,nst
         do k = 1,l
            if(k.eq.l)cycle
            cmplxcoh1(1,:) = c2(nindx(k,l),:) + (0.d0,1.d0)*c2(nindx(l,k),:)
            cmplxcoh2(1,:) = c2(nindx(k,l),:) - (0.d0,1.d0)*c2(nindx(l,k),:)
            c2(nindx(k,l),:) = cmplxcoh1(1,:)
            c2(nindx(l,k),:) = cmplxcoh2(1,:)
         enddo
      enddo
!
      end subroutine obe_transform_rhs_mat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine obe_Doppler_av_st(rhovec_st_av,urms)
!
!  Returns the steady state density matrix with Doppler averaging done
!  semi-analytically, using the Faddeeva function.
!
      implicit none
      integer, parameter :: neqs = nst*nst
      integer, parameter :: lwork = 2*(neqs-1)
      integer, parameter :: lrwork = 8*(neqs-1)
      real(kd), dimension (neqs), intent(out) :: rhovec_st_av
      double precision, optional, intent(in) :: urms
      real(kd), dimension (:,:), allocatable :: rhs_mat
      complex(kd), dimension (:,:), allocatable :: c2
      complex(kd), dimension (:,:), allocatable :: amat,bmat,c20,c21,vl,vr
      complex(kd), dimension (:), allocatable :: alpha,beta,cvec,cnotnul,cb
      integer, dimension (:), allocatable :: ipiv
      complex(kd), dimension (:), allocatable :: work
      real(kd), dimension (:), allocatable :: rwork
      real(kd) :: v
      integer, dimension (nst,nst) :: nindx
      integer :: i,j,k,l,n,m,mre,mim
      integer :: j_state,ivel
      integer :: info,lda,ldb,ldvl,ldvr,nrhs,istat
      character*1 :: jobvl,jobvr
      complex(kd) :: csum,cmu,csumn,csumd,coef
      double precision, save :: pi,sqrtpi
      logical, save :: fl_init = .true.
      external :: cgesv,cggev,zgesv,zggev
!
      if(nst.lt.2)then
         print*,'obe_Doppler_av_st: unsuitable values of nst ',nst
         call obe_error
      endif
!
      if(fl_init)then
         pi = 4.d0*datan(1.d0)
         sqrtpi = sqrt(pi)
         fl_init = .false.
      endif
      allocate(rhs_mat(neqs,neqs),c2(neqs,neqs),amat(neqs-1,neqs-1), &
         bmat(neqs-1,neqs-1),c20(neqs-1,neqs-1),c21(neqs-1,neqs-1),  &
         vl(neqs-1,neqs-1),vr(neqs-1,neqs-1),alpha(neqs-1),          &
         beta(neqs-1),cvec(neqs-1),cnotnul(neqs-1),cb(neqs-1),       &
         ipiv(neqs-1),work(lwork),rwork(lrwork),stat=istat)
         if(istat.gt.0)then
            print*,'obe_Doppler_av_st: error at the allocation of'
            print*,'the arrays.'
            call obe_error
         endif
!
!  Prepare the indexes labelling the elements of the density matrix.
!
      nindx = 0
!  Put the population of state 1 first.
      n = 1
      nindx(1,1) = n
      do k = 1,nst
         do l = 1,nst
            if(k.eq.1 .and. l.eq.1)cycle ! Already indexed.
            n = n + 1
            nindx(k,l) = n
         enddo
      enddo
      if(n.ne.neqs)then
         print*,'obe_Doppler_av_st: logic error.',n,neqs
         call obe_error
      endif
!
!  Big loop. At the first run, calculate the relevant arrays for a zero
!  velocity, at the second for a velocity of 1 m/s.
!
      do ivel = 0,1
!
         if(ivel.eq.0)then
            v = 0.0_kd
         else
            v = 1.0_kd  ! Set to 1 m/s
         endif
!
!  Get the matrix representing the right-hand sides in the form
!  of a complex matrix acting on populations and complex coherences.
         call obe_set_rhsmat(rhs_mat,v)
         call obe_transform_rhs_mat(nindx,rhs_mat,c2)
!
!  Set the population of state 1 to 1, hence the corresponding row and
!  the corresponding column of c2 go and the column vector cb is
!  the negative of the latter... 
!
         j_state = 1
         if(nindx(j_state,j_state).ne.1)then
            print*,'obe_Doppler_av_st: j_state must be 1'
            print*,'in the current state of development of the code.'
            call obe_error
         endif
         cb(1:neqs-1) = -c2(2:neqs,1)
!
!  Subtract the relevant quantity from the entries of c2 corresponding
!  to populations.
!
!  The following code applies only to the case where
!  nindx(j_state,j_state).eq.1. This is first (re)checked.
         if(nindx(j_state,j_state).ne.1)then
            print*,'obe_Doppler_av_st: j_state must be 1'
            print*,'in the current state of development of the code.'
            call obe_error
         endif
!
         do k = 1,nst
            if(k.eq.j_state)cycle
            c2(2:neqs,nindx(k,k)) = c2(2:neqs,nindx(k,k)) + cb(1:neqs-1)
         enddo
!
!  Store the result
         if(ivel.eq.0)c20(1:neqs-1,1:neqs-1) = c2(2:neqs,2:neqs)
         if(ivel.eq.1)c21(1:neqs-1,1:neqs-1) = c2(2:neqs,2:neqs)
!
!  End of the big loop
!
      enddo
!
      c21 = c21 - c20
!
!  Next step: Find the null space and the generalized eigenvalues.
!
      jobvl = 'V'
      jobvr = 'V'
      n = neqs - 1
      lda = neqs - 1
      ldb = neqs - 1
      ldvl = neqs - 1
      ldvr = neqs - 1
      amat = c20
      bmat = c21
      if(kd.eq.kind(1.0))then
         call cggev(jobvl,jobvr,n,amat,lda,bmat,ldb,alpha,beta,      &
                    vl,ldvl,vr,ldvr,work,lwork,rwork,info)
      elseif(kd.eq.kind(1.d0))then
         call zggev(jobvl,jobvr,n,amat,lda,bmat,ldb,alpha,beta,      &
                    vl,ldvl,vr,ldvr,work,lwork,rwork,info)
      else
         print*,'Logic error detected in obe_Doppler_av_st'
         call obe_error
      endif
      if(info.ne.0)then
         print*,'obe_Doppler_av_st: info is not zero at the'
         print*,'return from c/zggev.',info
         call obe_error
      endif
!
!  Calculate the coefficients of the expansion of the b-vector (the
!  right-hand side) on the generalized eigenvectors corresponding
!  to a finite eigenvalue (i.e., those eigenvectors orthogonal to the
!  null space of c20), and proceed with the calculation.
      cnotnul = (0.0_kd,0.0_kd)
      do n = 1,neqs-1
!  The determination of the null space effected by the following line of
!  code may need to be revisited if inappropriate for the case at hand.
         if(abs(beta(n)).le.1.e-10_kd)cycle
         cmu = alpha(n)/beta(n)
         cvec = matmul(c21,vr(:,n))
         csumn = (0.0_kd,0.0_kd)
         csumd = (0.0_kd,0.0_kd)
         do m = 1,neqs-1
            csumn = csumn + conjg(vl(m,n))*cb(m)
            csumd = csumd + conjg(vl(m,n))*cvec(m)
         enddo
         call obe_calc_coef_D ! This internal subroutine calculates 
                              ! the variable coef used below.
         cnotnul = cnotnul + coef*vr(:,n)
         cb = cb - (csumn/csumd)*cvec
      enddo
!
!  Contribution of the null space: we solve L_0 y = b - ...
      n = neqs - 1
      nrhs = 1
      amat = c20
      lda = neqs - 1
      ldb = neqs - 1
      if(kd.eq.kind(1.0))then
         call cgesv(n,nrhs,amat,lda,ipiv,cb,ldb,info)
      elseif(kd.eq.kind(1.d0))then
         call zgesv(n,nrhs,amat,lda,ipiv,cb,ldb,info)
      else
         print*,'Logic error 2 detected in obe_Doppler_av_st'
         call obe_error
      endif
      if(info.ne.0)then
         print*,'obe_Doppler_av_st: info is not zero at the'
         print*,'return from c/zgesv.',info
         call obe_error
      endif
!
!  Solution vector: the sum of the contribution of the null space,
!  which is in cb, with the rest, which is in cnotnul.
      cb = cb + cnotnul
!
!  Arrange the populations and the coherences in the required format 
!
      n = 1
      rhovec_st_av(1) = 1.0_kd
      do j = 1,nst
         do i = 1,j
            if(i.eq.1 .and. j.eq.1)then
               cycle
            elseif(i.ne.j)then
               n = n + 1
               rhovec_st_av(n) = real(cb(nindx(i,j)-1),kd)
               n = n + 1
               rhovec_st_av(n) = dimag(1.d0*cb(nindx(i,j)-1))
            else
               n = n + 1
               rhovec_st_av(n) = real(cb(nindx(i,j)-1),kd)
               rhovec_st_av(1) = rhovec_st_av(1) - rhovec_st_av(n)
            endif
         enddo
      enddo
!
      deallocate(rhs_mat,c2,amat,bmat,c20,c21,vl,vr,alpha, &
         beta,cvec,cnotnul,cb,ipiv,work,rwork)
!
      contains

!!!!!!!!!!!!!!!!!!

      subroutine obe_calc_coef_D
!
      implicit none
      double complex :: eta,w,complexint
      double precision :: x,y,u,v
!
      eta = -cmu/urms
!
      if(dimag(eta).gt.0.d0)then
         x = dreal(eta)
         y = dimag(eta)
         call obe_Faddeeva(x,y,u,v)
         w = dcmplx(u,v)
         complexint = (0.d0,1.d0)*pi*w
      elseif(dimag(eta).lt.0.d0)then
         x = dreal(eta)
         y = -dimag(eta)
         call obe_Faddeeva(x,y,u,v)
         w = dcmplx(u,v)
         complexint = dconjg((0.d0,1.d0)*pi*w)
      else
         print*,'obe_calc_coef_D: Zero imaginary part for eta'
         print*,cmu,eta
         call obe_error
      endif
! 
      coef = complexint/(urms*sqrtpi)
!
      coef = csumn * coef / csumd
!
      end subroutine obe_calc_coef_D

!!!!!!!!!!!!!!!!!!

      end subroutine obe_Doppler_av_st

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine obe_steadystate_ladder(iDoppler,popfinal,rho_vec,urms)
!
!  Calculation of the steady state density matrix for a ladder system in
!  the weak probe approximation. This is done within the rate equation
!  approximation - which is not an approximation here as for a ladder
!  system in the weak probe approximation the populations are
!  invariable: the lower energy states are populated, the higher
!  energy states have a zero population, and there is no optical
!  pumping. The populations are passed to this subroutine through
!  the array popfinal.
!
      implicit none
      integer, parameter :: neqs = nst*nst
      real(kd), dimension (neqs), intent(out) :: rho_vec
      double precision, dimension (nmn:nmx), intent(in) :: popfinal
      double precision, optional, intent(in) :: urms
      integer, intent(in) :: iDoppler
      real(kd), dimension (:,:), allocatable, save :: arhs
      integer, dimension (:), allocatable, save :: ipiv
      integer :: istat
      logical, save :: fl_init_rate = .true.
!
      if(.not.fl_wkprbapp)then
         print*,'obe_steadystate_ladder can be used only for calculations'
         print*,'in the weak probe approximation.'
         call obe_error
      endif
!
!  Call the appropriate internal subprograms to obtain the coherences from
!  the populations:
      if(iDoppler.eq.0)then
         call obe_ladder_noDoppler
      elseif(iDoppler.eq.1)then
         if(.not.present(urms))then
            print*,'obe_steadystate_ladder: The optional argument urms'
            print*,'should be specified if iDoppler = 1.'
            call obe_error
         endif
         call obe_ladder_withDoppler
      else
         print*,'obe_steadystate_ladder: illegal value of iDoppler ',iDoppler
         call obe_error
      endif
!
      contains
 
!!!!!!!!!!!!!!!!!!!

      subroutine obe_ladder_noDoppler
!
!  Given the populations, calculate the coherences in the rate equation
!  approximation. No Doppler averaging is done.
!
      implicit none
      integer, parameter :: neqs = nst*nst
      real(kd), dimension (:,:), allocatable :: arhs
      integer, dimension (:), allocatable :: ipiv
      logical, dimension (:), allocatable :: flzeros
      integer :: i,info,lda,ldb,m,numb,nrhs,istat
      external dgesv,sgesv
!
!  Prepare the rate equation calculation and obtain the
!  matrix used to calculate the right-hand sides.
      if(allocated(arhs))deallocate(arhs)
      if(allocated(ipiv))deallocate(ipiv)
      if(allocated(flzeros))deallocate(flzeros)
      allocate(arhs(neqs,neqs),ipiv(neqs),flzeros(neqs),stat=istat)
      if(istat.gt.0)then
         print*,'obe_ladder_noDoppler: Error at the allocation of arhs etc.'
         call obe_error
      endif
      call obe_set_rhsmat(arhs)
!
!  First calculate the rhs of the system of linear equations. Start by
!  setting all the elements of rho_vec to zero, except the populations.
      rho_vec = 0.0_kd
      do i = nmn,nmx
         call obe_pop_index(i,m)
         rho_vec(m) = popfinal(i)
      enddo
!  Calculate the rhs of the system of linear equations.
      rho_vec = -matmul(arhs,rho_vec)
!
!  Before modifying arhs, take note of where it has rows of zeroes.
      do i = 1,neqs
         flzeros(i) = all(arhs(i,:) .eq. 0.0_kd)
      enddo
!
!  The populations should not be recalculated (which amounts to putting 1
!  on the corresponding diagonal elements of the matrix of the system, and
!  0 in the corresponding off-diagonal elements.
      do i = 1,nst
         call ldbl_pop_index(i,m)
         arhs(m,:) = 0.0_kd
         arhs(:,m) = 0.0_kd
         arhs(m,m) = 1.0_kd
      enddo
!
!  Check for rows of zeroes. For such rows, the corresponding element of the
!  density matrix remains constant under time propagation. Hence, do not
!  recalculate them.
      do i = 1,neqs
         if(flzeros(i))arhs(i,i) = 1.0_kd     
      enddo
!
!  Check for columns of zeroes. The corresponding elements of rho_vec should
!  not be recalculated - hence put 1 on the corresponding diagonal element
!  of the matrix of the system and 0 on the corresponding off-diagonal elements.
      do i = 1,neqs
         if(all(arhs(:,i) .eq. 0.0_kd))then
            arhs(i,:) = 0.0_kd
            arhs(i,i) = 1.0_kd
         endif
      enddo
!
!  Solve the system of linear equations with rho_vec as their rhs.
!  Variables passed to s/dgesv (solution of system of linear equations,
!  from lapack):
!    numb: the number of equations in the linear system to be solved
!    nrhs: number of right-hand sides (here one only)
!    lda: leading dimension of the array arhs
!    rho_vec: the column vector(s) forming the right-hand side(s) of the system
!       (here there is only one right-hand side; this array is overwritten
!       by s/dgesv)
!    ldb: the leading dimension of the array rho_vec
!    info: an integer variable whose value is set to zero by s/dgesv in case
!       of successful completion or to another value in case of problem
!       (see online information about s/dgesv for details)
!
      numb = neqs
      nrhs = 1
      lda = neqs
      ldb = neqs
!
      if(kd.eq.kind(1.0))then
         call sgesv(numb,nrhs,arhs,lda,ipiv,rho_vec,ldb,info)
      elseif(kd.eq.kind(1.d0))then
         call dgesv(numb,nrhs,arhs,lda,ipiv,rho_vec,ldb,info)
      else
         print*,'Logic error detected in obe_ladder_noDoppler'
         call obe_error
      endif
!
!  info.ne.0 : error condition
!  info.eq.0 : normal return from s/dgesv.
!
      if(info.ne.0)then
         print*,'obe_ladder_noDoppler: s/dgesv returns with info = ',info
         call obe_error
      endif
!
!  At this stage, rho_vec contains the correct coherences (within the rate
!  equation approximation) but not the correct populations. Replace the
!  (wrong) populations by those passed to obe_calc_coher:
      do i = nmn,nmx
         call obe_pop_index(i,m)
         rho_vec(m) = popfinal(i)
      enddo
!
!  That's all...
      deallocate(arhs,ipiv,flzeros)
!
      end subroutine obe_ladder_noDoppler

!!!!!!!!!!!!!!!!!!!

      subroutine obe_ladder_withDoppler
!
!  Given the populations, calculate the coherences in the rate equation
!  approximation. Doppler averaging is done semi-analytically.
!
      implicit none
      integer, parameter :: neqs = nst*nst
      integer, parameter :: lwork = 2*neqs
      integer, parameter :: lrwork = 8*neqs
      real(kd), dimension (:,:), allocatable :: rhs_mat
      complex(kd), dimension (:,:), allocatable :: c2,c20,c21
      complex(kd), dimension (:), allocatable :: cb
      complex(kd), dimension (:,:), allocatable :: amat,bmat,vl,vr
      complex(kd), dimension (:), allocatable :: alpha,beta,cvec,cnotnul
      complex(kd), dimension (:), allocatable :: work
      real(kd), dimension (:), allocatable :: rwork
      integer, dimension (:), allocatable :: ipiv
      logical, dimension (:), allocatable :: flzeros
      integer, dimension (nst,nst) :: nindx
      complex(kd) :: cmu,csumn,csumd,coef
      double complex :: eta,w,complexint
      real(kd) :: vel
      double precision :: x,y,u,v
      double precision, save :: pi,sqrtpi
      integer :: i,info,lda,ldb,m,numb,nrhs,istat,ivel,j,k,l,n,ldvl,ldvr
      logical, save :: fl_init = .true.
      character*1 :: jobvr,jobvl
      external :: cgesv,cggev,zgesv,zggev
!
      if(fl_init)then
         pi = 4.d0*datan(1.d0)
         sqrtpi = sqrt(pi)
         fl_init = .false.
      endif
!
!  Prepare the rate equation calculation and obtain the
!  matrix used to calculate the right-hand sides.
      if(allocated(rhs_mat))deallocate(rhs_mat)
      if(allocated(c2))deallocate(c2)
      if(allocated(c20))deallocate(c20)
      if(allocated(c21))deallocate(c21)
      if(allocated(cb))deallocate(cb)
      if(allocated(ipiv))deallocate(ipiv)
      if(allocated(flzeros))deallocate(flzeros)
      allocate(rhs_mat(neqs,neqs),c2(neqs,neqs),c20(neqs,neqs),c21(neqs,neqs),  &
               cb(neqs),ipiv(neqs),flzeros(neqs),amat(neqs,neqs),               &
               bmat(neqs,neqs),vl(neqs,neqs),vr(neqs,neqs),alpha(neqs),         &
               beta(neqs),cvec(neqs),cnotnul(neqs),work(lwork),rwork(lrwork),   &
               stat=istat)
      if(istat.gt.0)then
         print*,'obe_ladder_withDoppler: Error at the allocation of the arrays.'
         call obe_error
      endif
!
!  Prepare the indexes labelling the elements of the density matrix.
!
      nindx = 0
!  Put the population of state 1 first.
      n = 1
      nindx(1,1) = n
      do k = 1,nst
         do l = 1,nst
            if(k.eq.1 .and. l.eq.1)cycle ! Already indexed.
            n = n + 1
            nindx(k,l) = n
         enddo
      enddo
      if(n.ne.neqs)then
         print*,'obe_ladder_withDoppler: logic error.',n,neqs
         call obe_error
      endif
!
!  Big loop. At the first run, calculate the relevant arrays for a zero
!  velocity, at the second for a velocity of 1 m/s.
!
      do ivel = 0,1
!
         if(ivel.eq.0)then
            vel = 0.0_kd
         else
            vel = 1.0_kd  ! Set to 1 m/s
         endif
!
!  Get the matrix representing the right-hand sides in the form
!  of a complex matrix acting on populations and complex coherences.
         call obe_set_rhsmat(rhs_mat,vel)
         call obe_transform_rhs_mat(nindx,rhs_mat,c2)
!
!  First calculate the rhs of the system of linear equations. Start by
!  setting all the elements of the equivalent of rho_vec to zero, except
!  the populations. It is assumed that the results do not depend on the
!  atom's velocity; whether this is actually the case is not checked.
         if(ivel.eq.0)then
            cb = 0.0_kd
            do i = nmn,nmx
               cb(nindx(i+noffset,i+noffset)) = popfinal(i)
            enddo
!  Calculate the rhs of the system of linear equations.
            cb = -matmul(c2,cb)
         endif
!
!  Before modifying c2, take note of where it has rows of zeroes.
!  it is assumed that the non-zero elements of c2 are nornally of
!  order unity, so that any element less than 1e-12 in absolute
!  magnitude can safely be assumed to be zero.
         flzeros = .false.
         do k = 1,neqs
            flzeros(k) = all(abs(c2(k,:)) .lt. 1.e-12_kd)
         enddo
!
!  The populations should not be recalculated (which amounts to putting 1
!  on the corresponding diagonal elements of the matrix of the system, and
!  0 in the corresponding off-diagonal elements.
         do i = 1,nst
            c2(nindx(i,i),:) = (0.0_kd,0.0_kd)
            c2(:,nindx(i,i)) = (0.0_kd,0.0_kd)
            c2(nindx(i,i),nindx(i,i)) = (1.0_kd,0.0_kd)
         enddo
!
!  Check for rows of zeroes. For such rows, the corresponding element of the
!  density matrix remains constant under time propagation. Hence, do not
!  recalculate them.
         do k = 1,neqs
            if(flzeros(k))c2(k,k) = (1.0_kd,0.0_kd)     
         enddo
!
!  Check for columns of zeroes. The corresponding elements of rho_vec should
!  not be recalculated - hence put 1 on the corresponding diagonal element
!  of the matrix of the system and 0 on the corresponding off-diagonal elements.
!  See comment above about the 1.e-12.
         do k = 1,neqs
            if(all(abs(c2(:,k)) .lt. 1.e-12_kd))then
               c2(k,:) = (0.0_kd,0.0_kd)
               c2(k,k) = (1.0_kd,0.0_kd)
            endif
         enddo
!
!  Store the result
!
      if(ivel.eq.0)c20 = c2
      if(ivel.eq.1)c21 = c2
!
!  End of the big loop
!
      enddo
!
      c21 = c21 - c20
!
!  Next step: Find the null space and the generalized eigenvalues.
!
      jobvl = 'V'
      jobvr = 'V'
      n = neqs
      lda = neqs
      ldb = neqs
      ldvl = neqs
      ldvr = neqs
      amat = c20
      bmat = c21
      if(kd.eq.kind(1.0))then
         call cggev(jobvl,jobvr,n,amat,lda,bmat,ldb,alpha,beta,      &
                    vl,ldvl,vr,ldvr,work,lwork,rwork,info)
      elseif(kd.eq.kind(1.d0))then
         call zggev(jobvl,jobvr,n,amat,lda,bmat,ldb,alpha,beta,      &
                    vl,ldvl,vr,ldvr,work,lwork,rwork,info)
      else
         print*,'Logic error detected in obe_ladder_withDoppler'
         call obe_error
      endif
      if(info.ne.0)then
         print*,'obe_ladder_withDoppler: info is not zero at the'
         print*,'return from c/zggev.',info
         call obe_error
      endif
!
!  Calculate the coefficients of the expansion of the b-vector (the
!  right-hand side) on the generalized eigenvectors corresponding
!  to a finite eigenvalue (i.e., those eigenvectors orthogonal to the
!  null space of c20), and proceed with the calculation.
!
      cnotnul = (0.0_kd,0.0_kd)
      do n = 1,neqs
!  The determination of the null space effected by the following line of
!  code may need to be revisited if inappropriate for the case at hand.
         if(abs(beta(n)).le.1.e-10_kd)cycle
         cmu = alpha(n)/beta(n)
         cvec = matmul(c21,vr(:,n))
         csumn = (0.0_kd,0.0_kd)
         csumd = (0.0_kd,0.0_kd)
         do m = 1,neqs
            csumn = csumn + conjg(vl(m,n))*cb(m)
            csumd = csumd + conjg(vl(m,n))*cvec(m)
         enddo
!
         eta = -cmu/urms
         if(dimag(eta).gt.0.d0)then
            x = dreal(eta)
            y = dimag(eta)
            call obe_Faddeeva(x,y,u,v)
            w = dcmplx(u,v)
            complexint = (0.d0,1.d0)*pi*w
         elseif(dimag(eta).lt.0.d0)then
            x = dreal(eta)
            y = -dimag(eta)
            call obe_Faddeeva(x,y,u,v)
            w = dcmplx(u,v)
            complexint = dconjg((0.d0,1.d0)*pi*w)
         else
            print*,'obe_ladder_withDoppler: Zero imaginary part for eta'
            print*,cmu,eta
            call obe_error
         endif
         coef = complexint/(urms*sqrtpi)
         coef = csumn * coef / csumd
!
         cnotnul = cnotnul + coef*vr(:,n)
         cb = cb - (csumn/csumd)*cvec
!
      enddo
!
!  Contribution of the null space: we solve L_0 y = b - ...
      n = neqs
      nrhs = 1
      amat = c20
      lda = neqs
      ldb = neqs
!
      if(kd.eq.kind(1.0))then
         call cgesv(n,nrhs,amat,lda,ipiv,cb,ldb,info)
      elseif(kd.eq.kind(1.d0))then
         call zgesv(n,nrhs,amat,lda,ipiv,cb,ldb,info)
      else
         print*,'Logic error 2 detected in obe_ladder_withDoppler'
         call obe_error
      endif
!
      if(info.ne.0)then
         print*,'obe_ladder_withDoppler: s/dgesv returns with info = ',info
         call obe_error
      endif
!
!  Solution vector: the sum of the contribution of the null space,
!  which is in cb, with the rest, which is in cnotnul.
      cb = cb + cnotnul
!
!  At this stage, cb contains the correct coherences (within the rate
!  equation approximation) but not the correct populations. Replace the
!  (wrong) populations by those passed to obe_calc_coher:
      do i = nmn,nmx
         cb(nindx(i+noffset,i+noffset)) = popfinal(i)
      enddo
!  Note: Since popfinal does not depend on the velocity, Doppler
!  averaging its elements amounts to multiplying them by 1. Hence,
!  there is no harm with adding them to cb at this stage.
!
!  Arrange the populations and the coherences in the required format
!
      n = 0
      do j = 1,nst
         do i = 1,j
            if(i.ne.j)then
               n = n + 1
               rho_vec(n) = real(cb(nindx(i,j)),kd)
               n = n + 1
               rho_vec(n) = dimag(1.d0*cb(nindx(i,j)))
            else
               n = n + 1
               rho_vec(n) = real(cb(nindx(i,j)),kd)
            endif
         enddo
      enddo
!
!  That's all...
      deallocate(rhs_mat,c2,c20,c21,cb,ipiv,flzeros,amat,bmat,   &
                 vl,vr,alpha,beta,cvec,cnotnul,work,rwork)
!
      end subroutine obe_ladder_withDoppler

!!!!!!!!!!!!!!!!!!!

      end subroutine obe_steadystate_ladder

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine obe_steadystate_onefld(iprep,iDoppler,rhovec_st,urms)
!
!  Returns the steady state density matrix for a single field, with or
!  without Doppler averaging. This averaging, if required, is done
!  semi-analytically, using the Faddeeva function. 
!
!  iprep: The calculation is prepared, including the (possibly)
!  time-consuming diagonalisation of the matrix, before the density
!  matrix is calculated, if iprep = 1. If iprep = 0, the density matrix is
!  calculated using the eigenvalues and eigenvectors calculated in the
!  last call with iprep = 1. 
!
      implicit none
      integer, parameter :: neqs = nst*nst
      integer, parameter :: lwork = 2*(neqs-1)
      integer, parameter :: lrwork = 8*(neqs-1)
      integer, intent(in) :: iDoppler,iprep
      real(kd), dimension (neqs), intent(out) :: rhovec_st
      double precision, optional, intent(in) :: urms
      real(kd), dimension (neqs,neqs) :: rhs_mat
      complex(kd), dimension (neqs,neqs) :: c2
      complex(kd), dimension (neqs-1,neqs-1) :: amat,bmat,c20,c21,vl,vr
      complex(kd), dimension (neqs-1) :: alpha,beta,cvec,cnotnul,cb,cbstore
      real(kd) :: Delta_store
      integer, dimension (nst,nst) :: nindx
      integer :: i,j,k,l,n,m,mre,mim
      integer :: j_state,idetuning
      integer, dimension (neqs-1) :: ipiv
      complex(kd), dimension (lwork) :: work
      real(kd), dimension (lrwork) :: rwork
      integer :: info,lda,ldb,ldvl,ldvr,nrhs
      character*1 :: jobvl,jobvr
      double complex :: csum,cmu,csumn,csumd,coef
      double precision :: ak
      double precision :: pi,sqrtpi
      logical :: fl_init = .true.
      external :: cgetrf,cgetrs,cggev,zgetrf,zgetrs,zggev
      save
!
      if(iprep.ne.0 .and. iprep.ne.1)then
         print*,'obe_steadystate_onefld: iprep must be 0 or 1. ',iprep
         call obe_error
      endif
      if(nst.lt.2 .or. nflds.ne.1)then
         print*,'obe_steadystate_onefld: unsuitable values of nst'
         print*,'or nflds: ',nst,nflds
         call obe_error
      endif
!
      if(iDoppler.eq.1)then
         if(.not.present(urms))then
            print*,'obe_steadystate_onefld: the optional argument'
            print*,'urms must be specified for a calculation'
            print*,'involving a Doppler averaging.'
            call obe_error
         endif
         if(fl_init)then
            pi = 4.d0*datan(1.d0)
            sqrtpi = sqrt(pi)
            fl_init = .false.
         endif
         if(wavel_nm(1).eq.0.0_kd)then  ! wavel_nm is a global variable
            print*,'obe_steadystate_onefld: the wavelength of field 1'
            print*,'must be non-zero for a calculation'
            print*,'involving a Doppler averaging.'
            call obe_error
         endif
         ak = 2.0d0 * pi / (wavel_nm(1)*1.0d-9) ! wave number
      elseif(iDoppler.ne.0)then
         print*,'obe_steadystate_onefld: illegal value of iDoppler.'
         print*,iDoppler
         call obe_error
      endif
!
!  Saving the current value of the detuning
      Delta_store = Delta(1)
!
!  Skip all the preparation if iprep = 0. There is no check that this is
!  reasonable...
!
      if(iprep.eq.0)go to 100
!
!  Prepare the indexes labelling the elements of the density matrix.
!
      nindx = 0
!  Put the population of state 1 first.
      n = 1
      nindx(1,1) = n
      do k = 1,nst
         do l = 1,nst
            if(k.eq.1 .and. l.eq.1)cycle ! Already indexed.
            n = n + 1
            nindx(k,l) = n
         enddo
      enddo
      if(n.ne.neqs)then
         print*,'obe_steadystate_onefld: logic error.',n,neqs
         call obe_error
      endif
!
!  Big loop. At the first run, calculate the relevant arrays for a zero
!  detuning, at the second for an angular frequency detuning of 1e-6 s^-1.
!
      do idetuning = 0,1
!
         if(idetuning.eq.0)then
            Delta(1) = 0.0_kd
         else
            Delta(1) = 1.0_kd  ! Set to 1 MHz in frequency.
         endif
!
!  Get the matrix representing the right-hand sides in the form
!  of a complex matrix acting on populations and complex coherences.
!  Starting from its real representation returned by obe_set_rhsmat 
!  ensures that it will have already been modified as necessary if 
!  calculation is done within the weak probe approximation.
!  There is no reason to go through obe_set_rhsmat if the weak
!  probe approximation is not invoked, though, in which case
!  we get the necessary complex matrix directly from ldbl.
!  Programming note: the first method yields a matrix which is the
!  transpose of that obtained in the second method. This is taken care
!  of later on, where necessary.
         if(fl_wkprbapp)then
            call obe_set_rhsmat(rhs_mat)
            call obe_transform_rhs_mat(nindx,rhs_mat,c2)
         else
            call ldbl_set_rhsmat(rhs_cmat=c2)
         endif
!
!  Set the population of state 1 to 1, hence the corresponding row and
!  the corresponding column of c2 go and the column vector cb is
!  the negative of the latter... 
!
         j_state = 1
         if(nindx(j_state,j_state).ne.1)then
            print*,'obe_steadystate_onefld: j_state must be 1'
            print*,'in the current state of development of the code.'
            call obe_error
         endif
         cb(1:neqs-1) = -c2(2:neqs,1)
!
!  Subtract the relevant quantity from the entries of c2 corresponding
!  to populations.
!
!  The following code applies only to the case where
!  nindx(j_state,j_state).eq.1. This is first (re)checked.
         if(nindx(j_state,j_state).ne.1)then
            print*,'obe_steadystate_onefld: j_state must be 1'
            print*,'in the current state of development of the code.'
            call obe_error
         endif
!
         do k = 1,nst
            if(k.eq.j_state)cycle
            c2(2:neqs,nindx(k,k)) = c2(2:neqs,nindx(k,k)) + cb(1:neqs-1)
         enddo
!
!  Store the result
         if(idetuning.eq.0)c20(1:neqs-1,1:neqs-1) = c2(2:neqs,2:neqs)
         if(idetuning.eq.1)c21(1:neqs-1,1:neqs-1) = c2(2:neqs,2:neqs)
!
!  End of the big loop
!
      enddo
!
      c21 = c21 - c20
!
!  Next step: Find the null space and the generalized eigenvalues.
!
      jobvl = 'V'
      jobvr = 'V'
      n = neqs-1
      lda = neqs - 1
      ldb = neqs - 1
      ldvl = neqs - 1
      ldvr = neqs - 1
      amat = c20
      bmat = c21
      if(kd.eq.kind(1.0))then
         call cggev(jobvl,jobvr,n,amat,lda,bmat,ldb,alpha,beta,      &
                    vl,ldvl,vr,ldvr,work,lwork,rwork,info)
      elseif(kd.eq.kind(1.d0))then
         call zggev(jobvl,jobvr,n,amat,lda,bmat,ldb,alpha,beta,      &
                    vl,ldvl,vr,ldvr,work,lwork,rwork,info)
      else
         print*,'Logic error detected in obe_steadystate_onefld'
         call obe_error
      endif
      if(info.ne.0)then
         print*,'obe_steadystate_onefld: info is not zero at the'
         print*,'return from zggev.',info
         call obe_error
      endif
!
!  Contribution of the null space: we need to solve L_0 y = b - ...
!  We prepare this calculation my LU-factorizing L0:
      n = neqs - 1
      amat = c20
      lda = neqs - 1
      if(kd.eq.kind(1.0))then
         call cgetrf(n,n,amat,lda,ipiv,info)
      elseif(kd.eq.kind(1.d0))then
         call zgetrf(n,n,amat,lda,ipiv,info)
      else
         print*,'Logic error 2 detected in obe_steadystate_onefld'
         call obe_error
      endif
      if(info.ne.0)then
         print*,'obe_steadystate_onefld: info is not zero at the'
         print*,'return from cgetrf or zgetrf.',info
         call obe_error
      endif
!
      cbstore = cb
!
!  The calculation with iprep = 0 resumes here at the 100 continue
!  below.
!
 100  continue 
!
      cb = cbstore
!
!  Calculate the coefficients of the expansion of the b-vector (the
!  right-hand side) on the generalized eigenvectors corresponding
!  to a finite eigenvalue (i.e., those eigenvectors orthogonal to the
!  null space of c20), and proceed with the calculation.
      cnotnul = (0.0_kd,0.0_kd)
      do n = 1,neqs-1
!  The determination of the null space effected by the following line of
!  code may need to be revisited if inappropriate for the case at hand.
         if(cdabs(1.d0*beta(n)).le.1.d-10)cycle
         cmu = alpha(n)/beta(n)
         if(dimag(cmu).eq.0.d0)then
            print*,'obe_steadystate_onefld: cmu is real, not complex... ',cmu
            call obe_error
         endif
         cvec = matmul(c21,vr(:,n))
         csumn = (0.0_kd,0.0_kd)
         csumd = (0.0_kd,0.0_kd)
         do m = 1,neqs-1
            csumn = csumn + conjg(vl(m,n))*cb(m)
            csumd = csumd + conjg(vl(m,n))*cvec(m)
         enddo
         call obe_calc_coef   ! This internal subroutine calculates 
                              ! the variable coef used below.
         cnotnul = cnotnul + coef*vr(:,n)
         cb = cb - (csumn/csumd)*cvec
      enddo
!
!  Contribution of the null space: we solve L_0 y = b - ...
      n = neqs - 1
      nrhs = 1
      lda = neqs - 1
      ldb = neqs - 1
      if(kd.eq.kind(1.0))then
         call cgetrs('N',n,nrhs,amat,lda,ipiv,cb,ldb,info)
      elseif(kd.eq.kind(1.d0))then
         call zgetrs('N',n,nrhs,amat,lda,ipiv,cb,ldb,info)
      else
         print*,'Logic error 3 detected in obe_steadystate_onefld'
         call obe_error
      endif
      if(info.ne.0)then
         print*,'obe_steadystate_onefld: info is not zero at the'
         print*,'return from cgetrs or zgetrs.',info
         call obe_error
      endif
!
!  Solution vector: the sum of the contribution of the null space,
!  which is in cb, with the rest, which is in cnotnul.
      cb = cb + cnotnul
!
!  Arrange the populations and the coherences in the required format 
!
      n = 1
      rhovec_st(1) = 1.0_kd
      do j = 1,nst
         do i = 1,j
            if(i.eq.1 .and. j.eq.1)then
               cycle
            elseif(i.ne.j)then
               n = n + 1
               if(fl_wkprbapp)then
                  rhovec_st(n) = real(cb(nindx(i,j)-1),kd)
               else
                  rhovec_st(n) = real(cb(nindx(j,i)-1),kd)
               endif
               n = n + 1
               if(fl_wkprbapp)then
                  rhovec_st(n) = dimag(1.d0*cb(nindx(i,j)-1))
               else
                  rhovec_st(n) = dimag(1.d0*cb(nindx(j,i)-1))
               endif
            else
               n = n + 1
               if(fl_wkprbapp)then
                  rhovec_st(n) = real(cb(nindx(i,j)-1),kd)
               else
                  rhovec_st(n) = real(cb(nindx(j,i)-1),kd)
               endif
               rhovec_st(1) = rhovec_st(1) - rhovec_st(n)
            endif
         enddo
      enddo
!
!  Reset the detuning to its original value and exit.
      Delta(1) = Delta_store
!
      contains

!!!!!!!!!!!!!!!!!!

      subroutine obe_calc_coef
!
      implicit none
      double complex :: eta,w,complexint
      double precision :: x,y,u,v
!
      if(iDoppler.eq.0)then
!
         coef = csumn/((cmu+(1.d0,0.d0)*Delta_store)*csumd)
!
      else
!
         eta = (cmu+(1.d0,0.d0)*Delta_store)/(ak*urms)
!  Multiply by 2pi times 10^6 because cmu and Delta_store are expressed
!  in MHz (cmu because Delta(1) was set to 1 MHz, not to 2 pi 1e6 s^-1).
         eta = 2.d0*pi * 1.d6 * eta
!
         if(dimag(eta).gt.0.d0)then
            x = dreal(eta)
            y = dimag(eta)
            call obe_Faddeeva(x,y,u,v)
            w = dcmplx(u,v)
            complexint = (0.d0,1.d0)*pi*w
         elseif(dimag(eta).lt.0.d0)then
            x = dreal(eta)
            y = -dimag(eta)
            call obe_Faddeeva(x,y,u,v)
            w = dcmplx(u,v)
            complexint = dconjg((0.d0,1.d0)*pi*w)
         else
            print*,'obe_calc_coef: Zero imaginary part for eta'
            print*,cmu,eta
            call obe_error
         endif
 
         coef = -complexint/(ak*urms*sqrtpi)
!  Multiply by 2pi times 1d6 because the value of coef obtained at the
!  previous line amounts to an average value of 1/(Delta + cmu) where
!  the demominator is expressed in 1/s rather than in MHz.
         coef = coef * 2.d0*pi*1.d6
         coef = csumn * coef / csumd
!
      endif
!
      end subroutine obe_calc_coef

!!!!!!!!!!!!!!!!!!

      end subroutine obe_steadystate_onefld

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine obe_steadystate_onefld_weakprb(iDoppler,popinit,rhovec_st,urms)
!
!  Returns the steady state density matrix for a single field in the
!  weak probe approximation, with or without Doppler averaging.
!  This averaging, if required, is done semi-analytically, using the
!  Faddeeva function. 
!
!  It is assumed that all the non-zero elements of popinit are
!  populations of low energy states - i.e., of states that are coupled
!  to states of higher energy by the field.
!
      implicit none
      integer, parameter :: neqs = nst*nst
      integer, intent(in) :: iDoppler
      real(kd), dimension (neqs), intent(out) :: rhovec_st
      double precision, dimension(nmn:nmx), intent(in) :: popinit
      double precision, optional, intent(in) :: urms
      double complex :: cfact,coher,w
      double precision, dimension(nst) :: Gammatot_decay_f
      double precision :: ak,delta_trans,denom,gamma_tr,u,v,x,y
      double precision, save :: conv,pi,sqrtpi
      integer :: i,j,m,mre,mim
      logical, save :: fl_init = .true.
!
      if(fl_init)then
         pi = 4.d0*datan(1.d0)
         conv = 1.0e-6_kd/(2.d0*pi*hbar)
         sqrtpi = sqrt(pi)
         fl_init = .false.
      endif
!
!  The following prevents this subroutine to work if obe_setcsts
!  has been called with iweakprb = 0. It is not strictly necessary,
!  since nothing in the calculations done here depend on the value of
!  iweakprb. However, allowing this subroutine to perfom weak probe
!  calculations even when iweakpb = 0 is potentially confusing, hence
!  this restriction.
      if(.not.fl_wkprbapp)then
         print*,'obe_steadystate_onefld_weakprb can be used only for'
         print*,'calculations in the weak probe approximation.'
         call obe_error
      endif
!
      if(nst.lt.2 .or. nflds.ne.1)then
         print*,'obe_steadystate_onefld_weakprb: unsuitable values of nst'
         print*,'or nflds: ',nst,nflds
         call obe_error
      endif
!
      if(iDoppler.eq.1)then
         if(.not.present(urms))then
            print*,'obe_steadystate_onefld_weakprb: the optional argument'
            print*,'urms must be specified for a calculation'
            print*,'involving a Doppler averaging.'
            call obe_error
         endif
         if(wavel_nm(1).eq.0.0_kd)then  ! wavel_nm is a global variable
            print*,'obe_steadystate_onefld_weakprb: the wavelength of field 1'
            print*,'must be non-zero for a calculation'
            print*,'involving a Doppler averaging.'
            call obe_error
         endif
         ak = 2.0d0 * pi / (wavel_nm(1)*1.0d-9) ! wave number
         denom = 1.d0/(1.d-6*ak*urms/(2.d0*pi))
      elseif(iDoppler.ne.0)then
         print*,'obe_steadystate_onefld_weakprb: illegal value of iDoppler.'
         print*,iDoppler
         call obe_error
      endif
!
      if(n_collapse_ops .ne. -1)then
         print*,'Using obe_steadystate_onefld_weakprb is incompatible'
         print*,'with specifying collapse operators explicitly.'
         call obe_error
      endif
!
!  Calculate the total spontaneous decay rate of the upper states
      Gammatot_decay_f = 0.d0
      do j = 1,nst
         do i = 1,nst
            Gammatot_decay_f(j) = Gammatot_decay_f(j) + Gamma_decay_f(i,j,1)
         enddo
         if(Gammatot_decay_f(j).ne.0.d0 .and. popinit(j-noffset).ne.0.d0)then
            print*,'obe_steadystate_onefld_weakprb finds a non-zero decay rate'
            print*,'for state ',j-noffset,' which has also a non-zero initial'
            print*,'population. This is inconsistent since only the lower'
            print*,'energy states should be populated.'
            call obe_error
         endif
      enddo
!
      if(fl_Rabi)then
         cfact = (0.d0,1.d0)/2.d0
      else
         cfact = (0.d0,1.d0)*campl(1)*conv/2.d0
      endif
!
      rhovec_st = 0.0_kd
!
      do j = 1,nst
         call ldbl_pop_index(j,m)
         rhovec_st(m) = popinit(j-noffset)
         if(popinit(j-noffset).le.0.d0)cycle
         do i = 1,nst
            if(fl_Rabi)then
               if(Rabifreq(i,j,1).eq.(0.d0,0.d0))cycle
            else
               if(dipmom(i,j,1).eq.(0.d0,0.d0))cycle
            endif
!
            if(popinit(i-noffset).gt.0.d0)then
               print*,'obe_steadystate_onefld_weakprb: the populations'
               print*,'of the upper states must be zero.'
               print*,i-noffset,popinit(i-noffset)
               call obe_error
            endif
!  Homogeneous broadening for this pair of states:
            gamma_tr = Gammatot_decay_f(i)/2.d0
            if(fl_add_dephas)gamma_tr = gamma_tr + add_dephas_f(i,j)
!
            delta_trans = energ_f(j) - energ_f(i) + Delta(1)
!
            if(iDoppler.eq.1)then
!  Coherence with Doppler averaging:
               x = delta_trans*denom
               y = gamma_tr*denom
               call obe_Faddeeva(x,y,u,v)
               w = dcmplx(u,v)
               if(fl_Rabi)then
                  coher = cfact*Rabifreq(i,j,1)*popinit(j-noffset)*    &
                                                    sqrtpi*w*denom
               else
                  coher = cfact*dipmom(i,j,1)*popinit(j-noffset)*      &
                                                    sqrtpi*w*denom
               endif
            else
!  Coherence without Doppler averaging:
               if(fl_Rabi)then
                  coher = cfact*Rabifreq(i,j,1)*popinit(j-noffset)/    &
                            (gamma_tr-(0.d0,1.d0)*delta_trans)
               else
                  coher = cfact*dipmom(i,j,1)*popinit(j-noffset)/      &
                            (gamma_tr-(0.d0,1.d0)*delta_trans)
               endif
            endif
!
            if(j.gt.i)then
               call ldbl_coher_index(i,j,mre,mim)
               rhovec_st(mre) = dble(coher)
               rhovec_st(mim) = dimag(coher)
            elseif(i.gt.j)then
               call ldbl_coher_index(j,i,mre,mim)
               rhovec_st(mre) = dble(coher)
               rhovec_st(mim) = -dimag(coher)
            else
               print*,'obe_steadystate_onefld_weakprb founds a non-zero'
               print*,'diagonal dipole moment.',i
               call obe_error
            endif
!
          enddo
      enddo
!
      end subroutine obe_steadystate_onefld_weakprb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine obe_steadystate_onefld_powerbr(iDoppler,popinit,rhovec_st,urms)
!
!  Returns the steady state density matrix for a single field in the
!  weak probe approximation, with or without Doppler averaging, with
!  power broadening included (but not optical pumping). Doppler
!  averaging, if required, is done semi-analytically, using the
!  Faddeeva function. The only difference with the subroutine
!  obe_steadystate_onefld_weakprb is the power broadening element.
!
!  It is assumed that all the non-zero elements of popinit are
!  populations of low energy states - i.e., of states that are coupled
!  to states of higher energy by the field.
!
      implicit none
      integer, parameter :: neqs = nst*nst
      integer, intent(in) :: iDoppler
      real(kd), dimension (neqs), intent(out) :: rhovec_st
      double precision, dimension(nmn:nmx), intent(in) :: popinit
      double precision, optional, intent(in) :: urms
      double complex :: cfact,coher,rabi,w
      double precision, dimension(nst) :: Gammatot_decay_f
      double precision :: ak,delta_trans,denom,gamma_tr,r,u,v,x,y
      double precision, save :: conv,pi,sqrtpi
      integer :: i,j,m,mre,mim
      logical, save :: fl_init = .true.
!
      if(fl_init)then
         pi = 4.d0*datan(1.d0)
         conv = 1.0e-6_kd/(2.d0*pi*hbar)
         sqrtpi = sqrt(pi)
         fl_init = .false.
      endif
!
!  The following prevents this subroutine to work if obe_setcsts
!  has been called with iweakprb = 0. It is not strictly necessary,
!  since nothing in the calculations done here depend on the value of
!  iweakprb. However, allowing this subroutine to perfom weak probe
!  calculations even when iweakpb = 0 is potentially confusing, hence
!  this restriction.
      if(.not.fl_wkprbapp)then
         print*,'obe_steadystate_onefld_powerbr can be used only for'
         print*,'calculations in the weak probe approximation.'
         call obe_error
      endif
!
      if(nst.lt.2 .or. nflds.ne.1)then
         print*,'obe_steadystate_onefld_powerbr: unsuitable values of nst'
         print*,'or nflds: ',nst,nflds
         call obe_error
      endif
!
      if(iDoppler.eq.1)then
         if(.not.present(urms))then
            print*,'obe_steadystate_onefld_powerbr: the optional argument'
            print*,'urms must be specified for a calculation'
            print*,'involving a Doppler averaging.'
            call obe_error
         endif
         if(wavel_nm(1).eq.0.0_kd)then  ! wavel_nm is a global variable
            print*,'obe_steadystate_onefld_powerbr: the wavelength of field 1'
            print*,'must be non-zero for a calculation'
            print*,'involving a Doppler averaging.'
            call obe_error
         endif
         ak = 2.0d0 * pi / (wavel_nm(1)*1.0d-9) ! wave number
         denom = 1.d0/(1.d-6*ak*urms/(2.d0*pi))
      elseif(iDoppler.ne.0)then
         print*,'obe_steadystate_onefld_powerbr: illegal value of iDoppler.'
         print*,iDoppler
         call obe_error
      endif
!
      if(n_collapse_ops .ne. -1)then
         print*,'Using obe_steadystate_onefld_powerbr is incompatible'
         print*,'with specifying collapse operators explicitly.'
         call obe_error
      endif
!
!  Calculate the total spontaneous decay rate of the upper states
      Gammatot_decay_f = 0.d0
      do j = 1,nst
         do i = 1,nst
            Gammatot_decay_f(j) = Gammatot_decay_f(j) + Gamma_decay_f(i,j,1)
         enddo
         if(Gammatot_decay_f(j).ne.0.d0 .and. popinit(j-noffset).ne.0.d0)then
            print*,'obe_steadystate_onefld_powerbr finds a non-zero decay rate'
            print*,'for state ',j-noffset,' which has also a non-zero initial'
            print*,'population. This is inconsistent since only the lower'
            print*,'energy states should be populated.'
            call obe_error
         endif
      enddo
!
      if(fl_Rabi)then
         cfact = (0.d0,1.d0)/2.d0
      else
         cfact = (0.d0,1.d0)*campl(1)*conv/2.d0
      endif
!
      rhovec_st = 0.0_kd
!
      do j = 1,nst
         call ldbl_pop_index(j,m)
         rhovec_st(m) = popinit(j-noffset)
         if(popinit(j-noffset).le.0.d0)cycle
         do i = 1,nst
            if(fl_Rabi)then
               if(Rabifreq(i,j,1).eq.(0.d0,0.d0))cycle
            else
               if(dipmom(i,j,1).eq.(0.d0,0.d0))cycle
            endif
!
            if(popinit(i-noffset).gt.0.d0)then
               print*,'obe_steadystate_onefld_powerbr: the populations'
               print*,'of the upper states must be zero.'
               print*,i-noffset,popinit(i-noffset)
               call obe_error
            endif
!  Homogeneous broadening for this pair of states:
            gamma_tr = Gammatot_decay_f(i)/2.d0
            if(fl_add_dephas)gamma_tr = gamma_tr + add_dephas_f(i,j)
            if(fl_Rabi)then
               rabi = Rabifreq(i,j,1)
            else
               rabi = dipmom(i,j,1)*campl(1)*conv
            endif
            r = dsqrt(gamma_tr**2+cdabs(rabi)**2*gamma_tr/    &
                                        Gammatot_decay_f(i))
!
            delta_trans = energ_f(j) - energ_f(i) + Delta(1)
!
            if(iDoppler.eq.1)then
!  Coherence with Doppler averaging:
               x = delta_trans*denom
               y = r*denom
               call obe_Faddeeva(x,y,u,v)
               w = dcmplx(u,v)
               coher = (gamma_tr+r)/(2.d0*r)*sqrtpi*w*denom
               x = -delta_trans*denom
               y = r*denom
               call obe_Faddeeva(x,y,u,v)
               w = dcmplx(u,v)
               coher = coher + (gamma_tr-r)/(2.d0*r)*sqrtpi*w*denom
            else
!  Coherence without Doppler averaging:
               coher = (gamma_tr+r)/(2.d0*r)/(r-(0.d0,1.d0)*delta_trans) + &
                       (gamma_tr-r)/(2.d0*r)/(r+(0.d0,1.d0)*delta_trans)
            endif
            if(fl_Rabi)then
               coher = cfact*Rabifreq(i,j,1)*popinit(j-noffset)*coher
            else
               coher = cfact*dipmom(i,j,1)*popinit(j-noffset)*coher
            endif
!
            if(j.gt.i)then
               call ldbl_coher_index(i,j,mre,mim)
               rhovec_st(mre) = dble(coher)
               rhovec_st(mim) = dimag(coher)
            elseif(i.gt.j)then
               call ldbl_coher_index(j,i,mre,mim)
               rhovec_st(mre) = dble(coher)
               rhovec_st(mim) = -dimag(coher)
            else
               print*,'obe_steadystate_onefld_powerbr founds a non-zero'
               print*,'diagonal dipole moment.',i
               call obe_error
            endif
!
          enddo
      enddo
!
      end subroutine obe_steadystate_onefld_powerbr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end module obe

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine ext_setsys(h_mat,decrates_mat,t2dephas_mat,ncoll)
!
!  For use by ldbl_set_rhsmat (part of the ldbl module).
!  ext_setsys simply gets the arrays h_mat, decrates_mat and t2dephas_mat
!  from obe_setsys (part of the obe module) and pass them to ldbl.
!
      use general_settings
      use obe
      implicit none
      complex(kd), dimension (nst,nst), intent(out) :: h_mat
      real(kd), dimension (:,:,:), allocatable, intent(out) :: decrates_mat
      complex(kd), dimension (nst,nst), intent(out) :: t2dephas_mat
      integer, intent(out) :: ncoll
!
      call obe_setsys(h_mat,decrates_mat,t2dephas_mat,ncoll)
!
      end subroutine ext_setsys
