      module mbe
!
!  The program units contained in this module integrate the Maxwell-Bloch
!  equations. It makes use of obe, the optical-Bloch equations module.
!
      use general_settings
      use ldbl
      use obe
      implicit none
!
!  The following statement makes the contents of mbe "invisible" to other
!  parts of the programs, with the exceptions of the variables and subprograms
!  declared below as public.
!
      private
!
!------------------------------------------------------------------------------
!
!  The following subprograms contained in mbe may be called from outside the
!  module:
!
      public mbe_propagate_2, mbe_set_envlp, mbe_set_tdfields_A,        &
             mbe_set_tdfields_B, mbe_tdint_2, mbe_propagate_1,          &
             mbe_tdint_1
!
!  All the other subprograms contained in mbe can be accessed only from mbe.
!
!------------------------------------------------------------------------------
!
!  The following variables will be "known" by all the subprograms
!  contained in this module but are not directly accessible from the outside.
!
!  It is assumed within mbe that Field_p_re, Field_p_im, Field_c_re and
!  Field_c_im are the real and imaginary parts of the electric field 
!  components of, respectively, field 1 (Field_p) and field 2
!  (Field_c).
      real(kd) :: Field_p_re, Field_p_im, Field_c_re, Field_c_im
!
      real(kd), dimension(:), allocatable :: Field_p_re_vec
      real(kd), dimension(:), allocatable :: Field_p_im_vec
      real(kd), dimension(:), allocatable :: Field_c_re_vec
      real(kd), dimension(:), allocatable :: Field_c_im_vec
!
!  Array initialized by mbe_set_rhsmat_0 and used by mbe_set_rhsmat
      real(kd) :: rhs_mat_0(nst*nst,nst*nst)
!
!  Variables used in the time integration (allocated by mbe_propagate_2
!  or mbe_tdint_2 or their 1-field version):
      integer :: ktime
      real(kd), dimension(:), allocatable :: Fd_pre2
      real(kd), dimension(:), allocatable :: Fd_pim2
      real(kd), dimension(:), allocatable :: Fd_cre2
      real(kd), dimension(:), allocatable :: Fd_cim2
!
!  Variables allocated by mbe_set_envlp:
      double precision :: tw_p,t0_p,t1_p,tw_c,t0_c,t1_c
      character*2 :: pulse_type_p='xx',pulse_type_c='xx'
      logical :: force0_p,force0_c
      logical :: fl_p_env_notset=.true.,fl_c_env_notset=.true.
!
!  Variables allocated and/or set by mbe_tdfields_A or mbe_tdfields_B:
      integer :: n_time_steps
      double complex, dimension(:), allocatable :: Field_p_vec, Field_c_vec
      real(kd), dimension(:), allocatable :: t_mesh
      double complex :: cfield0_p,cfield0_c
      logical :: fl_notinit_fields = .true.
!
!  Variables relevant for the Doppler integration
      real(kd) :: current_vel, urms
      real(kd), dimension(:), allocatable :: vmesh, fMvweight
      integer :: n_v_values
!
!  Variables relevant for calculations in the rate equation approximation
      real(kd), dimension (:,:), allocatable :: arhs
      integer, dimension (:), allocatable :: ipiv_r
!
!  Variables relevant for the output of the results
      integer, dimension (:,:), allocatable :: iunits_arr
      logical :: fl_outputfiles = .false.
!
!  Logical variables used to check or control the job flow
      logical :: fl_Doppler = .false.
      logical :: fl_init_ommats = .true.
      logical :: fl_rhsmat0_not_init = .true.
      logical :: fl_Doppler_notinit = .true.
      logical :: fl_twofields = .true.
      logical :: fl_spline = .true.
!
!  Variables used to communicate with the dop853 ODE solver
      double precision, save :: rtol_dop853,atol_dop853
      logical, save :: fl_set_tol_dop853 = .false.
!
      contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine mbe_error
!
!  Control may be passed to this subprogram from elsewhere within mbe
!  following detection of an error.
!
!  At the moment mbe_error simply stops the execution without saving files.
!  The user may want to replace this by a more sophisticated error handling
!  procedure.
!
      stop 'Error in the mbe module...'
!
      end subroutine mbe_error

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine mbe_set_outputfiles
!
!  Set the global variable fl_outputfiles, and also, if fl_outputfiles is
!  .true., the global variable iunits_arr, by obtaining this information
!  from the obe module through a call to obe_get_iunits. The array
!  iunits_arr is allocatable. It is allocated by obe_get_iunits only if
!  fl_outputfiles is .true., and won't therefore be allocated if
!  fl_outputfiles is .false.
!
      implicit none
!
      call obe_get_iunits(fl_outputfiles,iunits_arr)
!
      end subroutine mbe_set_outputfiles

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine mbe_set_rhsmat_0
!
!  Calculate the matrix representing the right-hand sides of the optical
!  Bloch equations for zero fields, rhs_mat_0. This array is global within
!  mbe and is used by mbe_set_rhsmat.
!
!  It is assumed that mbe_set_rhsmat_0 has been called at an appropriate
!  point, prior the use of rhsmat_0, and with the appropriate settings passed
!  on to the obe module. mbe_set_rhsmat_0 is called by mbe_propagate_2 and
!  mbe_tdint_2 at the start of the calculation.
!
!  mbe_set_rhsmat_0 also checks that two and only two fields are set in obe.
!
      implicit none
      integer :: nflds_in
      double complex :: cfp_store,cfc_store
      double complex, dimension(:), allocatable :: cf
!
!  Check that everything is as expected
      call obe_get_campl(nflds_in,cf)
      if(fl_twofields .and. nflds_in.ne.2)then
         print*,'mbe_set_rhsmat_0 detects an error: nflds is not 2 in obe'
         print*,'but mbe is used for a two-field calculation.'
         call mbe_error
      endif
      if(.not.fl_twofields .and. nflds_in.ne.1)then
         print*,'mbe_set_rhsmat_0 detects an error: nflds is not 1 in obe'
         print*,'but mbe is used for a single field calculation.'
         call mbe_error
      endif
!
!  Set the complex field amplitudes of the two fields to zero for calculating
!  rhsmat_0:
      if(allocated(cf))then
         cfp_store = cf(1)
         cf(1) = (0.d0,0.d0)
         call obe_reset_campl(1,cf(1))
         if(fl_twofields)then
            cfc_store = cf(2)
            cf(2) = (0.d0,0.d0)
            call obe_reset_campl(2,cf(2))
         endif
      else
         print*,'mbe_set_rhsmat_0: the array fields is not allocated...'
         call mbe_error
      endif
!
!  Calculate rhs_mat_0
      call ldbl_set_rhsmat(rhs_mat_0)
!
!  Reset the complex field amplitudes of the two fields to their original 
!  values:
      cf(1) = cfp_store
      call obe_reset_campl(1,cf(1))
      if(fl_twofields)then
         cf(2) = cfc_store
         call obe_reset_campl(2,cf(2))
      endif
!
!  All done...
      fl_rhsmat0_not_init = .false.
!
      end subroutine mbe_set_rhsmat_0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine mbe_set_rhsmat2(rhs_mat,k)
!
!  This subroutine comes in two versions: mbe_set_rhsmat2 for two fields
!  and mbe_set_rhsmat1 for a single field calculation. This is mbe_set_rhsmat2.
!
!  Calculate the matrix representing the right-hand sides of the optical
!  Bloch equations, rhs_mat, using
!    rhs_mat(:,:) = omeg0_mat(:,:) +
!      sum_j [ Re(omega_j) omegrl_mat(:,:,j) + Im(omega_j) omegim_mat(:,:,j)],
!  or equivalently for a Doppler averaging.
!
!  The arrays appearing in this equation are calculated within obe by
!  obe_set_ommats. They depend on the detunings, on the dephasing rates and
!  possibly also on other parameters of the field (apart from the complex field
!  amplitudes) and must be re-calculated if the value of any of these is
!  changed.
!
!  A value of .true. for the flag fl_init_ommats is the signal to invoke
!  obe_set_ommats to recalculate these arrays. (This flag is reset to .false.
!  by this subroutine.)
!
!  If k.eq.0, the zero-field array is taken to be rhs_mat_0, not omeg0_mat.
!
!  Note: Only the probe field is taken into account in part of the calculation
!  when the global logical variable fl_twofields is .false.
!
      implicit none
      integer, parameter :: neqs = nst*nst
      integer, parameter :: nflds = 2
      integer, intent(in) :: k
      real(kd), dimension (neqs,neqs), save  :: omeg0_mat
      real(kd), dimension (:,:,:), allocatable, save  :: omeg0_v_mat
      real(kd), dimension (neqs,neqs,nflds), save  :: omegrl_mat
      real(kd), dimension (neqs,neqs,nflds), save  :: omegim_mat
      real(kd), dimension (neqs,neqs), intent(out) :: rhs_mat
      integer :: istat
      logical, save :: fl_need_init_varrays = .true.
!
      if(.not.fl_twofields)then
         print*,'mbe_set_rhsmat2 should not be used in a 1-field calculation.'
         call mbe_error
      endif
!
      if(fl_init_ommats)then
!  Get the necessary arrays and reset fl_init_ommats and fl_need_init_varrays.
!  First, signal that omeg0_v_mat should also be (re)-initialized if needed.
         fl_need_init_varrays = .true.
         if(fl_Doppler)then
            if(allocated(omeg0_v_mat))deallocate(omeg0_v_mat)
            allocate(omeg0_v_mat(neqs,neqs,1),stat=istat)
            if(istat.gt.0)then
               print*,'mbe_set_rhsmat: Error at the allocation of omeg0_v_mat'
               call mbe_error
            endif
            call obe_set_ommats(omeg0_mat,omegrl_mat,omegim_mat,omeg0_v_mat,'v')
            fl_need_init_varrays = .false.
         else
            call obe_set_ommats(omeg0_mat,omegrl_mat,omegim_mat)
         endif
         fl_init_ommats = .false.
      endif
!
!  Calculate the first order matrices, changing the fields one by one
      if(k.eq.0)then
         if(fl_rhsmat0_not_init)then
            print*,'mbe_set_rhsmat2: Error.'
            print*,'Attempt at using rhs_mat_0 before initialization.'
            call mbe_error
         endif
         rhs_mat = rhs_mat_0 + &
            Field_p_re*omegrl_mat(:,:,1) + &
            Field_c_re*omegrl_mat(:,:,2) + &
            Field_p_im*omegim_mat(:,:,1) + &
            Field_c_im*omegim_mat(:,:,2) 
      else
         rhs_mat = omeg0_mat + &
            Field_p_re*omegrl_mat(:,:,1) + &
            Field_c_re*omegrl_mat(:,:,2) + &
            Field_p_im*omegim_mat(:,:,1) + &
            Field_c_im*omegim_mat(:,:,2) 
      endif
!
! Take the Doppler shift into account, if necessary, for the velocity defined
! by the value of the global variable current_vel.
!
! IMPORTANT: It is here assumed that the correction is the same whether
! k = 0 or not. I.e., it is assumed that the v-dependent terms in the array
! representing the rhs of the optical Bloch equations are the same whether
! this array represents the full system (k > 0) or the simplified system used
! when k = 0. This assumption is justified in so far the latter differs from
! the full system only through the neglect of decay to a sink state, but may
! be incorrect in other cases.
!
      if(fl_Doppler)then
         if(fl_need_init_varrays)then
            print*,'mbe_set_rhsmat2 has detected a logic error'
            call mbe_error
         endif
         rhs_mat = rhs_mat + current_vel*omeg0_v_mat(:,:,1)
      endif
!
      end subroutine mbe_set_rhsmat2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine mbe_set_rhsmat1(rhs_mat,k)
!
!  The same as mbe_set_rhsmat2 but for a single field.
!
      implicit none
      integer, parameter :: neqs = nst*nst
      integer, parameter :: nflds = 1
      integer, intent(in) :: k
      real(kd), dimension (neqs,neqs), save  :: omeg0_mat
      real(kd), dimension (:,:,:), allocatable, save  :: omeg0_v_mat
      real(kd), dimension (neqs,neqs,nflds), save  :: omegrl_mat
      real(kd), dimension (neqs,neqs,nflds), save  :: omegim_mat
      real(kd), dimension (neqs,neqs), intent(out) :: rhs_mat
      integer :: istat
      logical, save :: fl_need_init_varrays = .true.
!
      if(fl_twofields)then
         print*,'mbe_set_rhsmat1 should not be used in a 2-field calculation.'
         call mbe_error
      endif
!
      if(fl_init_ommats)then
!  Get the necessary arrays and reset fl_init_ommats and fl_need_init_varrays.
!  First, signal that omeg0_v_mat should also be (re)-initialized if needed.
         fl_need_init_varrays = .true.
         if(fl_Doppler)then
            if(allocated(omeg0_v_mat))deallocate(omeg0_v_mat)
            allocate(omeg0_v_mat(neqs,neqs,1),stat=istat)
            if(istat.gt.0)then
               print*,'mbe_set_rhsmat1: Error at the allocation of omeg0_v_mat'
               call mbe_error
            endif
            call obe_set_ommats(omeg0_mat,omegrl_mat,omegim_mat,omeg0_v_mat,'v')
            fl_need_init_varrays = .false.
         else
            call obe_set_ommats(omeg0_mat,omegrl_mat,omegim_mat)
         endif
         fl_init_ommats = .false.
      endif
!
!  Calculate the first order matrices, changing the fields one by one
      if(k.eq.0)then
         if(fl_rhsmat0_not_init)then
            print*,'mbe_set_rhsmat1: Error.'
            print*,'Attempt at using rhs_mat_0 before initialization.'
            call mbe_error
         endif
         rhs_mat = rhs_mat_0 + &
            Field_p_re*omegrl_mat(:,:,1) + &
            Field_p_im*omegim_mat(:,:,1)
      else
         rhs_mat = omeg0_mat + &
            Field_p_re*omegrl_mat(:,:,1) + &
            Field_p_im*omegim_mat(:,:,1)
      endif
!
! Take the Doppler shift into account, if necessary, for the velocity defined
! by the value of the global variable current_vel.
!
! IMPORTANT: It is here assumed that the correction is the same whether
! k = 0 or not. I.e., it is assumed that the v-dependent terms in the array
! representing the rhs of the optical Bloch equations are the same whether
! this array represents the full system (k > 0) or the simplified system used
! when k = 0. This assumption is justified in so far the latter differs from
! the full system only through the neglect of decay to a sink state, but may
! be incorrect in other cases.
!
      if(fl_Doppler)then
         if(fl_need_init_varrays)then
            print*,'mbe_set_rhsmat1 has detected a logic error'
            call mbe_error
         endif
         rhs_mat = rhs_mat + current_vel*omeg0_v_mat(:,:,1)
      endif
!
      end subroutine mbe_set_rhsmat1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
      subroutine mbe_spline_prep(y_in,spl)
!
!  Prepare a natural spline. The arrays y_in and spl are expected to have 
!  the dimension 0:n_time_steps (n_time_steps is a global variable set by
!  mbe_propagate_2, as is the array t_mesh).
!
!  Based on the subroutine SPLINE of W H Press, B P Flannery, 
!  S A Teukolsky and W T Vetterling, "Numerical Recipes" (CUP,
!  Cambridge, 1986).
!
      real(kd), dimension(0:), intent(in) :: y_in
      real(kd), dimension(0:), intent(out) :: spl
      real(kd), dimension(:), allocatable :: upd
      real(kd) :: denom,ratio1,ratio2,ratio2pr
      integer :: i,istat
!
      allocate(upd(0:n_time_steps),stat=istat)
         if(istat.gt.0)then
         print*,'mbe_spline_prep: error at the allocation of the array upd.'
         call mbe_error
      endif
!
      spl(0) = 0.0_kd
      upd(0) = 0.0_kd
      ratio2pr = (y_in(1) - y_in(0))/(t_mesh(1) - t_mesh(0))
      do i = 1,n_time_steps-1
         ratio1 = (t_mesh(i) - t_mesh(i-1))/(t_mesh(i+1) - t_mesh(i-1))
         denom = ratio1*spl(i-1)+2.0_kd
         spl(i) = (ratio1-1.0_kd)/denom
         ratio2 = (y_in(i+1) - y_in(i))/(t_mesh(i+1) - t_mesh(i))
         upd(i) = 6.0_kd*(ratio2 - ratio2pr)/(t_mesh(i+1) - t_mesh(i-1)) &
                  - ratio1*upd(i-1)
         upd(i) = upd(i)/denom
         ratio2pr = ratio2
      enddo
      spl(n_time_steps) = 0.0_kd
      do i = n_time_steps-1,0,-1
         spl(i) = spl(i)*spl(i+1)+upd(i)
      enddo
!
      deallocate(upd)
!     
      end subroutine mbe_spline_prep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine mbe_spline_int(t,res,y_in,spl)
!
!  Spline interpolation of the arrays y_in, using the array spl previously
!  prepared by mbe_spline_prep.
!
!  The global variable ktime is set by mbe_propagate_2 or mbe_tdint_2.
!
!  Based on the subroutine SPLINT of W H Press, B P Flannery, 
!  S A Teukolsky and W T Vetterling, "Numerical Recipes" (CUP,
!  Cambridge, 1986).
!
      real(kd), intent(in) :: t
      real(kd), dimension(0:), intent(in) :: y_in,spl
      real(kd), intent(out) :: res
!     real(kd) :: a,b,h
      real(kd) :: d1,d2,f,h
      integer :: kp1
!
      kp1 = ktime + 1
!
!  Check the values t and of ktime:
      if(ktime.lt.0 .or. kp1.gt.n_time_steps)then
         print*,'mbe_spline_int: Bad value of ktime. ',ktime
         call mbe_error
      elseif(t.lt.t_mesh(ktime) .or. t.gt.t_mesh(kp1))then
         print*,'mbe_spline_int: t and ktime are inconsistent. '
         print*,ktime,t_mesh(ktime),t,t_mesh(kp1)
         call mbe_error
      endif
!
      h = t_mesh(kp1) - t_mesh(ktime)
      f = h*h/6.0_kd
      d1 = (t_mesh(kp1) - t)/h
      d2 = (t - t_mesh(ktime))/h
!     res = a*y(ktime) + b*y(kp1) + &
!              ((a**3-a)*y2(ktime)+(b**3-b)*y2(kp1))*(h**2)/6.0_kd
      res = d1 * (y_in(ktime) + ((d1*d1-1.0_kd)*spl(ktime))*f ) +  &
            d2 * (y_in(kp1)   + ((d2*d2-1.0_kd)*spl(kp1)  )*f ) 
!
      end subroutine mbe_spline_int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
      subroutine mbe_calc_flds(t)
!
!  mbe_calc_flds calculates the amplitude of the applied field(s) at
!  time t and puts the results in the global variables Field_p_re,
!  Field_p_im, Field_c_re and Field_c_im. It is assumed that the
!  arrays Field_p_re_vec, Fd_pre2, Field_p_im_vec, Fd_pim2,
!  Field_c_re_vec, Fd_cre2, Field_c_im_vec and Fd_cim2 are ready to be
!  used.
!
      implicit none
      real(kd), intent(in) :: t
      double complex :: cfieldt
      double precision :: tdble
!
      if(fl_spline)then
         if(fl_twofields)then
            call mbe_spline_int(t,Field_p_re,Field_p_re_vec(:),Fd_pre2)
            call mbe_spline_int(t,Field_p_im,Field_p_im_vec(:),Fd_pim2)
            call mbe_spline_int(t,Field_c_re,Field_c_re_vec(:),Fd_cre2)
            call mbe_spline_int(t,Field_c_im,Field_c_im_vec(:),Fd_cim2)
         else
            call mbe_spline_int(t,Field_p_re,Field_p_re_vec(:),Fd_pre2)
            call mbe_spline_int(t,Field_p_im,Field_p_im_vec(:),Fd_pim2)
         endif
      else
         tdble = dble(t)
         if(fl_twofields)then
            call mbe_calc_applfld(tdble,tw_p,pulse_type_p,t0_p,t1_p,   &
                                  cfield0_p,cfieldt)
            Field_p_re = real(dreal(cfieldt),kd)
            Field_p_im = real(dimag(cfieldt),kd)
            call mbe_calc_applfld(tdble,tw_c,pulse_type_c,t0_c,t1_c,   &
                                  cfield0_c,cfieldt)
            Field_c_re = real(dreal(cfieldt),kd)
            Field_c_im = real(dimag(cfieldt),kd)
         else
            call mbe_calc_applfld(tdble,tw_p,pulse_type_p,t0_p,t1_p,   &
                                  cfield0_p,cfieldt)
            Field_p_re = real(dreal(cfieldt),kd)
            Field_p_im = real(dimag(cfieldt),kd)
         endif
      endif
!
      end subroutine mbe_calc_flds
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
      subroutine mbe_intobe(imethod,t1,t2,nsteps,rho_vec)
!
!  Integrate the Optical Bloch equations over time from t1 to t2,
!  starting from the initial values passed to the subprogram through
!  the vector rho_vec, and return the result to the calling program
!  through the same vector (which is thus overwritten by obe_intobe).
!  This vector contains the populations and the real and imaginary parts
!  of the (i,j) coherences for j > i.
!
!  Currently, the interval [t1,t2] is divided into nsteps steps of duration
!  tstep, and the solution is propagated over each step using either the
!  4th order Runge-Kutta method (when imethod.eq.1) or Butcher's fifth order
!  Runge-Kutta method (imethod.eq.5) or by dop853 (imethod.eq.4).
!
      implicit none
      integer, parameter :: neqs = nst*nst
      real(kd), dimension (neqs), intent(inout) :: rho_vec
      real(kd), dimension (neqs,neqs) :: rhs_mat
      real(kd), intent(in) :: t1,t2
      integer, intent(in) :: imethod,nsteps
      real(kd) :: t,tstep
      integer :: n,j
      logical :: fl_recalc
!
!  Calculate tstep...
      tstep = (t2-t1)/real(nsteps,kd)
!
!  ... and integrate, starting with the populations and coherences passed
!  to mbe_intobe through rho_vec. This initial density matrix is
!  propagated by mbe_rk4 or rk5, the Runge-Kutta solvers, or by dop853.
!  At each step, rho_vec is overwritten with the propagated solution.
      t = t1
      do n = 1,nsteps
!  fl_recalc is used by mbe_rhs, which is called by mbe_rk4 and mbe_rk5.
         if(n.eq.1)then
            fl_recalc = .true.
         else
            fl_recalc = .false.
         endif
         if(imethod.eq.1)then
            call mbe_rk4(rho_vec,t,tstep,rho_vec)
         elseif(imethod.eq.2)then
            print*,'mbe_intobe: Option 2 for imethod is currently disabled.'
            call mbe_error
         elseif(imethod.eq.4)then
            call mbe_dop853(rho_vec,t,tstep,rho_vec)
         elseif(imethod.eq.5)then
            call mbe_rk5(rho_vec,t,tstep,rho_vec)
         else
            print*,'mbe_intobe: Illegal value of imethod ',imethod
            call mbe_error
         endif
         t = t1 + real(n,kd)*tstep
      enddo
!
      contains
 
!!!!!!!!!!!!!!!!!!!
 
      subroutine mbe_rhs(t,rho_vec,rhodot_vec)
!
!  mbe_rhs calculates the right-hand side of the master equation and returns
!  the result in rhodot_vec.
!
!  The calculation assumes that the matrix representing the righ-hand side
!  of the system of equations does not change between calls to mbe_rhs for
!  the same value of t. I.e., this matrix is not recalculated if t has the
!  same value as at the previous call. The matrix is also recalculated if
!  the logical variable fl_recalc .eqv. .true. (this variable is set in
!  mbe_intobe).
!
      implicit none
      integer, parameter :: neqs = nst*nst
      real(kd), dimension (neqs), intent(in) :: rho_vec
      real(kd), dimension (neqs), intent(out) :: rhodot_vec
      real(kd), intent(in) :: t
      real(kd), dimension (neqs,neqs), save :: rhs_mat
      real(kd), save :: told = huge(1.0_kd)
      real(kd) :: t_used
      double complex :: cfieldt
      double precision :: tdble
!
!  Obtain the matrix used to calculate the right-hand sides for the 
!  current values of the complex field amplitudes, if necessary.
!  The time assumed in the calculation is normally the value of t,
!  unless this value if larger than the upper end of the integration
!  range (t2), in which case it is taken to be t2 (overshooting t2 is
!  normally due to truncation errors).
      if(fl_recalc .or. t.ne.told)then
         if(t.gt.t2)then
            t_used = t2
         else
            t_used = t
         endif
         if(fl_spline)then
            if(fl_twofields)then
               call mbe_spline_int(t_used,Field_p_re,Field_p_re_vec(:),Fd_pre2)
               call mbe_spline_int(t_used,Field_p_im,Field_p_im_vec(:),Fd_pim2)
               call mbe_spline_int(t_used,Field_c_re,Field_c_re_vec(:),Fd_cre2)
               call mbe_spline_int(t_used,Field_c_im,Field_c_im_vec(:),Fd_cim2)
               call mbe_set_rhsmat2(rhs_mat,1)
            else
               call mbe_spline_int(t_used,Field_p_re,Field_p_re_vec(:),Fd_pre2)
               call mbe_spline_int(t_used,Field_p_im,Field_p_im_vec(:),Fd_pim2)
               call mbe_set_rhsmat1(rhs_mat,1)
            endif
         else
            tdble = dble(t_used)
            if(fl_twofields)then
               call mbe_calc_applfld(tdble,tw_p,pulse_type_p,t0_p,t1_p,   &
                                     cfield0_p,cfieldt)
               Field_p_re = real(dreal(cfieldt),kd)
               Field_p_im = real(dimag(cfieldt),kd)
               call mbe_calc_applfld(tdble,tw_c,pulse_type_c,t0_c,t1_c,   &
                                     cfield0_c,cfieldt)
               Field_c_re = real(dreal(cfieldt),kd)
               Field_c_im = real(dimag(cfieldt),kd)
               call mbe_set_rhsmat2(rhs_mat,1)
            else
               call mbe_calc_applfld(tdble,tw_p,pulse_type_p,t0_p,t1_p,   &
                                     cfield0_p,cfieldt)
               Field_p_re = real(dreal(cfieldt),kd)
               Field_p_im = real(dimag(cfieldt),kd)
               call mbe_set_rhsmat1(rhs_mat,1)
            endif
         endif
!
         fl_recalc = .false.
         told = t_used
      endif
!
      rhodot_vec = matmul(rhs_mat,rho_vec)
!
      end subroutine mbe_rhs
 
!!!!!!!!!!!!!!!!!!!
 
      subroutine mbe_rk4(rho,tstart,tstep,rho_out)
!
!  Propagates the array rho over one time step of duration tstep using
!  the fourth order Runge-Kutta method. The derivatives at time t are
!  calculated in the subroutine mbe_rhs(t,rho,rhodot).
!
!  rho: At entry the array rho at t = tstart. At return, the array rho
!       at t = tstart + tstep
!  tstart: The initial time. No change at return in this version of this
!          subroutine.
!  rho_out: At return, the array rho at t = tstart + tstep. (This array
!           may be the same as rho.)
!
      implicit none
      integer, parameter :: neqs = nst*nst
      integer, parameter :: n = neqs
      real(kd), dimension(n), intent(in) :: rho
      real(kd), dimension(n), intent(out) :: rho_out
      real(kd), intent(in) :: tstart
      real(kd), intent(in) :: tstep
      real(kd) :: half_tstep,thalfstep
      real(kd), dimension(n) :: rhodot,rhok,rhokdot1,rhokdot2
      integer :: istat
!
!  Integrate the equations. Start by calculating the derivatives at the
!  initial value of t:
      call mbe_rhs(tstart,rho,rhodot)
!
      half_tstep = 0.5_kd*tstep
      thalfstep = tstart + half_tstep
!
      rhok = rho + half_tstep*rhodot           ! rhok is now rho + k_1 / 2
      call mbe_rhs(thalfstep,rhok,rhokdot1)
      rhok = rho + half_tstep*rhokdot1         ! rhok is now rho + k_2 / 2
      call mbe_rhs(thalfstep,rhok,rhokdot2)
      rhodot = rhodot + 2.0_kd*(rhokdot1+rhokdot2)
      rhok = rho + tstep*rhokdot2              ! rhok is now rho + k_3
      call mbe_rhs(tstart+tstep,rhok,rhokdot1)
      rhodot = rhodot + rhokdot1
!
      rho_out = rho + (tstep/6.0_kd)*rhodot
!
      end subroutine mbe_rk4

!!!!!!!!!!!!!!!!!!!
 
      subroutine mbe_rk5(y,x,h,yout)
!
!  Fifth order Runge Kutta ("Butcher's formula"). See mbe_rk4 for further information.
!
      implicit none
      integer, parameter :: neqs = nst*nst
      integer, parameter :: n = neqs
      real(kd), dimension(n), intent(in) :: y
      real(kd), dimension(n), intent(out) :: yout
      real(kd), intent(in) :: x,h
      real(kd), dimension(n) :: yst,dy1,dy2,dy3,dy4,dy5,dy6
      real(kd) :: hh
      integer :: i
!
      call mbe_rhs(x,y,dy1)
      yst = y + h*dy1/4.0_kd
      call mbe_rhs(x+h/4.0_kd,yst,dy2)
      yst = y + h*(dy1+dy2)/8.0_kd
      call mbe_rhs(x+h/4.0_kd,yst,dy3)
      yst = y + h*(-dy2/2.0_kd+dy3)
      call mbe_rhs(x+h/2.0_kd,yst,dy4)
      yst = y + h*(3.0_kd*dy1+9.0_kd*dy4)/16.0_kd
      call mbe_rhs(x+3.0_kd*h/4.0_kd,yst,dy5)
      yst = y + h*(-3.0_kd*dy1+2.0_kd*dy2+12.0_kd*dy3-12.0_kd*dy4+ &
                                                8.0_kd*dy5)/7.0_kd
      call mbe_rhs(x+h,yst,dy6)
      yout = y + h*(7.0_kd*dy1+32.0_kd*dy3+12.0_kd*dy4+32.0_kd*dy5+ &
                                                7.0_kd*dy6)/90.0_kd
!
      end subroutine mbe_rk5
 
!!!!!!!!!!!!!!!!!!!

      subroutine mbe_dop853(rho,tstart,tstep,rho_out)
!
!  Propagates the array rho over one time step of duration tstep using
!  the dop853 subroutines, which mbe_dop853 is merely an interface
!  with.
!
!  The arguments of mbe_dop853 are the same as those of mbe_rk4. See
!  mbe_rk4 for further information.
!
      implicit none
      integer, parameter :: neqs = nst*nst
      integer, parameter :: n = neqs
      integer, parameter :: nrdens = 0  ! No dense output required.
      integer, parameter :: lwork = 11*n + 8*nrdens + 21
      integer, parameter :: liwork = nrdens + 21
      real(kd), dimension(n), intent(in) :: rho
      real(kd), dimension(n), intent(out) :: rho_out
      real(kd), intent(in) :: tstart
      real(kd), intent(in) :: tstep
      double precision, dimension(lwork) :: work
      integer, dimension(liwork) :: iwork
      double precision, dimension(1) :: rtol,atol,rpar
      double precision :: uround
      integer, dimension(1) :: ipar
      integer :: itol,iout,idid,istat,llwork,lliwork
      external fcn_dummy,solout_dummy
      save
!
!  Check that the module is programmed in double precision, to avoid
!  a mismatch of variable type with dop853
      if(kd.ne.kind(1.d0))then
         print*,'mbe_dop853 can be used only if the rest of the module'
         print*,'is programmed in double precision.'
         call mbe_error
      endif
!
!  fl_set_tol_dop853, rtol_top853 and atol_top853 are global parameters.
      if(.not.fl_set_tol_dop853)then
         print*,'mbe_dop853: the tolerance parameters rtol and atol'
         print*,'should have been initialized by a call to'
         print*,'obe_set_tol_dop853 before the first call to'
         print*,'mbe_dop853.'
         call mbe_error
      endif
      rtol = rtol_dop853
      atol = atol_dop853
      itol = 0
      iout = 0
      uround = spacing(1.d0)
!
      work = 0.d0
      work(1) = uround
      iwork = 0
      llwork = lwork
      lliwork = liwork
!  The following arguments are not used by dop853 in the present
!  context, but need to be included in the calling list of consistency:
!      fcn_dummy, solout_dummy, rpar, ipar.
!  fcn_dummy and solout_dummy are the names of dummy external subroutines
!  which do not do anything. These two subroutines do not belong to
!  this module but have been appended at the end of this file.
!
      rho_out = rho
      call dop853(neqs,fcn_dummy,t,rho_out,t+tstep,          &
         rtol,atol,itol,solout_dummy,iout,work,llwork,       &
         iwork,lliwork,rpar,ipar,idid)
!
      if(idid.lt.0)then
         print*,'mbe_dop853: dop853 returns with an error condition.'
         print*,idid
         call mbe_error
      endif
!
      end subroutine mbe_dop853

!!!!!!!!!!!!!!!!!!!

! *********************************************************************
!  The following is a Fortran 90 translation of Hairer and Wanner's
!  subroutine dop853. There are no substantive changes from the
!  original Fortran 77 code, apart from the disabling of the calls to
!  the external subroutine SOLOUT and the issuing of an error message
!  and error condition if DOP853 is called with IOUT > 0. These changes
!  can be easily reversed if this would become desirable in a later
!  stage of development of the module.
!
!  Original Fortran 77 code retrieved on the 15th of September 2022
!  from http://www.unige.ch/~hairer/prog/nonstiff/dop853.f.
!
!  The original code is subject to the following copyright notice:
!
!  Copyright (c) 2004, UNIGE
!
!  Redistribution and use in source and binary forms, with or without 
!  modification, are permitted provided that the following conditions are 
!  met:
!  
!  - Redistributions of source code must retain the above copyright 
!  notice, this list of conditions and the following disclaimer.
!  
!  - Redistributions in binary form must reproduce the above copyright 
!  notice, this list of conditions and the following disclaimer in the 
!  documentation and/or other materials provided with the distribution.
!  
!  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS 
!  IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED 
!  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
!  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR 
!  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
!  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
!  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
!  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
!  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
!  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
!  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! *********************************************************************
!
      SUBROUTINE DOP853(N,FCN,X,Y,XEND,                          &
     &                  RTOL,ATOL,ITOL,                          &
     &                  SOLOUT,IOUT,                             &
     &                  WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)
! ----------------------------------------------------------
!     NUMERICAL SOLUTION OF A SYSTEM OF FIRST 0RDER
!     ORDINARY DIFFERENTIAL EQUATIONS  Y'=F(X,Y).
!     THIS IS AN EXPLICIT RUNGE-KUTTA METHOD OF ORDER 8(5,3)  
!     DUE TO DORMAND & PRINCE (WITH STEPSIZE CONTROL AND
!     DENSE OUTPUT)
!
!     AUTHORS: E. HAIRER AND G. WANNER
!              UNIVERSITE DE GENEVE, DEPT. DE MATHEMATIQUES
!              CH-1211 GENEVE 24, SWITZERLAND 
!              E-MAIL:  Ernst.Hairer@unige.ch
!                       Gerhard.Wanner@unige.ch
!     
!     THIS CODE IS DESCRIBED IN:
!         E. HAIRER, S.P. NORSETT AND G. WANNER, SOLVING ORDINARY
!         DIFFERENTIAL EQUATIONS I. NONSTIFF PROBLEMS. 2ND EDITION. 
!         SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS, 
!         SPRINGER-VERLAG (1993)
!      
!     VERSION OF OCTOBER 11, 2009
!      (new option IOUT=3 for sparse dense output)
!
!     INPUT PARAMETERS  
!     ----------------  
!     N           DIMENSION OF THE SYSTEM 
!
!     FCN         NAME (EXTERNAL) OF SUBROUTINE COMPUTING THE
!                 VALUE OF F(X,Y):
!                    SUBROUTINE FCN(N,X,Y,F,RPAR,IPAR)
!                    DOUBLE PRECISION X,Y(N),F(N)
!                    F(1)=...   ETC.
!
!     X           INITIAL X-VALUE
!
!     Y(N)        INITIAL VALUES FOR Y
!
!     XEND        FINAL X-VALUE (XEND-X MAY BE POSITIVE OR NEGATIVE)
!
!     RTOL,ATOL   RELATIVE AND ABSOLUTE ERROR TOLERANCES. THEY
!                 CAN BE BOTH SCALARS OR ELSE BOTH VECTORS OF LENGTH N.
!                 ATOL SHOULD BE STRICTLY POSITIVE (POSSIBLY VERY SMALL)
!
!     ITOL        SWITCH FOR RTOL AND ATOL:
!                   ITOL=0: BOTH RTOL AND ATOL ARE SCALARS.
!                     THE CODE KEEPS, ROUGHLY, THE LOCAL ERROR OF
!                     Y(I) BELOW RTOL*ABS(Y(I))+ATOL
!                   ITOL=1: BOTH RTOL AND ATOL ARE VECTORS.
!                     THE CODE KEEPS THE LOCAL ERROR OF Y(I) BELOW
!                     RTOL(I)*ABS(Y(I))+ATOL(I).
!
!     SOLOUT      NAME (EXTERNAL) OF SUBROUTINE PROVIDING THE
!                 NUMERICAL SOLUTION DURING INTEGRATION. 
!                 IF IOUT.GE.1, IT IS CALLED DURING INTEGRATION.
!                 SUPPLY A DUMMY SUBROUTINE IF IOUT=0. 
!                 IT MUST HAVE THE FORM
!                    SUBROUTINE SOLOUT (NR,XOLD,X,Y,N,CON,ICOMP,ND,
!                                       RPAR,IPAR,IRTRN,XOUT)
!                    DIMENSION Y(N),CON(8*ND),ICOMP(ND)
!                    ....  
!                 SOLOUT FURNISHES THE SOLUTION "Y" AT THE NR-TH
!                    GRID-POINT "X" (THEREBY THE INITIAL VALUE IS
!                    THE FIRST GRID-POINT).
!                 "XOLD" IS THE PRECEEDING GRID-POINT.
!                 "IRTRN" SERVES TO INTERRUPT THE INTEGRATION. IF IRTRN
!                    IS SET <0, DOP853 WILL RETURN TO THE CALLING
!                    PROGRAM.
!                    IF THE NUMERICAL SOLUTION IS ALTERED IN SOLOUT,
!                    SET  IRTRN = 2
!                 "XOUT" CAN BE USED FOR EFFICIENT INTERMEDIATE OUTPUT
!                    IF ONE PUTS IOUT=3
!                    WHEN NR=1 DEFINE THE FIRST OUTPUT POINT XOUT IN
!                    SOLOUT.
!                      THE SUBROUTINE SOLOUT WILL BE CALLED ONLY WHEN 
!                      XOUT IS IN THE INTERVAL [XOLD,X]; DURING THIS
!                      CALL
!                      A NEW VALUE FOR XOUT CAN BE DEFINED, ETC.
!           
!          -----  CONTINUOUS OUTPUT: -----
!                 DURING CALLS TO "SOLOUT", A CONTINUOUS SOLUTION
!                 FOR THE INTERVAL [XOLD,X] IS AVAILABLE THROUGH
!                 THE FUNCTION
!                        >>>   CONTD8(I,S,CON,ICOMP,ND)   <<<
!                 WHICH PROVIDES AN APPROXIMATION TO THE I-TH
!                 COMPONENT OF THE SOLUTION AT THE POINT S. THE VALUE
!                 S SHOULD LIE IN THE INTERVAL [XOLD,X].
!           
!     IOUT        SWITCH FOR CALLING THE SUBROUTINE SOLOUT:
!                    IOUT=0: SUBROUTINE IS NEVER CALLED
!                    IOUT=1: SUBROUTINE IS CALLED AFTER EVERY SUCCESSFUL
!                    STEP
!                    IOUT=2: DENSE OUTPUT IS PERFORMED AFTER EVERY
!                    SUCCESSFUL STEP
!                            (IN THIS CASE IWORK(5) MUST BE SPECIFIED)
!                    IOUT=3: DENSE OUTPUT IS PERFORMED IN STEPS DEFINED
!                    BY THE USER
!                            (SEE "XOUT" ABOVE)
!
!     WORK        ARRAY OF WORKING SPACE OF LENGTH "LWORK".
!                 WORK(1),...,WORK(20) SERVE AS PARAMETERS FOR THE CODE.
!                 FOR STANDARD USE, SET THEM TO ZERO BEFORE CALLING.
!                 "LWORK" MUST BE AT LEAST  11*N+8*NRDENS+21
!                 WHERE  NRDENS = IWORK(5)
!
!     LWORK       DECLARED LENGHT OF ARRAY "WORK".
!
!     IWORK       INTEGER WORKING SPACE OF LENGHT "LIWORK".
!                 IWORK(1),...,IWORK(20) SERVE AS PARAMETERS FOR THE
!                 CODE.
!                 FOR STANDARD USE, SET THEM TO ZERO BEFORE CALLING.
!                 "LIWORK" MUST BE AT LEAST NRDENS+21 .
!
!     LIWORK      DECLARED LENGHT OF ARRAY "IWORK".
!
!     RPAR, IPAR  REAL AND INTEGER PARAMETERS (OR PARAMETER ARRAYS)
!     WHICH  
!                 CAN BE USED FOR COMMUNICATION BETWEEN YOUR CALLING
!                 PROGRAM AND THE FCN, JAC, MAS, SOLOUT SUBROUTINES. 
!
!-----------------------------------------------------------------------
! 
!     SOPHISTICATED SETTING OF PARAMETERS
!     -----------------------------------
!              SEVERAL PARAMETERS (WORK(1),...,IWORK(1),...) ALLOW
!              TO ADAPT THE CODE TO THE PROBLEM AND TO THE NEEDS OF
!              THE USER. FOR ZERO INPUT, THE CODE CHOOSES DEFAULT
!              VALUES.
!
!    WORK(1)   UROUND, THE ROUNDING UNIT, DEFAULT 2.3D-16.
!
!    WORK(2)   THE SAFETY FACTOR IN STEP SIZE PREDICTION,
!              DEFAULT 0.9D0.
!
!    WORK(3), WORK(4)   PARAMETERS FOR STEP SIZE SELECTION
!              THE NEW STEP SIZE IS CHOSEN SUBJECT TO THE RESTRICTION
!                 WORK(3) <= HNEW/HOLD <= WORK(4)
!              DEFAULT VALUES: WORK(3)=0.333D0, WORK(4)=6.D0
!
!    WORK(5)   IS THE "BETA" FOR STABILIZED STEP SIZE CONTROL
!              (SEE SECTION IV.2). POSITIVE VALUES OF BETA ( <= 0.04 )
!              MAKE THE STEP SIZE CONTROL MORE STABLE.
!              NEGATIVE WORK(5) PROVOKE BETA=0.
!              DEFAULT 0.0D0.
!
!    WORK(6)   MAXIMAL STEP SIZE, DEFAULT XEND-X.
!
!    WORK(7)   INITIAL STEP SIZE, FOR WORK(7)=0.D0 AN INITIAL GUESS
!              IS COMPUTED WITH HELP OF THE FUNCTION HINIT
!
!    IWORK(1)  THIS IS THE MAXIMAL NUMBER OF ALLOWED STEPS.
!              THE DEFAULT VALUE (FOR IWORK(1)=0) IS 100000.
!
!    IWORK(2)  SWITCH FOR THE CHOICE OF THE COEFFICIENTS
!              IF IWORK(2).EQ.1  METHOD DOP853 OF DORMAND AND PRINCE
!              (SECTION II.6).
!              THE DEFAULT VALUE (FOR IWORK(2)=0) IS IWORK(2)=1.
!
!    IWORK(3)  SWITCH FOR PRINTING ERROR MESSAGES
!              IF IWORK(3).LT.0 NO MESSAGES ARE BEING PRINTED
!              IF IWORK(3).GT.0 MESSAGES ARE PRINTED WITH
!              WRITE (IWORK(3),*) ...  
!              DEFAULT VALUE (FOR IWORK(3)=0) IS IWORK(3)=6
!
!    IWORK(4)  TEST FOR STIFFNESS IS ACTIVATED AFTER STEP NUMBER
!              J*IWORK(4) (J INTEGER), PROVIDED IWORK(4).GT.0.
!              FOR NEGATIVE IWORK(4) THE STIFFNESS TEST IS
!              NEVER ACTIVATED; DEFAULT VALUE IS IWORK(4)=1000
!
!    IWORK(5)  = NRDENS = NUMBER OF COMPONENTS, FOR WHICH DENSE OUTPUT
!              IS REQUIRED; DEFAULT VALUE IS IWORK(5)=0;
!              FOR   0 < NRDENS < N   THE COMPONENTS (FOR WHICH DENSE
!              OUTPUT IS REQUIRED) HAVE TO BE SPECIFIED IN
!              IWORK(21),...,IWORK(NRDENS+20);
!              FOR  NRDENS=N  THIS IS DONE BY THE CODE.
!
!----------------------------------------------------------------------
!
!     OUTPUT PARAMETERS 
!     ----------------- 
!     X           X-VALUE FOR WHICH THE SOLUTION HAS BEEN COMPUTED
!                 (AFTER SUCCESSFUL RETURN X=XEND).
!
!     Y(N)        NUMERICAL SOLUTION AT X
! 
!     H           PREDICTED STEP SIZE OF THE LAST ACCEPTED STEP
!
!     IDID        REPORTS ON SUCCESSFULNESS UPON RETURN:
!                   IDID= 1  COMPUTATION SUCCESSFUL,
!                   IDID= 2  COMPUT. SUCCESSFUL (INTERRUPTED BY SOLOUT)
!                   IDID=-1  INPUT IS NOT CONSISTENT,
!                   IDID=-2  LARGER NMAX IS NEEDED,
!                   IDID=-3  STEP SIZE BECOMES TOO SMALL.
!                   IDID=-4  PROBLEM IS PROBABLY STIFF (INTERRUPTED).
!
!   IWORK(17)  NFCN    NUMBER OF FUNCTION EVALUATIONS
!   IWORK(18)  NSTEP   NUMBER OF COMPUTED STEPS
!   IWORK(19)  NACCPT  NUMBER OF ACCEPTED STEPS
!   IWORK(20)  NREJCT  NUMBER OF REJECTED STEPS (DUE TO ERROR TEST),
!                      (STEP REJECTIONS IN THE FIRST STEP ARE NOT
!                      COUNTED)
!-----------------------------------------------------------------------
! *** *** *** *** *** *** *** *** *** *** *** *** ***
!          DECLARATIONS 
! *** *** *** *** *** *** *** *** *** *** *** *** ***
!!!   IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!!!   IMPLICIT INTEGER (I-N)
      implicit none
      DIMENSION Y(N),ATOL(*),RTOL(*),WORK(LWORK),IWORK(LIWORK)
      DIMENSION RPAR(*),IPAR(*)
      double precision :: x,y,xend,rtol,atol,work,rpar
      integer :: n,itol,iout,lwork,iwork,liwork,ipar,idid
      double precision :: uround,safe,fac1,fac2,beta,hmax,h
      integer :: i,nfcn,nstep,naccpt,nrejct,iprint,nmax,meth,  &
                 nstiff,nrdens,iek1,iek2,iek3,iek4,iek5,iek6,  &
                 iek7,iek8,iek9,iek10,iey1,ieco,istore,icomp
      LOGICAL ARRET
      EXTERNAL FCN,SOLOUT
! *** *** *** *** *** *** ***
!        SETTING THE PARAMETERS 
! *** *** *** *** *** *** ***
      NFCN=0
      NSTEP=0
      NACCPT=0
      NREJCT=0
      ARRET=.FALSE.
! -------- IPRINT FOR MONITORING THE PRINTING
      IF(IWORK(3).EQ.0)THEN
         IPRINT=6
      ELSE
         IPRINT=IWORK(3)
      END IF
! -------- NMAX , THE MAXIMAL NUMBER OF STEPS ----- 
      IF(IWORK(1).EQ.0)THEN
         NMAX=100000
      ELSE
         NMAX=IWORK(1)
         IF(NMAX.LE.0)THEN
            IF (IPRINT.GT.0) WRITE(IPRINT,*)        &
     &          ' WRONG INPUT IWORK(1)=',IWORK(1)
            ARRET=.TRUE.
         END IF
      END IF
! -------- METH   COEFFICIENTS OF THE METHOD
      IF(IWORK(2).EQ.0)THEN
         METH=1
      ELSE
         METH=IWORK(2)
         IF(METH.LE.0.OR.METH.GE.4)THEN
            IF (IPRINT.GT.0) WRITE(IPRINT,*)        &
     &          ' CURIOUS INPUT IWORK(2)=',IWORK(2)
            ARRET=.TRUE.
         END IF
      END IF  
! -------- NSTIFF   PARAMETER FOR STIFFNESS DETECTION  
      NSTIFF=IWORK(4) 
      IF (NSTIFF.EQ.0) NSTIFF=1000
      IF (NSTIFF.LT.0) NSTIFF=NMAX+10
! -------- NRDENS   NUMBER OF DENSE OUTPUT COMPONENTS
      NRDENS=IWORK(5)
      IF(NRDENS.LT.0.OR.NRDENS.GT.N)THEN
         IF (IPRINT.GT.0) WRITE(IPRINT,*)        &
     &           ' CURIOUS INPUT IWORK(5)=',IWORK(5)
         ARRET=.TRUE.
      ELSE
         IF(NRDENS.GT.0.AND.IOUT.LT.2)THEN
            IF (IPRINT.GT.0) WRITE(IPRINT,*)        &
     &       ' WARNING: PUT IOUT=2 OR IOUT=3 FOR DENSE OUTPUT '
         END IF 
         IF (NRDENS.EQ.N) THEN
            DO I=1,NRDENS
               IWORK(I+20)=I
            END DO
         END IF
      END IF       
! -------- UROUND   SMALLEST NUMBER SATISFYING 1.D0+UROUND>1.D0  
      IF(WORK(1).EQ.0.D0)THEN
         UROUND=2.3D-16
      ELSE
         UROUND=WORK(1)
         IF(UROUND.LE.1.D-35.OR.UROUND.GE.1.D0)THEN
            IF (IPRINT.GT.0) WRITE(IPRINT,*)        &
     &        ' WHICH MACHINE DO YOU HAVE? YOUR UROUND WAS:',WORK(1)
            ARRET=.TRUE.
         END IF
      END IF
! -------  SAFETY FACTOR -------------
      IF(WORK(2).EQ.0.D0)THEN
         SAFE=0.9D0
      ELSE
         SAFE=WORK(2)
         IF(SAFE.GE.1.D0.OR.SAFE.LE.1.D-4)THEN
            IF (IPRINT.GT.0) WRITE(IPRINT,*)        &
     &          ' CURIOUS INPUT FOR SAFETY FACTOR WORK(2)=',WORK(2)
            ARRET=.TRUE.
         END IF
      END IF
! -------  FAC1,FAC2     PARAMETERS FOR STEP SIZE SELECTION
      IF(WORK(3).EQ.0.D0)THEN
         FAC1=0.333D0
      ELSE
         FAC1=WORK(3)
      END IF
      IF(WORK(4).EQ.0.D0)THEN
         FAC2=6.D0
      ELSE
         FAC2=WORK(4)
      END IF
! --------- BETA FOR STEP CONTROL STABILIZATION -----------
      IF(WORK(5).EQ.0.D0)THEN
         BETA=0.0D0
      ELSE
         IF(WORK(5).LT.0.D0)THEN
            BETA=0.D0
         ELSE
            BETA=WORK(5)
            IF(BETA.GT.0.2D0)THEN
               IF (IPRINT.GT.0) WRITE(IPRINT,*)        &
     &          ' CURIOUS INPUT FOR BETA: WORK(5)=',WORK(5)
            ARRET=.TRUE.
         END IF
         END IF
      END IF
! -------- MAXIMAL STEP SIZE
      IF(WORK(6).EQ.0.D0)THEN
         HMAX=XEND-X
      ELSE
         HMAX=WORK(6)
      END IF
! -------- INITIAL STEP SIZE
      H=WORK(7)
! ------- PREPARE THE ENTRY-POINTS FOR THE ARRAYS IN WORK -----
      IEK1=21
      IEK2=IEK1+N
      IEK3=IEK2+N
      IEK4=IEK3+N
      IEK5=IEK4+N
      IEK6=IEK5+N
      IEK7=IEK6+N
      IEK8=IEK7+N
      IEK9=IEK8+N
      IEK10=IEK9+N
      IEY1=IEK10+N
      IECO=IEY1+N
! ------ TOTAL STORAGE REQUIREMENT -----------
      ISTORE=IECO+8*NRDENS-1
      IF(ISTORE.GT.LWORK)THEN
        IF (IPRINT.GT.0) WRITE(IPRINT,*)        &
     &   ' INSUFFICIENT STORAGE FOR WORK, MIN. LWORK=',ISTORE
        ARRET=.TRUE.
      END IF
      ICOMP=21
      ISTORE=ICOMP+NRDENS-1
      IF(ISTORE.GT.LIWORK)THEN
        IF (IPRINT.GT.0) WRITE(IPRINT,*)        &
     &   ' INSUFFICIENT STORAGE FOR IWORK, MIN. LIWORK=',ISTORE
        ARRET=.TRUE.
      END IF
! -------- WHEN A FAIL HAS OCCURED, WE RETURN WITH IDID=-1
      IF (ARRET) THEN
         IDID=-1
         RETURN
      END IF
! -------- CALL TO CORE INTEGRATOR ------------
      CALL DP86CO(N,FCN,X,Y,XEND,HMAX,H,RTOL,ATOL,ITOL,IPRINT,         &
     &   SOLOUT,IOUT,IDID,NMAX,UROUND,METH,NSTIFF,SAFE,BETA,FAC1,FAC2, &
     &   WORK(IEK1),WORK(IEK2),WORK(IEK3),WORK(IEK4),WORK(IEK5),       &
     &   WORK(IEK6),WORK(IEK7),WORK(IEK8),WORK(IEK9),WORK(IEK10),      &
     &   WORK(IEY1),WORK(IECO),IWORK(ICOMP),NRDENS,RPAR,IPAR,          &
     &   NFCN,NSTEP,NACCPT,NREJCT)
      WORK(7)=H
      IWORK(17)=NFCN
      IWORK(18)=NSTEP
      IWORK(19)=NACCPT
      IWORK(20)=NREJCT
! ----------- RETURN -----------
      RETURN
      END
!
!
!
!  ----- ... AND HERE IS THE CORE INTEGRATOR  ----------
!
      SUBROUTINE DP86CO(N,FCN,X,Y,XEND,HMAX,H,RTOL,ATOL,ITOL,IPRINT,  &
     &   SOLOUT,IOUT,IDID,NMAX,UROUND,METH,NSTIFF,SAFE,BETA,FAC1,FAC2,&
     &   K1,K2,K3,K4,K5,K6,K7,K8,K9,K10,Y1,CONT,ICOMP,NRD,RPAR,IPAR,  &
     &   NFCN,NSTEP,NACCPT,NREJCT)
! ----------------------------------------------------------
!     CORE INTEGRATOR FOR DOP853
!     PARAMETERS SAME AS IN DOP853 WITH WORKSPACE ADDED 
! ---------------------------------------------------------- 
!         DECLARATIONS 
! ---------------------------------------------------------- 
!!!   IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!!!   IMPLICIT INTEGER (I-N)
      implicit none
      double precision, parameter ::                 &
     &  c2  = 0.526001519587677318785587544488D-01,  & 
     &  c3  = 0.789002279381515978178381316732D-01,  &
     &  c4  = 0.118350341907227396726757197510D+00,  &
     &  c5  = 0.281649658092772603273242802490D+00,  &
     &  c6  = 0.333333333333333333333333333333D+00,  &
     &  c7  = 0.25D+00,                              &
     &  c8  = 0.307692307692307692307692307692D+00,  &
     &  c9  = 0.651282051282051282051282051282D+00,  &
     &  c10 = 0.6D+00,                               &
     &  c11 = 0.857142857142857142857142857142D+00,  &
     &  c14 = 0.1D+00,                               &
     &  c15 = 0.2D+00,                               &
     &  c16 = 0.777777777777777777777777777778D+00
      double precision, parameter ::                 &
     &  b1 =   5.42937341165687622380535766363D-2,   &
     &  b6 =   4.45031289275240888144113950566D0,    &
     &  b7 =   1.89151789931450038304281599044D0,    &
     &  b8 =  -5.8012039600105847814672114227D0,     &
     &  b9 =   3.1116436695781989440891606237D-1,    &
     &  b10 = -1.52160949662516078556178806805D-1,   &
     &  b11 =  2.01365400804030348374776537501D-1,   &
     &  b12 =  4.47106157277725905176885569043D-2
      double precision, parameter ::                  &
     &  bhh1 = 0.244094488188976377952755905512D+00,  &
     &  bhh2 = 0.733846688281611857341361741547D+00,  &
     &  bhh3 = 0.220588235294117647058823529412D-01   
      double precision, parameter ::                 &
     &  er1 =  0.1312004499419488073250102996D-01,   &
     &  er6 = -0.1225156446376204440720569753D+01,   &
     &  er7 = -0.4957589496572501915214079952D+00,   &
     &  er8 =  0.1664377182454986536961530415D+01,   &
     &  er9 = -0.3503288487499736816886487290D+00,   &
     &  er10 =  0.3341791187130174790297318841D+00,   &
     &  er11 =  0.8192320648511571246570742613D-01,   &
     &  er12 = -0.2235530786388629525884427845D-01
      double precision, parameter ::                   &
     &  a21 =    5.26001519587677318785587544488D-2,   &
     &  a31 =    1.97250569845378994544595329183D-2,   &
     &  a32 =    5.91751709536136983633785987549D-2,   &
     &  a41 =    2.95875854768068491816892993775D-2,   &
     &  a43 =    8.87627564304205475450678981324D-2,   &
     &  a51 =    2.41365134159266685502369798665D-1,   &
     &  a53 =   -8.84549479328286085344864962717D-1,   &
     &  a54 =    9.24834003261792003115737966543D-1,   &
     &  a61 =    3.7037037037037037037037037037D-2,    &
     &  a64 =    1.70828608729473871279604482173D-1,   &
     &  a65 =    1.25467687566822425016691814123D-1,   &
     &  a71 =    3.7109375D-2,                         &
     &  a74 =    1.70252211019544039314978060272D-1,   &
     &  a75 =    6.02165389804559606850219397283D-2,   &
     &  a76 =   -1.7578125D-2
      double precision, parameter ::                   &
     &  a81 =    3.70920001185047927108779319836D-2,   &
     &  a84 =    1.70383925712239993810214054705D-1,   &
     &  a85 =    1.07262030446373284651809199168D-1,   &
     &  a86 =   -1.53194377486244017527936158236D-2,   &
     &  a87 =    8.27378916381402288758473766002D-3,   &
     &  a91 =    6.24110958716075717114429577812D-1,   &
     &  a94 =   -3.36089262944694129406857109825D0,    &
     &  a95 =   -8.68219346841726006818189891453D-1,   &
     &  a96 =    2.75920996994467083049415600797D1,    &
     &  a97 =    2.01540675504778934086186788979D1,    &
     &  a98 =   -4.34898841810699588477366255144D1,    &
     &  a101 =   4.77662536438264365890433908527D-1,   &
     &  a104 =  -2.48811461997166764192642586468D0,    &
     &  a105 =  -5.90290826836842996371446475743D-1,   &
     &  a106 =   2.12300514481811942347288949897D1,    &
     &  a107 =   1.52792336328824235832596922938D1,    &
     &  a108 =  -3.32882109689848629194453265587D1,    &
     &  a109 =  -2.03312017085086261358222928593D-2
      double precision, parameter ::                  &
     &  a111 =  -9.3714243008598732571704021658D-1,   &
     &  a114 =   5.18637242884406370830023853209D0,   &
     &  a115 =   1.09143734899672957818500254654D0,   &
     &  a116 =  -8.14978701074692612513997267357D0,   &
     &  a117 =  -1.85200656599969598641566180701D1,   &
     &  a118 =   2.27394870993505042818970056734D1,   &
     &  a119 =   2.49360555267965238987089396762D0,   &
     &  a1110 = -3.0467644718982195003823669022D0,    &
     &  a121 =   2.27331014751653820792359768449D0,   &
     &  a124 =  -1.05344954667372501984066689879D1,   &
     &  a125 =  -2.00087205822486249909675718444D0,   &
     &  a126 =  -1.79589318631187989172765950534D1,   &
     &  a127 =   2.79488845294199600508499808837D1,   &
     &  a128 =  -2.85899827713502369474065508674D0,   &
     &  a129 =  -8.87285693353062954433549289258D0,   &
     &  a1210 =  1.23605671757943030647266201528D1,   &
     &  a1211 =  6.43392746015763530355970484046D-1
      double precision, parameter ::                  &
     &  a141 =  5.61675022830479523392909219681D-2,   &
     &  a147 =  2.53500210216624811088794765333D-1,   &
     &  a148 = -2.46239037470802489917441475441D-1,   &
     &  a149 = -1.24191423263816360469010140626D-1,   &
     &  a1410 =  1.5329179827876569731206322685D-1,   &
     &  a1411 =  8.20105229563468988491666602057D-3,  &
     &  a1412 =  7.56789766054569976138603589584D-3,  &
     &  a1413 = -8.298D-3
      double precision, parameter ::                   &
     &  a151 =  3.18346481635021405060768473261D-2,    &
     &  a156 =  2.83009096723667755288322961402D-2,    &
     &  a157 =  5.35419883074385676223797384372D-2,    &
     &  a158 = -5.49237485713909884646569340306D-2,    &
     &  a1511 = -1.08347328697249322858509316994D-4,   &
     &  a1512 =  3.82571090835658412954920192323D-4,   &
     &  a1513 = -3.40465008687404560802977114492D-4,   &
     &  a1514 =  1.41312443674632500278074618366D-1,   &
     &  a161 = -4.28896301583791923408573538692D-1,    &
     &  a166 = -4.69762141536116384314449447206D0,     &
     &  a167 =  7.68342119606259904184240953878D0,     &
     &  a168 =  4.06898981839711007970213554331D0,     &
     &  a169 =  3.56727187455281109270669543021D-1,    &
     &  a1613 = -1.39902416515901462129418009734D-3,   &
     &  a1614 =  2.9475147891527723389556272149D0,     &
     &  a1615 = -9.15095847217987001081870187138D0
      double precision, parameter ::                   &
     &  d41  = -0.84289382761090128651353491142D+01,   &
     &  d46  =  0.56671495351937776962531783590D+00,   &
     &  d47  = -0.30689499459498916912797304727D+01,   &
     &  d48  =  0.23846676565120698287728149680D+01,   &
     &  d49  =  0.21170345824450282767155149946D+01,   &
     &  d410 = -0.87139158377797299206789907490D+00,   &
     &  d411 =  0.22404374302607882758541771650D+01,   &
     &  d412 =  0.63157877876946881815570249290D+00,   &
     &  d413 = -0.88990336451333310820698117400D-01,   &
     &  d414 =  0.18148505520854727256656404962D+02,   &
     &  d415 = -0.91946323924783554000451984436D+01,   &
     &  d416 = -0.44360363875948939664310572000D+01
      double precision, parameter ::                   &
     &  d51  =  0.10427508642579134603413151009D+02,   &
     &  d56  =  0.24228349177525818288430175319D+03,   &
     &  d57  =  0.16520045171727028198505394887D+03,   &
     &  d58  = -0.37454675472269020279518312152D+03,   &
     &  d59  = -0.22113666853125306036270938578D+02,   &
     &  d510 =  0.77334326684722638389603898808D+01,   &
     &  d511 = -0.30674084731089398182061213626D+02,   &
     &  d512 = -0.93321305264302278729567221706D+01,   &
     &  d513 =  0.15697238121770843886131091075D+02,   &
     &  d514 = -0.31139403219565177677282850411D+02,   &
     &  d515 = -0.93529243588444783865713862664D+01,   &
     &  d516 =  0.35816841486394083752465898540D+02
      double precision, parameter ::                   &
     &  d61 =  0.19985053242002433820987653617D+02,    &
     &  d66 = -0.38703730874935176555105901742D+03,    &
     &  d67 = -0.18917813819516756882830838328D+03,    &
     &  d68 =  0.52780815920542364900561016686D+03,    &
     &  d69 = -0.11573902539959630126141871134D+02,    &
     &  d610 =  0.68812326946963000169666922661D+01,   &
     &  d611 = -0.10006050966910838403183860980D+01,   &
     &  d612 =  0.77771377980534432092869265740D+00,   &
     &  d613 = -0.27782057523535084065932004339D+01,   &
     &  d614 = -0.60196695231264120758267380846D+02,   &
     &  d615 =  0.84320405506677161018159903784D+02,   &
     &  d616 =  0.11992291136182789328035130030D+02
      double precision, parameter ::                   &
     &  d71  = -0.25693933462703749003312586129D+02,   &
     &  d76  = -0.15418974869023643374053993627D+03,   &
     &  d77  = -0.23152937917604549567536039109D+03,   &
     &  d78  =  0.35763911791061412378285349910D+03,   &
     &  d79  =  0.93405324183624310003907691704D+02,   &
     &  d710 = -0.37458323136451633156875139351D+02,   &
     &  d711 =  0.10409964950896230045147246184D+03,   &
     &  d712 =  0.29840293426660503123344363579D+02,   &
     &  d713 = -0.43533456590011143754432175058D+02,   &
     &  d714 =  0.96324553959188282948394950600D+02,   &
     &  d715 = -0.39177261675615439165231486172D+02,   &
     &  d716 = -0.14972683625798562581422125276D+03
      double precision :: x,xend,work,rpar
      integer :: n,itol,iout,lwork,iwork,liwork,ipar,idid,nrd
      double precision :: uround,safe,fac1,fac2,beta,hmax,h,cont
      integer :: nfcn,nstep,naccpt,nrejct,iprint,nmax,meth,  &
                 nstiff,nrdens,iek1,iek2,iek3,iek4,iek5,iek6,&
                 iek7,iek8,iek9,iek10,iey1,ieco,istore,icomp
      DOUBLE PRECISION Y(N),Y1(N),K1(N),K2(N),K3(N),K4(N),K5(N),K6(N)
      DOUBLE PRECISION K7(N),K8(N),K9(N),K10(N),ATOL(*),RTOL(*)     
      DIMENSION CONT(8*NRD),ICOMP(NRD),RPAR(*),IPAR(*)
      double precision :: facold,expo1,facc1,facc2,posneg,atoli, &
         rtoli,hlamb,xold,xph,err,err2,sk,erri,deno,fac11,fac,   &
         hnew,stnum,stden,ydiff,bspl,hout,xout
      integer :: i,j,iasti,iord,irtrn,nonsti
      LOGICAL REJECT,LAST,EVENT
      EXTERNAL FCN
      external solout !!! Added by RMP, 15/09/2022
      COMMON /CONDO8/XOLD,HOUT
! *** *** *** *** *** *** ***
!  INITIALISATIONS
! *** *** *** *** *** *** *** 
      FACOLD=1.D-4  
      EXPO1=1.d0/8.d0-BETA*0.2D0
      FACC1=1.D0/FAC1
      FACC2=1.D0/FAC2
      POSNEG=SIGN(1.D0,XEND-X) 
! --- INITIAL PREPARATIONS   
      ATOLI=ATOL(1)
      RTOLI=RTOL(1)    
      LAST=.FALSE. 
      HLAMB=0.D0
      IASTI=0
      CALL mbe_rhs(X,Y,K1) !!! CALL FCN(N,X,Y,K1,RPAR,IPAR)
      HMAX=ABS(HMAX)     
      IORD=8  
      IF (H.EQ.0.D0) H=HINIT(N,FCN,X,Y,XEND,POSNEG,K1,K2,K3,IORD,    &
     &                       HMAX,ATOL,RTOL,ITOL,RPAR,IPAR)
      NFCN=NFCN+2
      REJECT=.FALSE.
      XOLD=X
      IF (IOUT.NE.0) THEN 
!!!!!
!!! Original code:
!!        IRTRN=1 
!!        HOUT=1.D0
!!        CALL SOLOUT(NACCPT+1,XOLD,X,Y,N,CONT,ICOMP,NRD,            &
!!   &                RPAR,IPAR,IRTRN,XOUT)
!!        IF (IRTRN.LT.0) GOTO 79
!!! Current version, for consistency with the current version of the
!!! mbe module:
          print*,'DP86CO was called for a non-zero value of IOUT'
          call mbe_error
!!!!!
      END IF
! --- BASIC INTEGRATION STEP  
   1  CONTINUE
      IF (NSTEP.GT.NMAX) GOTO 78
      IF (0.1D0*ABS(H).LE.ABS(X)*UROUND)GOTO 77
      IF ((X+1.01D0*H-XEND)*POSNEG.GT.0.D0) THEN
         H=XEND-X 
         LAST=.TRUE.
      END IF
      NSTEP=NSTEP+1
! --- THE TWELVE STAGES
      IF (IRTRN.GE.2) THEN
         CALL mbe_rhs(X,Y,K1) !!! CALL FCN(N,X,Y,K1,RPAR,IPAR)
      END IF
      DO 22 I=1,N 
      Y1(I)=Y(I)+H*A21*K1(I)  
  22  continue
      CALL mbe_rhs(X+C2*H,Y1,K2) !!! CALL FCN(N,X+C2*H,Y1,K2,RPAR,IPAR)
      DO 23 I=1,N 
      Y1(I)=Y(I)+H*(A31*K1(I)+A32*K2(I))  
  23  continue
      CALL mbe_rhs(X+C3*H,Y1,K3) !!! CALL FCN(N,X+C3*H,Y1,K3,RPAR,IPAR)
      DO 24 I=1,N 
      Y1(I)=Y(I)+H*(A41*K1(I)+A43*K3(I))  
  24  continue
      CALL mbe_rhs(X+C4*H,Y1,K4) !!! CALL FCN(N,X+C4*H,Y1,K4,RPAR,IPAR)
      DO 25 I=1,N 
      Y1(I)=Y(I)+H*(A51*K1(I)+A53*K3(I)+A54*K4(I))
  25  continue
      CALL mbe_rhs(X+C5*H,Y1,K5) !!! CALL FCN(N,X+C5*H,Y1,K5,RPAR,IPAR)
      DO 26 I=1,N 
      Y1(I)=Y(I)+H*(A61*K1(I)+A64*K4(I)+A65*K5(I))
  26  continue
      CALL mbe_rhs(X+C6*H,Y1,K6) !!! CALL FCN(N,X+C6*H,Y1,K6,RPAR,IPAR)
      DO 27 I=1,N 
      Y1(I)=Y(I)+H*(A71*K1(I)+A74*K4(I)+A75*K5(I)+A76*K6(I))
  27  continue
      CALL mbe_rhs(X+C7*H,Y1,K7) !!! CALL FCN(N,X+C7*H,Y1,K7,RPAR,IPAR)
      DO 28 I=1,N 
      Y1(I)=Y(I)+H*(A81*K1(I)+A84*K4(I)+A85*K5(I)+A86*K6(I)+A87*K7(I))  
  28  continue
      CALL mbe_rhs(X+C8*H,Y1,K8) !!! CALL FCN(N,X+C8*H,Y1,K8,RPAR,IPAR)
      DO 29 I=1,N 
      Y1(I)=Y(I)+H*(A91*K1(I)+A94*K4(I)+A95*K5(I)+A96*K6(I)+A97*K7(I)  &
     &   +A98*K8(I))
  29  continue
      CALL mbe_rhs(X+C9*H,Y1,K9) !!! CALL FCN(N,X+C9*H,Y1,K9,RPAR,IPAR)
      DO 30 I=1,N 
      Y1(I)=Y(I)+H*(A101*K1(I)+A104*K4(I)+A105*K5(I)+A106*K6(I)        &
     &   +A107*K7(I)+A108*K8(I)+A109*K9(I))
  30  continue
      CALL mbe_rhs(X+C10*H,Y1,K10) !!! CALL FCN(N,X+C10*H,Y1,K10,RPAR,IPAR)
      DO 31 I=1,N 
      Y1(I)=Y(I)+H*(A111*K1(I)+A114*K4(I)+A115*K5(I)+A116*K6(I)        &
     &   +A117*K7(I)+A118*K8(I)+A119*K9(I)+A1110*K10(I))
  31  continue
      CALL mbe_rhs(X+C11*H,Y1,K2) !!! CALL FCN(N,X+C11*H,Y1,K2,RPAR,IPAR)
      XPH=X+H
      DO 32 I=1,N 
      Y1(I)=Y(I)+H*(A121*K1(I)+A124*K4(I)+A125*K5(I)+A126*K6(I)        &
     &   +A127*K7(I)+A128*K8(I)+A129*K9(I)+A1210*K10(I)+A1211*K2(I))
  32  continue
      CALL mbe_rhs(XPH,Y1,K3) !!! CALL FCN(N,XPH,Y1,K3,RPAR,IPAR)
      NFCN=NFCN+11
      DO 35 I=1,N 
      K4(I)=B1*K1(I)+B6*K6(I)+B7*K7(I)+B8*K8(I)+B9*K9(I)               &
     &   +B10*K10(I)+B11*K2(I)+B12*K3(I)
      K5(I)=Y(I)+H*K4(I)
  35  continue
! --- ERROR ESTIMATION  
      ERR=0.D0
      ERR2=0.D0
      IF (ITOL.EQ.0) THEN   
        DO 41 I=1,N 
        SK=ATOLI+RTOLI*MAX(ABS(Y(I)),ABS(K5(I)))
        ERRI=K4(I)-BHH1*K1(I)-BHH2*K9(I)-BHH3*K3(I)
        ERR2=ERR2+(ERRI/SK)**2
        ERRI=ER1*K1(I)+ER6*K6(I)+ER7*K7(I)+ER8*K8(I)+ER9*K9(I)        &
     &      +ER10*K10(I)+ER11*K2(I)+ER12*K3(I)
        ERR=ERR+(ERRI/SK)**2
  41    continue
      ELSE
        DO 42 I=1,N 
        SK=ATOL(I)+RTOL(I)*MAX(ABS(Y(I)),ABS(K5(I)))
        ERRI=K4(I)-BHH1*K1(I)-BHH2*K9(I)-BHH3*K3(I)
        ERR2=ERR2+(ERRI/SK)**2
        ERRI=ER1*K1(I)+ER6*K6(I)+ER7*K7(I)+ER8*K8(I)+ER9*K9(I)        &
     &      +ER10*K10(I)+ER11*K2(I)+ER12*K3(I)
        ERR=ERR+(ERRI/SK)**2
  42    continue
      END IF  
      DENO=ERR+0.01D0*ERR2
      IF (DENO.LE.0.D0) DENO=1.D0
      ERR=ABS(H)*ERR*SQRT(1.D0/(N*DENO))
! --- COMPUTATION OF HNEW
      FAC11=ERR**EXPO1
! --- LUND-STABILIZATION
      FAC=FAC11/FACOLD**BETA
! --- WE REQUIRE  FAC1 <= HNEW/H <= FAC2
      FAC=MAX(FACC2,MIN(FACC1,FAC/SAFE))
      HNEW=H/FAC  
      IF(ERR.LE.1.D0)THEN
! --- STEP IS ACCEPTED  
         FACOLD=MAX(ERR,1.0D-4)
         NACCPT=NACCPT+1
         CALL mbe_rhs(XPH,K5,K4) !!! CALL FCN(N,XPH,K5,K4,RPAR,IPAR)
         NFCN=NFCN+1
! ------- STIFFNESS DETECTION
         IF (MOD(NACCPT,NSTIFF).EQ.0.OR.IASTI.GT.0) THEN
            STNUM=0.D0
            STDEN=0.D0
            DO 64 I=1,N 
               STNUM=STNUM+(K4(I)-K3(I))**2
               STDEN=STDEN+(K5(I)-Y1(I))**2
 64         CONTINUE  
            IF (STDEN.GT.0.D0) HLAMB=ABS(H)*SQRT(STNUM/STDEN) 
            IF (HLAMB.GT.6.1D0) THEN
               NONSTI=0
               IASTI=IASTI+1  
               IF (IASTI.EQ.15) THEN
                  IF (IPRINT.GT.0) WRITE (IPRINT,*)                   &
     &               ' THE PROBLEM SEEMS TO BECOME STIFF AT X = ',X   
                  IF (IPRINT.LE.0) GOTO 76
               END IF
            ELSE
               NONSTI=NONSTI+1  
               IF (NONSTI.EQ.6) IASTI=0
            END IF
         END IF 
! ------- FINAL PREPARATION FOR DENSE OUTPUT
         EVENT=(IOUT.EQ.3).AND.(XOUT.LE.XPH)
         IF (IOUT.EQ.2.OR.EVENT) THEN
! ----    SAVE THE FIRST FUNCTION EVALUATIONS   
            DO 62 J=1,NRD
               I=ICOMP(J)
               CONT(J)=Y(I)
               YDIFF=K5(I)-Y(I)
               CONT(J+NRD)=YDIFF
               BSPL=H*K1(I)-YDIFF
               CONT(J+NRD*2)=BSPL
               CONT(J+NRD*3)=YDIFF-H*K4(I)-BSPL
               CONT(J+NRD*4)=D41*K1(I)+D46*K6(I)+D47*K7(I)+D48*K8(I)  &
     &                  +D49*K9(I)+D410*K10(I)+D411*K2(I)+D412*K3(I)
               CONT(J+NRD*5)=D51*K1(I)+D56*K6(I)+D57*K7(I)+D58*K8(I)  &
     &                  +D59*K9(I)+D510*K10(I)+D511*K2(I)+D512*K3(I)
               CONT(J+NRD*6)=D61*K1(I)+D66*K6(I)+D67*K7(I)+D68*K8(I)  &
     &                  +D69*K9(I)+D610*K10(I)+D611*K2(I)+D612*K3(I)
               CONT(J+NRD*7)=D71*K1(I)+D76*K6(I)+D77*K7(I)+D78*K8(I)  &
     &                  +D79*K9(I)+D710*K10(I)+D711*K2(I)+D712*K3(I)
   62       CONTINUE 
! ---     THE NEXT THREE FUNCTION EVALUATIONS
            DO 51 I=1,N 
               Y1(I)=Y(I)+H*(A141*K1(I)+A147*K7(I)+A148*K8(I)         &
     &            +A149*K9(I)+A1410*K10(I)+A1411*K2(I)+A1412*K3(I)    &
     &            +A1413*K4(I))    
  51        continue
            CALL mbe_rhs(X+C14*H,Y1,K10)!!!CALL FCN(N,X+C14*H,Y1,K10,RPAR,IPAR)
            DO 52 I=1,N 
               Y1(I)=Y(I)+H*(A151*K1(I)+A156*K6(I)+A157*K7(I)         &
     &            +A158*K8(I)+A1511*K2(I)+A1512*K3(I)+A1513*K4(I)     &
     &            +A1514*K10(I))
  52        continue
            CALL mbe_rhs(X+C15*H,Y1,K2)!!!CALL FCN(N,X+C15*H,Y1,K2,RPAR,IPAR)
            DO 53 I=1,N 
               Y1(I)=Y(I)+H*(A161*K1(I)+A166*K6(I)+A167*K7(I)         &
     &            +A168*K8(I)+A169*K9(I)+A1613*K4(I)+A1614*K10(I)     &
     &            +A1615*K2(I))                                        
  53        continue
            CALL mbe_rhs(X+C16*H,Y1,K3)!!!CALL FCN(N,X+C16*H,Y1,K3,RPAR,IPAR)
            NFCN=NFCN+3 
! ---     FINAL PREPARATION
            DO 63 J=1,NRD
               I=ICOMP(J)
               CONT(J+NRD*4)=H*(CONT(J+NRD*4)+D413*K4(I)+D414*K10(I)  &
     &            +D415*K2(I)+D416*K3(I))
               CONT(J+NRD*5)=H*(CONT(J+NRD*5)+D513*K4(I)+D514*K10(I)  &
     &            +D515*K2(I)+D516*K3(I))
               CONT(J+NRD*6)=H*(CONT(J+NRD*6)+D613*K4(I)+D614*K10(I)  &
     &            +D615*K2(I)+D616*K3(I))
               CONT(J+NRD*7)=H*(CONT(J+NRD*7)+D713*K4(I)+D714*K10(I)  &
     &            +D715*K2(I)+D716*K3(I))
  63        CONTINUE
            HOUT=H
         END IF
         DO 67 I=1,N
         K1(I)=K4(I)
         Y(I)=K5(I)
  67     continue
         XOLD=X
         X=XPH
!!!!!
!!! The following part of the original code has been commented out
!!       IF (IOUT.EQ.1.OR.IOUT.EQ.2.OR.EVENT) THEN
!!          CALL SOLOUT(NACCPT+1,XOLD,X,Y,N,CONT,ICOMP,NRD,           &
!!   &                  RPAR,IPAR,IRTRN,XOUT)
!!          IF (IRTRN.LT.0) GOTO 79
!!       END IF 
!!!!!
! ------- NORMAL EXIT
         IF (LAST) THEN
            H=HNEW
            IDID=1
            RETURN
         END IF
         IF(ABS(HNEW).GT.HMAX)HNEW=POSNEG*HMAX  
         IF(REJECT)HNEW=POSNEG*MIN(ABS(HNEW),ABS(H))
         REJECT=.FALSE. 
      ELSE  
! --- STEP IS REJECTED   
         HNEW=H/MIN(FACC1,FAC11/SAFE)
         REJECT=.TRUE.  
         IF(NACCPT.GE.1)NREJCT=NREJCT+1   
         LAST=.FALSE.
      END IF
      H=HNEW
      GOTO 1
! --- FAIL EXIT
  76  CONTINUE
      IDID=-4
      RETURN
  77  CONTINUE
      IF (IPRINT.GT.0) WRITE(IPRINT,979)X   
      IF (IPRINT.GT.0) WRITE(IPRINT,*)' STEP SIZE TOO SMALL, H=',H
      IDID=-3
      RETURN
  78  CONTINUE
      IF (IPRINT.GT.0) WRITE(IPRINT,979)X   
      IF (IPRINT.GT.0) WRITE(IPRINT,*)                                &
     &     ' MORE THAN NMAX =',NMAX,'STEPS ARE NEEDED' 
      IDID=-2
      RETURN
  79  CONTINUE
      IF (IPRINT.GT.0) WRITE(IPRINT,979)X
 979  FORMAT(' EXIT OF DOP853 AT X=',E18.4) 
      IDID=2
      RETURN
      END
!
      FUNCTION HINIT(N,FCN,X,Y,XEND,POSNEG,F0,F1,Y1,IORD,             &
     &                       HMAX,ATOL,RTOL,ITOL,RPAR,IPAR)
! ----------------------------------------------------------
! ----  COMPUTATION OF AN INITIAL STEP SIZE GUESS
! ----------------------------------------------------------
!!!   IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!!!   IMPLICIT INTEGER (I-N)
      implicit none
      double precision :: x,y,xend,posneg,f0,f1,y1,hmax,atol,    &
                          rtol,rpar
      integer :: n,iord,itol,ipar
      DIMENSION Y(N),Y1(N),F0(N),F1(N),ATOL(*),RTOL(*)
      DIMENSION RPAR(*),IPAR(*)
      double precision :: dnf,dny,atoli,rtoli,sk,h,der2,h1,hinit,der12
      integer :: i
      external fcn !!! Added by RMP 15/09/2022
! ---- COMPUTE A FIRST GUESS FOR EXPLICIT EULER AS
! ----   H = 0.01 * NORM (Y0) / NORM (F0)
! ---- THE INCREMENT FOR EXPLICIT EULER IS SMALL
! ---- COMPARED TO THE SOLUTION
      DNF=0.0D0
      DNY=0.0D0 
      ATOLI=ATOL(1)
      RTOLI=RTOL(1)    
      IF (ITOL.EQ.0) THEN   
        DO 10 I=1,N 
        SK=ATOLI+RTOLI*ABS(Y(I))
        DNF=DNF+(F0(I)/SK)**2
        DNY=DNY+(Y(I)/SK)**2 
  10    continue
      ELSE
        DO 11 I=1,N 
        SK=ATOL(I)+RTOL(I)*ABS(Y(I))
        DNF=DNF+(F0(I)/SK)**2
        DNY=DNY+(Y(I)/SK)**2 
  11    continue
      END IF
      IF (DNF.LE.1.D-10.OR.DNY.LE.1.D-10) THEN
         H=1.0D-6
      ELSE
         H=SQRT(DNY/DNF)*0.01D0 
      END IF
      H=MIN(H,HMAX)
      H=SIGN(H,POSNEG) 
! ---- PERFORM AN EXPLICIT EULER STEP
      DO 12 I=1,N
      Y1(I)=Y(I)+H*F0(I)
  12  continue
      CALL mbe_rhs(X+H,Y1,F1) !!! CALL FCN(N,X+H,Y1,F1,RPAR,IPAR) 
! ---- ESTIMATE THE SECOND DERIVATIVE OF THE SOLUTION
      DER2=0.0D0
      IF (ITOL.EQ.0) THEN   
        DO 15 I=1,N 
        SK=ATOLI+RTOLI*ABS(Y(I))
        DER2=DER2+((F1(I)-F0(I))/SK)**2   
  15    continue
      ELSE
        DO 16 I=1,N 
        SK=ATOL(I)+RTOL(I)*ABS(Y(I))
        DER2=DER2+((F1(I)-F0(I))/SK)**2   
  16    continue
      END IF
      DER2=SQRT(DER2)/H
! ---- STEP SIZE IS COMPUTED SUCH THAT
! ----  H**IORD * MAX ( NORM (F0), NORM (DER2)) = 0.01
      DER12=MAX(ABS(DER2),SQRT(DNF))
      IF (DER12.LE.1.D-15) THEN
         H1=MAX(1.0D-6,ABS(H)*1.0D-3)
      ELSE
         H1=(0.01D0/DER12)**(1.D0/IORD) 
      END IF
      H=MIN(100*ABS(H),H1,HMAX)
      HINIT=SIGN(H,POSNEG)  
      RETURN
      END 
!
      FUNCTION CONTD8(II,X,CON,ICOMP,ND)
! ----------------------------------------------------------
!     THIS FUNCTION CAN BE USED FOR CONINUOUS OUTPUT IN CONNECTION
!     WITH THE OUTPUT-SUBROUTINE FOR DOP853. IT PROVIDES AN
!     APPROXIMATION TO THE II-TH COMPONENT OF THE SOLUTION AT X.
! ----------------------------------------------------------
!!!   IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!!!   IMPLICIT INTEGER (I-N)
      implicit none
      double precision :: x,con,xold,h
      integer :: ii,icomp,nd
      DIMENSION CON(8*ND),ICOMP(ND)
      double precision :: s,s1,conpar,contd8
      integer :: i,j
      COMMON /CONDO8/XOLD,H
! ----- COMPUTE PLACE OF II-TH COMPONENT 
      I=0 
      DO 5 J=1,ND 
      IF (ICOMP(J).EQ.II) I=J
   5  CONTINUE
      IF (I.EQ.0) THEN
         WRITE (6,*) ' NO DENSE OUTPUT AVAILABLE FOR COMP.',II 
         RETURN
      END IF  
      S=(X-XOLD)/H
      S1=1.D0-S
      CONPAR=CON(I+ND*4)+S*(CON(I+ND*5)+S1*(CON(I+ND*6)+S*CON(I+ND*7)))
      CONTD8=CON(I)+S*(CON(I+ND)+S1*(CON(I+ND*2)+S*(CON(I+ND*3)        &
     &        +S1*CONPAR)))
      RETURN
      END

!!!!!!!!!!!!!!!!!!!
 
      end subroutine mbe_intobe

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine mbe_propagate_2(istart,popinit,imethod,nsubsteps,         &
               izrule,zmax_in,n_z_steps,fieldpr,fieldcp,                   &
               density_in,iDoppler,nz_writeout,nt_writeout,output_pr,      &
               icont,nsb_var)
!
!  icont (optional) : If .eq. 1, simply re-use the existing ommats arrays.
!  If .ne. 1, reset fl_init_ommats to .true. to force the re-calculation of
!  these arrays.
!
!  The time mesh (in the global array t_mesh) is expected to be specified
!  in microseconds.
!
!  iDoppler.eq.1 is the signal that the calculation involved an averaging
!  over velocity classes. The calculation then uses the abcissas and weights
!  (multiplied by the Maxwell distribution) contained in the arrays
!  vmesh and fMvweight of the obe module. The variable iDoppler should be
!  set to 0 if no Doppler averaging is intended.
!
!  output_pr is the name of an external subroutine which mbe_propagate_2 calls
!  for writing out the results of the calculation (the density matrix and/or
!  the calculted fields). *** This subroutine should not alter the content of 
!  any of its arguments. *** The variables nz_writeout and nt_writeout control
!  the mesh points at which output_pr is called.
!
!  imethod determines whether the integration over time is done using
!  a fourth or a fifth order Runge-Kutta algorithm (respectively
!  imethod.eq.1 or imethod.eq.5).
!
!  The time-dependent applied fields must have been set by a call
!  either to mbe_set_tdfields_A (with iinterp.eq.1) or to mbe_set_tdfields_B
!  prior to a call to mbe_propagate_2.
!
!  nsb_var is an optional argument, which must be present if the value of
!  nsubsteps is 0. In this case, the number of intermediate steps in the time
!  integration between t(j-1) and t(j) is taken to be nsb_var(j). 
!  If nsubsteps is > 0, then the number of steps between t(j-1) and t(j) is
!  taken to be nsubsteps, whether nsb_var is present or not.
!
      use obe_constants
      implicit none
      integer, parameter :: neqs = nst*nst
      double precision, dimension (nst), intent(in) :: popinit
      double precision, intent(in) :: density_in,zmax_in
      integer, dimension (:), intent(in), optional :: nsb_var
      integer, intent(in) :: imethod,istart,izrule,nsubsteps,       &
         nt_writeout,nz_writeout, n_z_steps, iDoppler
      integer :: npr,ncp
      integer, intent(in), optional :: icont
      type(obecfield), intent(in) :: fieldpr,fieldcp
      real(kd), dimension(:), allocatable :: z_mesh
      real(kd) :: density,zmax
      real(kd) :: atten_p,atten_c,t,time_step,zmin
      real(kd), dimension (:), allocatable :: Fd_p_re_prev
      real(kd), dimension (:), allocatable :: Fd_p_im_prev
      real(kd), dimension (:), allocatable :: Fd_c_re_prev
      real(kd), dimension (:), allocatable :: Fd_c_im_prev
      real(kd), dimension (:), allocatable :: alpha_p,refr_index_p
      real(kd), dimension (:), allocatable :: alpha_c,refr_index_c
      real(kd), dimension (:,:), allocatable :: rhs_mat
      real(kd), dimension (:,:), allocatable :: rhsprre,rhsprim
      real(kd), dimension (:,:), allocatable :: rhscpre,rhscpim
      real(kd), dimension (:,:), allocatable :: sumrho_pre,sumrho_pim
      real(kd), dimension (:,:), allocatable :: sumrho_cre,sumrho_cim
      real(kd), dimension(:), allocatable :: factpr,factcp
      integer, dimension (:), allocatable :: idxprre,idxprim,idxcpre,idxcpim
      real(kd), dimension (neqs) :: rhovec
      real(kd) :: akpr,akcp,hzennglc,hzennglp,unit_l,pi,z_step
      integer, dimension (:), allocatable :: nsbsteps
      integer :: i,istat,iv,ivmax,j,jmin,k,kmin,m,mrl,mim,nvsize
      logical, save :: fl_Doppler_old = .false.
      external :: output_pr
!
      if(.not.fl_spline)then
         print*,'mbe_propagate_2 requires interpolation...'
         call mbe_error
      endif
!
!  Initialize variables controlling the integration, if dop853 is
!  to be used. rtol_dop853, atol_dop853 and fl_set_tol_dop853 are
!  global variables.
      if(imethod.eq.4)then
         call obe_get_tol_dop853(rtol_dop853,atol_dop853,fl_set_tol_dop853)
      endif
!
!  Set fl_twofields to .true. as this subroutine is for a pair of fields.
      fl_twofields = .true.
!
      if(istart.lt.0 .or. istart.gt.2)then
         print*,'mbe_propagate_2: illegal value of istart. ',istart
         call mbe_error
      endif
!
!  Chek that nsubsteps and nsb_var (if present) make sense.
      if(nsubsteps.eq.-1)then
         if(present(nsb_var))then
            if(size(nsb_var).lt.n_time_steps)then
               print*,'mbe_propagate_2: the size of nsb_var is inconsistent'
               print*,'with the value of n_time_steps passed'
               print*,'to mbe_set_tdfields.'
               call mbe_error
            else
               kmin = lbound(nsb_var,1)
            endif
         else
            print*,'mbe_propagate_2: nsubsteps = -1 but nsb_var is not present.'
            call mbe_error
         endif
      elseif(nsubsteps.lt.1)then
         print*,'mbe_propagate_2: nsubsteps is nonsensical. ',nsubsteps
         call mbe_error
      endif
      if(allocated(nsbsteps))deallocate(nsbsteps)
      allocate(nsbsteps(n_time_steps),stat=istat)
      if(istat.gt.0)then
         print*,'mbe_propagate_2: Error at the allocation of nsbsteps.'
         call mbe_error
      endif
      do k = 1,n_time_steps
         if(nsubsteps.gt.0)then
            nsbsteps(k) = nsubsteps
         else
            nsbsteps(k) = nsb_var(kmin + k - 1)
            if(nsbsteps(k).lt.1)then
               print*,'mbe_propagate_2: Nonsensical value in nsb_var.'
               print*,kmin,k,kmin+k-1,nsb_var(kmin+k-1)
               call mbe_error
            endif
         endif
      enddo
!
!  Convert to kd precision
      density = real(density_in,kd)
      zmax = real(zmax_in,kd)
!
!  Check that the two fields are meant to propagate in the positive
!  z-direction.
      if(fieldpr%idir.ne.1 .or. fieldcp%idir.ne.1)then
         print*,'mbe_propagate_2: both fields must propagate in the'
         print*,'positive z-direction.'
         print*,fieldpr%idir,fieldcp%idir
         call mbe_error
      endif
!
!  Initialise the numerical integration over the velocity classes
!  if necessary.
      if(iDoppler.eq.1)then
         call mbe_set_Doppler
         if(fl_Doppler_notinit)then
            print*,'mbe_propagate_2: attempt at a calculation involving'
            print*,'a Doppler averaging before a call to obe_set_Doppler.'
            call mbe_error
         else
            fl_Doppler = .true.
            ivmax = n_v_values
         endif
      elseif(iDoppler.eq.0)then
         fl_Doppler = .false.
         ivmax = 1
      else
         print*,'mbe_propagate_2: illegal value of iDoppler ',iDoppler
         call mbe_error
      endif
!
!  Initialise rhsmat_0
      call mbe_set_rhsmat_0
!
!  Define the factors used in the calculation
!
!  First for the probe field.
!
!  Start by finding the number of transitions which need to be taken into
!  account and preparing the idx arrays. The content of fieldpr%dip_mom
!  will already have been checked by obe_setfields, which avoids the 
!  possibility that a same transition would be counted twice (once as i,j,
!  once as j,i), but for safety we check again.
      npr = 0
      do i = 1,nst
         do j = 1,nst
            if(fieldpr%dip_mom(j,i).ne.(0.d0,0.d0))then
               if(fieldpr%dip_mom(i,j).ne.(0.d0,0.d0))then
                  print*,'mbe_propagate_2: error for the probe field,'
                  print*,'the array dip_mom is incorrectly formatted.'
                  print*,j,i,fieldpr%dip_mom(j,i),fieldpr%dip_mom(i,j)
                  call mbe_error
               endif
               npr = npr + 1
            endif
         enddo
      enddo
!
      if(allocated(idxprre))deallocate(idxprre)
      if(allocated(idxprim))deallocate(idxprim)
      allocate(idxprre(npr),idxprim(npr),factpr(npr),stat=istat)
      if(istat.gt.0)then
         print*,'mbe_propagate_2: Error at the allocation of the'
         print*,'idxppre, idxprim and factpr arrays.'
         call mbe_error
      endif
!
      npr = 0
      do i = 1,nst
         do j = 1,nst
            if(fieldpr%dip_mom(j,i).ne.(0.d0,0.d0))then
               npr = npr + 1
               if(i.ge.j)then
                  print*,'mbe_propagate_2: error in dip_mom for the probe'
                  print*,'field, or error in the labelling of the states.'
                  print*,'The subroutine will work only if dip_mom(j,i) is'
                  print*,'non-zero only if state j is higher in energy than'
                  print*,'state i and j > i.', i, j
                  call mbe_error
               endif
               call ldbl_coher_index(i,j,idxprre(npr),idxprim(npr))
               if(real(dimag(fieldpr%dip_mom(j,i)),kd).ne.0.0_kd)then
                  print*,'mbe_propagate_2 cannot handle complex'
                  print*,'transitio dipole matrix elements.'
                  print*,j,i,fieldpr%dip_mom(j,i),' (probe field)'
                  call mbe_error
               else
                  factpr(npr) = real(fieldpr%dip_mom(j,i),kd)
               endif
            endif
         enddo
      enddo
!
! Repeat for the coupling field
!
      ncp = 0
      do i = 1,nst
         do j = 1,nst
            if(fieldcp%dip_mom(j,i).ne.(0.d0,0.d0))then
               if(fieldcp%dip_mom(i,j).ne.(0.d0,0.d0))then
                  print*,'mbe_propagate_2: error for the coupling field,'
                  print*,'the array dip_mom is incorrectly formatted.'
                  print*,j,i,fieldcp%dip_mom(j,i),fieldcp%dip_mom(i,j)
                  call mbe_error
               endif
               ncp = ncp + 1
            endif
         enddo
      enddo
!
      if(allocated(idxcpre))deallocate(idxcpre)
      if(allocated(idxcpim))deallocate(idxcpim)
      allocate(idxcpre(ncp),idxcpim(ncp),factcp(ncp),stat=istat)
      if(istat.gt.0)then
         print*,'mbe_propagate_2: Error at the allocation of the'
         print*,'idxcpre, idxcpim and factcp arrays.'
         call mbe_error
      endif
!
      ncp = 0
      do i = 1,nst
         do j = 1,nst
            if(fieldcp%dip_mom(j,i).ne.(0.d0,0.d0))then
               ncp = ncp + 1
               if(i.ge.j)then
                  print*,'mbe_propagate_2: error in dip_mom for the coupling'
                  print*,'field, or error in the labelling of the states.'
                  print*,'The subroutine will work only if dip_mom(j,i) is'
                  print*,'non-zero only if state j is higher in energy than'
                  print*,'state i and j > i.', i, j
                  call mbe_error
               endif
               call ldbl_coher_index(i,j,idxcpre(npr),idxcpim(npr))
               if(real(dimag(fieldcp%dip_mom(j,i)),kd).ne.0.0_kd)then
                  print*,'mbe_propagate_2 cannot handle complex transition dipole'
                  print*,'matrix elements.'
                  print*,j,i,fieldcp%dip_mom(j,i),' (coupling field)'
                  call mbe_error
               else
                  factcp(ncp) = real(fieldcp%dip_mom(j,i),kd)
               endif
            endif
         enddo
      enddo
!
!  Wavenumbers (wavelengths are meant to be specified in nm):
      pi = 4.0_kd*atan(1.0_kd)
      if(fieldpr%wavelength.eq.0.0_kd)then
         print*,'mbe_propagate_2: The wavelength of the probe field is zero.'
         call mbe_error
      endif
      if(fieldcp%wavelength.eq.0.0_kd)then
         print*,'mbe_propagate_2: The wavelength of the coupling field is zero.'
         call mbe_error
      endif
      akpr = 2.0_kd * pi / (fieldpr%wavelength*1.0e-9_kd)
      akcp = 2.0_kd * pi / (fieldcp%wavelength*1.0e-9_kd)
!     
!  Signal to re-initialize the arrays set by mbe_set_rhsmat (the logical
!  variable fl_init_ommats is global within the module), unless we continue
!  with the arrays initialized at a previous run.
      fl_init_ommats = .true.
      if(present(icont))then
         if(icont.eq.1)then
            if(fl_Doppler.neqv.fl_Doppler_old)then
               print*,'mbe_propagate_2: a continuation run with a different'
               print*,'Doppler option is not allowed.'
               call mbe_error
            endif
            fl_init_ommats = .false.
         endif
      endif
      fl_Doppler_old = fl_Doppler
!
      if(nt_writeout.le.0 .or. nz_writeout.le.0)then
         print*,'mbe_propagate_2: the value of nz_writeout or nt_writeout is'
         print*,'nonsensical. ',nz_writeout,nt_writeout
         call mbe_error
      endif
!
!  Check that the complex field amplitudes have been provided (they must
!  be provided through a call to mbe_set_tdfields) and if they have been
!  transcribe them into the relevant arrays.
      if(fl_notinit_fields)then
         print*,'mbe_propagate_2 was called before a call to mbe_set_tdfields.'
         call mbe_error
      endif
!  Write the input fields into the relevant arrays.
      if(allocated(Field_p_re_vec))deallocate(Field_p_re_vec)
      if(allocated(Field_p_im_vec))deallocate(Field_p_im_vec)
      if(allocated(Field_c_re_vec))deallocate(Field_c_re_vec)
      if(allocated(Field_c_im_vec))deallocate(Field_c_im_vec)
      allocate(Field_p_re_vec(0:n_time_steps),                 &
               Field_p_im_vec(0:n_time_steps),                 &
               Field_c_re_vec(0:n_time_steps),                 &
               Field_c_im_vec(0:n_time_steps), stat = istat)
      if(istat.ne.0)then
         print*,'mbe_propagate_2: Error at the allocation of the arrays (1).'
         call mbe_error
      endif
      do k = 0, n_time_steps
         Field_p_re_vec(k) = real(dreal(Field_p_vec(k)),kd)
         Field_p_im_vec(k) = real(dimag(Field_p_vec(k)),kd)
         Field_c_re_vec(k) = real(dreal(Field_c_vec(k)),kd)
         Field_c_im_vec(k) = real(dimag(Field_c_vec(k)),kd)
      enddo
!    
!  Allocate the rest of the storage:
!
      allocate(Fd_p_re_prev(0:n_time_steps),  &
               Fd_p_im_prev(0:n_time_steps),  &
               Fd_c_re_prev(0:n_time_steps),  &
               Fd_c_im_prev(0:n_time_steps),  &
               Fd_pre2(0:n_time_steps),  &
               Fd_pim2(0:n_time_steps),  &
               Fd_cre2(0:n_time_steps),  &
               Fd_cim2(0:n_time_steps),  &
               rhs_mat(neqs,neqs),            &
               stat=istat)
      if(istat.gt.0)then
         print*,'mbe_propagate_2: Error at the allocation of the arrays (2).'
         call mbe_error
      endif
      if(izrule.eq.0)then
         allocate(z_mesh(-1:n_z_steps),           &
                  alpha_p(0:n_time_steps),        &
                  alpha_c(0:n_time_steps),        &
                  refr_index_p(0:n_time_steps),   &
                  refr_index_c(0:n_time_steps),   &
                  stat=istat)
      elseif(izrule.eq.2)then
         allocate(z_mesh(-1:n_z_steps),           &
                  rhsprre(0:n_time_steps,-2:-1),  &
                  rhsprim(0:n_time_steps,-2:-1),  &
                  rhscpre(0:n_time_steps,-2:-1),  &
                  rhscpim(0:n_time_steps,-2:-1),  &
                  stat=istat)
      elseif(izrule.eq.3)then
         allocate(z_mesh(-6:n_z_steps),          &
                  rhsprre(0:n_time_steps,-5:0),  &
                  rhsprim(0:n_time_steps,-5:0),  &
                  rhscpre(0:n_time_steps,-5:0),  &
                  rhscpim(0:n_time_steps,-5:0),  &
                  stat=istat)
      else
         print*,'mbe_propagate_2: Illegal value of izrule. ',izrule
         call mbe_error
      endif
      if(istat.gt.0)then
         print*,'mbe_propagate_2: Error at the allocation of the arrays (3).'
         call mbe_error
      endif
      allocate(sumrho_pre(0:n_time_steps,0:ivmax),  &
               sumrho_pim(0:n_time_steps,0:ivmax),  &
               sumrho_cre(0:n_time_steps,0:ivmax),  &
               sumrho_cim(0:n_time_steps,0:ivmax),  &
               stat=istat)
      if(istat.gt.0)then
         print*,'mbe_propagate_2: Error at the allocation of the arrays (4).'
         call mbe_error
      endif
!
!  Allocation of storage used in case of calculations using the
!  rate equation approximation for the initial time. The arrays arhs and
!  ipiv_r are global within this module.
!
      if(istart.eq.1)then
         if(allocated(arhs))deallocate(arhs)
         if(allocated(ipiv_r))deallocate(ipiv_r)
         allocate(arhs(neqs,neqs),ipiv_r(neqs),stat=istat)
         if(istat.gt.0)then
            print*,'mbe_propagate_2: Error at the allocation of arhs or ipiv_r.'
            call mbe_error
         endif
      endif
!
!  Space mesh. The unit of length is set to 1 mum for propagation.
!
      unit_l = 1.0e-6_kd
!
!  We start the loop over space at j = -1 or j = -6 to accommodate the
!  extra initial steps required by either the mid-point formula (j = -1)
!  or the 4th-order Runge-Kutta formula (j=-6) used between z = 0 and
!  either z = z_step or z = 2z_step.
!
      zmin = 0.0_kd
      z_step = (zmax - zmin)/real(n_z_steps,kd)
      if(izrule.eq.0 .or. izrule.eq.2)then
         jmin = -1
         z_mesh(-1) = zmin
         z_mesh(0) = zmin + z_step/2.0_kd
         do j = 1,n_z_steps
            z_mesh(j) = zmin+real(j,kd)*z_step
         enddo
      else
         if(n_z_steps.lt.2)then
            print*,'Bad value of n_z_steps in mbe_propagate_2'
            call mbe_error
         endif
         jmin = -6
         z_mesh(-6) = zmin
         z_mesh(-5) = zmin + z_step/2.0_kd
         z_mesh(-4) = zmin + z_step/2.0_kd
         z_mesh(-3) = zmin + z_step
         z_mesh(-2) = zmin + z_step
         z_mesh(-1) = zmin + z_step + z_step/2.0_kd
         z_mesh( 0) = zmin + z_step + z_step/2.0_kd
         z_mesh( 1) = zmin + z_step + z_step
         do j = 2,n_z_steps
            z_mesh(j) = zmin+real(j,kd)*z_step
         enddo
      endif
!
      if(izrule.eq.0)then
!
!  Initialization of the calculation without propagation: the absorption
!  coefficient and refractive index at the very front of the medium
!  are set to their vacuum values.
!
         alpha_p = 0.0_kd
         alpha_c = 0.0_kd
         refr_index_p = 1.0_kd
         refr_index_c = 1.0_kd
!
      else
!
!  Propagation parameters (the density is meant to be specified in
!  number of atoms / m^3).
!  
      hzennglp = z_step * (akpr * density / epsilon0) * unit_l
      hzennglc = z_step * (akcp * density / epsilon0) * unit_l
!
!  Minus sign to compensate for an erroneous minus sign in the
!  internal subprograms mbe_calc_Fd_2 and mbe_calc_Fd_3.
      hzennglp = -hzennglp
      hzennglc = -hzennglc
!
      endif
!
!
!  Start the loop over the space steps:
!
      do j = jmin,n_z_steps
!
!  Get the complex field amplitudes at that point in the medium. The
!  algorithm will depend on the value of izrule and whether we do a real
!  propagation calculation or only use Beer's Law (the latter is used when
!  izrule.eq.0).
!
         if(izrule.eq.0)then
!
!  This part of the code needs revisiting...
            print*,'mbe_propagate_2: calculations with izrule.eq.0 might be'
            print*,'incorrect in regards to the variation of the phase'
            print*,'of the field(s), as only the modulus of the amplitude'
            print*,'is changed, which amounts to positing that Re chi = 0.'
            print*,'Also, the calculation of the complex susceptibility might'
            print*,'need correction and re-alignment on obe_susceptiblity.'
            print*,'This part of the code is thus unsafe pending re-analysis.'
            call mbe_error
!  Beer's Law (no change for j.eq.jmin, which is not really a propagation step)
            if(j.gt.jmin)then
               do k = 0,n_time_steps
                  atten_p = exp(-alpha_p(k) * refr_index_p(k) *        &
                                (z_mesh(j)-z_mesh(j-1))*unit_l/2.0_kd)
                  Field_p_re_vec(k) = Field_p_re_vec(k) * atten_p
                  Field_p_im_vec(k) = Field_p_im_vec(k) * atten_p
                  atten_c = exp(-alpha_c(k) * refr_index_c(k) *        &
                                (z_mesh(j)-z_mesh(j-1))*unit_l/2.0_kd)
                  Field_c_re_vec(k) = Field_c_re_vec(k) * atten_c
                  Field_c_im_vec(k) = Field_c_im_vec(k) * atten_c
               enddo
            endif
         elseif(izrule.eq.2)then
            call mbe_calc_Fd_2
         else
            call mbe_calc_Fd_3
         endif
!
!  Prepare the spline interpolation used in the calculation of the density
!  matrix:
         call mbe_spline_prep(Field_p_re_vec,Fd_pre2)
         call mbe_spline_prep(Field_p_im_vec,Fd_pim2)
         call mbe_spline_prep(Field_c_re_vec,Fd_cre2)
         call mbe_spline_prep(Field_c_im_vec,Fd_cim2)
!
!  Loop over velocity (runs only if the calculation involves a Doppler average)
!
         do iv = 1,ivmax
            if(fl_Doppler)then
               current_vel = vmesh(iv)
            else
               current_vel = 0.0_kd ! The value of current_vel is relevant 
                                    ! only when fl_Doppler is .true.
            endif
!
!  Propagate the density matrix in time at the current point in the medium.
!  The initial condition is governed by the variable istart, which is passed to
!  mbe_propagate_2 and must be set to 0, 1 or 2 by the calling program. If istart
!  is 2, the atoms are assumed to be initially in the steady-state they evolve
!  into if starting in the state with a non-zero initial population. Otherwise,
!  the atoms are assumed to be initially in the mixed state defined by the
!  population array popinit: if istart is 0, the coherences are zero, whereas
!  if istart is 1 the coherences are calculated in the rate equations 
!  approximation.
!
            do k = 0,n_time_steps
!
               if(k.eq.0)then
                  rhovec(:) = 0.0_kd
                  do i = 1,nst
                     call ldbl_pop_index(i,m)
                     rhovec(m) = real(popinit(i),kd)
                  enddo
                  Field_p_re = Field_p_re_vec(0)
                  Field_p_im = Field_p_im_vec(0)
                  Field_c_re = Field_c_re_vec(0)
                  Field_c_im = Field_c_im_vec(0)
                  call mbe_set_rhsmat2(rhs_mat,0)
! Unless istart is 1 or 2, the coherences are initially zero.
                  if(istart.eq.2)then
! Steady state
                     call ldbl_steadystate(rhovec,2,rhs_mat)
                  elseif(istart.eq.1)then
! Rate approximation
                     call mbe_calc_coher(rhs_mat,rhovec,1)
                     call mbe_calc_coher(rhs_mat,rhovec,0)
                  endif
               else
!
                  ktime = k - 1
                  call mbe_intobe(imethod,t_mesh(k-1),t_mesh(k),   &
                                               nsbsteps(k),rhovec)
!
               endif
!
!  Sum the relevant coherences (for the propagation calculation)
               sumrho_pre(k,iv) = 0.0_kd
               sumrho_pim(k,iv) = 0.0_kd
               do i = 1,npr
                  sumrho_pre(k,iv) = sumrho_pre(k,iv) - factpr(i)* &
                                               rhovec(idxprre(i))
                  sumrho_pim(k,iv) = sumrho_pim(k,iv) + factpr(i)* &
                                               rhovec(idxprim(i))
               enddo
               sumrho_cre(k,iv) = 0.0_kd
               sumrho_cim(k,iv) = 0.0_kd
               do i = 1,ncp
                  sumrho_cre(k,iv) = sumrho_cre(k,iv) - factcp(i)* &
                                               rhovec(idxcpre(i))
                  sumrho_cim(k,iv) = sumrho_cim(k,iv) + factcp(i)* &
                                               rhovec(idxcpim(i))
               enddo
!
!  Write out the density matrix and the fields (except for the intermediate
!  step between zmin and zmin + z_step or zmin + 2 z_step required by the
!  rule used to start the spatial integration):
               if(      j.eq.jmin                                            &
                  .or. (izrule.eq.2 .and. j.eq. 1 .and. nz_writeout.eq.1)    &
                  .or. (izrule.eq.3 .and. j.eq.-2 .and. nz_writeout.eq.1)    &
                  .or. (j.gt.1 .and. mod(j,nz_writeout).eq.0) )then
                  if(mod(k,nt_writeout).eq.0)then
                     call output_pr(j,k,iv,z_mesh(j),t_mesh(k),rhovec,       &
                                    Field_p_re_vec(k),Field_p_im_vec(k),     &
                                    Field_c_re_vec(k),Field_c_im_vec(k))    
                  endif
               endif
!
            enddo
!
         enddo
!
         call mbe_Doppler
!
!  If "no propagation", calculate the absorption coefficient at this
!  point in the medium and skip the rest. Otherwise, save the
!  sum / integrated coherences for use in the next space step.
!
         if(izrule.eq.0 .and. j.lt.n_z_steps)then
            do k = 0,n_time_steps
               call mbe_calc_susc(Field_p_re_vec(k),Field_p_im_vec(k),       &
                                  sumrho_pre(k,0),sumrho_pim(k,0),           &
                                  akpr,density,refr_index_p(k),alpha_p(k))
               call mbe_calc_susc(Field_c_re_vec(k),Field_c_im_vec(k),    &
                                     sumrho_cre(k,0),sumrho_cim(k,0),        &
                                     akcp,density,refr_index_c(k),alpha_c(k))
            enddo
         elseif(izrule.eq.2)then
            if(j.ne.-1 .and. j.ne.1)then
               rhsprre(:,-2) = rhsprre(:,-1)
               rhsprim(:,-2) = rhsprim(:,-1)
               rhscpre(:,-2) = rhscpre(:,-1)
               rhscpim(:,-2) = rhscpim(:,-1)
            endif
            rhsprre(:,-1) = hzennglp*sumrho_pre(:,0)
            rhsprim(:,-1) = hzennglp*sumrho_pim(:,0)
            rhscpre(:,-1) = hzennglc*sumrho_cre(:,0)
            rhscpim(:,-1) = hzennglc*sumrho_cim(:,0)
         elseif(izrule.eq.3)then
            if(j.eq.-6)then
               rhsprre(:,-5) = hzennglp*sumrho_pre(:,0)
               rhsprim(:,-5) = hzennglp*sumrho_pim(:,0)
               rhscpre(:,-5) = hzennglc*sumrho_cre(:,0)
               rhscpim(:,-5) = hzennglc*sumrho_cim(:,0)
               rhsprre(:,-4) = rhsprre(:,-5)
               rhsprim(:,-4) = rhsprim(:,-5) 
               rhscpre(:,-4) = rhscpre(:,-5)
               rhscpim(:,-4) = rhscpim(:,-5)
            elseif(j.eq.-5 .or. j.eq.-1)then
               rhsprre(:,-3) = hzennglp*sumrho_pre(:,0)
               rhsprim(:,-3) = hzennglp*sumrho_pim(:,0)
               rhscpre(:,-3) = hzennglc*sumrho_cre(:,0)
               rhscpim(:,-3) = hzennglc*sumrho_cim(:,0)
            elseif(j.eq.-4 .or. j.eq.0)then
               rhsprre(:,-2) = hzennglp*sumrho_pre(:,0)
               rhsprim(:,-2) = hzennglp*sumrho_pim(:,0)
               rhscpre(:,-2) = hzennglc*sumrho_cre(:,0)
               rhscpim(:,-2) = hzennglc*sumrho_cim(:,0)
            elseif(j.eq.-3 .or. j.eq.1)then
               rhsprre(:,-1) = hzennglp*sumrho_pre(:,0)
               rhsprim(:,-1) = hzennglp*sumrho_pim(:,0)
               rhscpre(:,-1) = hzennglc*sumrho_cre(:,0)
               rhscpim(:,-1) = hzennglc*sumrho_cim(:,0)
            elseif(j.eq.-2)then
               rhsprre(:,-4) = hzennglp*sumrho_pre(:,0)
               rhsprim(:,-4) = hzennglp*sumrho_pim(:,0)
               rhscpre(:,-4) = hzennglc*sumrho_cre(:,0)
               rhscpim(:,-4) = hzennglc*sumrho_cim(:,0)
            elseif(j.eq.2)then
               rhsprre(:,-3) = rhsprre(:,-5)
               rhsprim(:,-3) = rhsprim(:,-5) 
               rhscpre(:,-3) = rhscpre(:,-5)
               rhscpim(:,-3) = rhscpim(:,-5)
               rhsprre(:,-2) = rhsprre(:,-4)
               rhsprim(:,-2) = rhsprim(:,-4) 
               rhscpre(:,-2) = rhscpre(:,-4)
               rhscpim(:,-2) = rhscpim(:,-4)
               rhsprre(:,-1) = hzennglp*sumrho_pre(:,0)
               rhsprim(:,-1) = hzennglp*sumrho_pim(:,0)
               rhscpre(:,-1) = hzennglc*sumrho_cre(:,0)
               rhscpim(:,-1) = hzennglc*sumrho_cim(:,0)
            elseif(j.gt.2)then
               rhsprre(:,-3) = rhsprre(:,-2)
               rhsprim(:,-3) = rhsprim(:,-2)
               rhscpre(:,-3) = rhscpre(:,-2)
               rhscpim(:,-3) = rhscpim(:,-2)
               rhsprre(:,-2) = rhsprre(:,-1)
               rhsprim(:,-2) = rhsprim(:,-1)
               rhscpre(:,-2) = rhscpre(:,-1)
               rhscpim(:,-2) = rhscpim(:,-1)
               rhsprre(:,-1) = hzennglp*sumrho_pre(:,0)
               rhsprim(:,-1) = hzennglp*sumrho_pim(:,0)
               rhscpre(:,-1) = hzennglc*sumrho_cre(:,0)
               rhscpim(:,-1) = hzennglc*sumrho_cim(:,0)
            endif
         endif
!
      enddo
!
      deallocate(Field_p_re_vec,Field_p_im_vec,Field_c_re_vec,Field_c_im_vec)
      deallocate(Fd_p_re_prev,Fd_p_im_prev,Fd_c_re_prev,Fd_c_im_prev,         &
         Fd_pre2,Fd_pim2,Fd_cre2,Fd_cim2,rhs_mat,z_mesh,                      &
         sumrho_pre,sumrho_pim,sumrho_cre,sumrho_cim)
      if(izrule.eq.0)then
         deallocate(alpha_p,alpha_c,refr_index_p,refr_index_c)
      else
         deallocate(rhsprre,rhsprim,rhscpre,rhscpim)
      endif
!
      contains
 
!!!!!!!!!!!!!!!!!!!
 
      subroutine mbe_calc_Fd_2
!
!  Warning: the right-hand sides are multiplied by an erroneous -1 factor.
!  This is compensated by a change of sign of hzennglp and hzennglc.
!
!  Calculate the new field amplitudes. Use the mid point formula for the 
!  first step and the 2nd order Adams-Bashforth for the rest.
!
!  Entrance of the medium (j = jmin): The fields are the applied fields.
!                                     mbe_propagate_2 expects to find the
!                                     corresponding complex fields 
!                                     already in Field_p_re_vec etc.
!  First step into the medium: Use the mid-point formula to calculate
!                              the fields.
!  Subsequent steps: Use the Adams-Bashforth method.
!
      if(j.eq.-1)then
!  Use the existing contents of Field_p_re_vec etc. Save these values
!  for later use in the midpoint formula.
         Fd_p_re_prev = Field_p_re_vec
         Fd_p_im_prev = Field_p_im_vec
         Fd_c_re_prev = Field_c_re_vec
         Fd_c_im_prev = Field_c_im_vec
      elseif(j.eq.0)then
!  Midpoint formula (first half step)
         Field_p_re_vec = Fd_p_re_prev - rhsprim(:,-1)*0.5_kd
         Field_p_im_vec = Fd_p_im_prev + rhsprre(:,-1)*0.5_kd
         Field_c_re_vec = Fd_c_re_prev - rhscpim(:,-1)*0.5_kd
         Field_c_im_vec = Fd_c_im_prev + rhscpre(:,-1)*0.5_kd
      elseif(j.eq.1)then
!  Midpoint formula (second half step)
         Field_p_re_vec = Fd_p_re_prev - rhsprim(:,-1)
         Field_p_im_vec = Fd_p_im_prev + rhsprre(:,-1)
         Field_c_re_vec = Fd_c_re_prev - rhscpim(:,-1)
         Field_c_im_vec = Fd_c_im_prev + rhscpre(:,-1)
      elseif(j.ge.2)then
!  Adams-Bashforth step
         Field_p_re_vec = Field_p_re_vec - rhsprim(:,-1)*1.5_kd &
                                         + rhsprim(:,-2)*0.5_kd
         Field_p_im_vec = Field_p_im_vec + rhsprre(:,-1)*1.5_kd &
                                         - rhsprre(:,-2)*0.5_kd
         Field_c_re_vec = Field_c_re_vec - rhscpim(:,-1)*1.5_kd &
                                         + rhscpim(:,-2)*0.5_kd
         Field_c_im_vec = Field_c_im_vec + rhscpre(:,-1)*1.5_kd &
                                         - rhscpre(:,-2)*0.5_kd
      else
!
         print*,'Logic error in mbe_propagate_2 or mbe_calc_Fd_2 ',j
         call mbe_error
!
      endif
!
      end subroutine mbe_calc_Fd_2

!!!!!!!!!!!!!!!!!!!
 
      subroutine mbe_calc_Fd_3
!
!  Warning: the right-hand sides are multiplied by an erroneous -1 factor.
!  This is compensated by a change of sign of hzennglp and hzennglc.
!
!  Calculate the new field amplitudes. Use RK4 for the first two steps and
!  the 3rd order Adams-Bashforth for the rest.
!
!  Entrance of the medium (j = jmin): The fields are the applied fields.
!                                     mbe_propagate_2 expects to find the
!                                     corresponding complex fields 
!                                     already in Field_p_re_vec etc.
!  First two steps into the medium: Use RK4 to calculate the fields.
!  Subsequent steps: Use the Adams-Bashforth method.
!
      if(j.eq.-6)then
!
         continue  ! Use the existing contents of Field_p_re_vec etc
!
      elseif(j.eq.-5 .or. j.eq.-1)then
!
         Fd_p_re_prev = Field_p_re_vec
         Fd_p_im_prev = Field_p_im_vec
         Fd_c_re_prev = Field_c_re_vec
         Fd_c_im_prev = Field_c_im_vec
!  RK4 (first step)
         Field_p_re_vec = Fd_p_re_prev - rhsprim(:,-4)
         Field_p_im_vec = Fd_p_im_prev + rhsprre(:,-4)
         Field_c_re_vec = Fd_c_re_prev - rhscpim(:,-4)
         Field_c_im_vec = Fd_c_im_prev + rhscpre(:,-4)
      elseif(j.eq.-4 .or. j.eq.0)then
!  RK4 (second step)
         Field_p_re_vec = Fd_p_re_prev - rhsprim(:,-3)*0.5_kd 
         Field_p_im_vec = Fd_p_im_prev + rhsprre(:,-3)*0.5_kd 
         Field_c_re_vec = Fd_c_re_prev - rhscpim(:,-3)*0.5_kd 
         Field_c_im_vec = Fd_c_im_prev + rhscpre(:,-3)*0.5_kd 
      elseif(j.eq.-3 .or. j.eq.1)then
!  RK4 (third step)
         Field_p_re_vec = Fd_p_re_prev - rhsprim(:,-2)*0.5_kd 
         Field_p_im_vec = Fd_p_im_prev + rhsprre(:,-2)*0.5_kd 
         Field_c_re_vec = Fd_c_re_prev - rhscpim(:,-2)*0.5_kd 
         Field_c_im_vec = Fd_c_im_prev + rhscpre(:,-2)*0.5_kd 
      elseif(j.eq.-2 .or. j.eq.2)then
!  RK4 (fourth step)
         Field_p_re_vec = Fd_p_re_prev - (       rhsprim(:,-4)+        &
                                          2.0_kd*rhsprim(:,-3)+        &
                                          2.0_kd*rhsprim(:,-2)+        &
                                                 rhsprim(:,-1)  )/6.0_kd
         Field_p_im_vec = Fd_p_im_prev + (       rhsprre(:,-4)+        &
                                          2.0_kd*rhsprre(:,-3)+        &
                                          2.0_kd*rhsprre(:,-2)+        &
                                                 rhsprre(:,-1)  )/6.0_kd
         Field_c_re_vec = Fd_c_re_prev - (       rhscpim(:,-4)+        &
                                          2.0_kd*rhscpim(:,-3)+        &
                                          2.0_kd*rhscpim(:,-2)+        &
                                                 rhscpim(:,-1)  )/6.0_kd
         Field_c_im_vec = Fd_c_im_prev + (       rhscpre(:,-4)+        &
                                          2.0_kd*rhscpre(:,-3)+        &
                                          2.0_kd*rhscpre(:,-2)+        &
                                                 rhscpre(:,-1)  )/6.0_kd
!
      elseif(j.ge.3)then
!
         Fd_p_re_prev = Field_p_re_vec
         Fd_p_im_prev = Field_p_im_vec
         Fd_c_re_prev = Field_c_re_vec
         Fd_c_im_prev = Field_c_im_vec
!
!  Adams-Bashforth step
!  3rd order
         Field_p_re_vec = Fd_p_re_prev - (23.0_kd*rhsprim(:,-1) - &
                                          16.0_kd*rhsprim(:,-2) + &
                                           5.0_kd*rhsprim(:,-3))/12.0_kd 
         Field_p_im_vec = Fd_p_im_prev + (23.0_kd*rhsprre(:,-1) - &
                                          16.0_kd*rhsprre(:,-2) + &
                                           5.0_kd*rhsprre(:,-3))/12.0_kd 
         Field_c_re_vec = Fd_c_re_prev - (23.0_kd*rhscpim(:,-1) - &
                                          16.0_kd*rhscpim(:,-2) + &
                                           5.0_kd*rhscpim(:,-3))/12.0_kd 
         Field_c_im_vec = Fd_c_im_prev + (23.0_kd*rhscpre(:,-1) - &
                                          16.0_kd*rhscpre(:,-2) + &
                                           5.0_kd*rhscpre(:,-3))/12.0_kd 
!  2nd order
!        Field_p_re_vec = Fd_p_re_prev - rhsprim(:,-1)*1.50_kd &
!                                      + rhsprim(:,-2)*0.50_kd
!        Field_p_im_vec = Fd_p_im_prev + rhsprre(:,-1)*1.50_kd &
!                                      - rhsprre(:,-2)*0.50_kd
!        Field_c_re_vec = Fd_c_re_prev - rhscpim(:,-1)*1.50_kd &
!                                      + rhscpim(:,-2)*0.50_kd
!        Field_c_im_vec = Fd_c_im_prev + rhscpre(:,-1)*1.50_kd &
!                                      - rhscpre(:,-2)*0.50_kd
!  1st order
!        Field_p_re_vec = Fd_p_re_prev - rhsprim(:,-1)
!        Field_p_im_vec = Fd_p_im_prev + rhsprre(:,-1)
!        Field_c_re_vec = Fd_c_re_prev - rhscpim(:,-1)
!        Field_c_im_vec = Fd_c_im_prev + rhscpre(:,-1)

!  Calculate the density matrix with the new field amplitudes:
!
!  Prepare the spline interpolation used in the calculation of the density
!  matrix:
         call mbe_spline_prep(Field_p_re_vec,Fd_pre2)
         call mbe_spline_prep(Field_p_im_vec,Fd_pim2)
         call mbe_spline_prep(Field_c_re_vec,Fd_cre2)
         call mbe_spline_prep(Field_c_im_vec,Fd_cim2)
!  Loop over velocity 
         do iv = 1,ivmax
            if(fl_Doppler)then
               current_vel = vmesh(iv)
            else
               current_vel = 0.0_kd
            endif
!  Propagate the density matrix in time at the current point in the medium.
            do k = 0,n_time_steps
               if(k.eq.0)then
                  rhovec(:) = 0.0_kd
                  do i = 1,nst
                     call ldbl_pop_index(i,m)
                     rhovec(m) = real(popinit(i),kd)
                  enddo
                  Field_p_re = Field_p_re_vec(0)
                  Field_p_im = Field_p_im_vec(0)
                  Field_c_re = Field_c_re_vec(0)
                  Field_c_im = Field_c_im_vec(0)
                  call mbe_set_rhsmat2(rhs_mat,0)
! Unless istart is 1 or 2, the coherences are initially zero.
                  if(istart.eq.2)then
! Steady state
                     call ldbl_steadystate(rhovec,2,rhs_mat)
                  elseif(istart.eq.1)then
! Rate approximation
                     call mbe_calc_coher(rhs_mat,rhovec,1)
                     call mbe_calc_coher(rhs_mat,rhovec,0)
                  endif
               else
                  ktime = k - 1
                  call mbe_intobe(imethod,t_mesh(k-1),t_mesh(k),   &
                                               nsbsteps(k),rhovec)
               endif
!  Sum the relevant coherences (for the propagation calculation)
               sumrho_pre(k,iv) = 0.0_kd
               sumrho_pim(k,iv) = 0.0_kd
               do i = 1,npr
                  sumrho_pre(k,iv) = sumrho_pre(k,iv) - factpr(i)* &
                                                rhovec(idxprre(i))
                  sumrho_pim(k,iv) = sumrho_pim(k,iv) + factpr(i)* &
                                                rhovec(idxprim(i))
               enddo
               sumrho_cre(k,iv) = 0.0_kd
               sumrho_cim(k,iv) = 0.0_kd
               do i = 1,ncp
                  sumrho_cre(k,iv) = sumrho_cre(k,iv) - factcp(i)* &
                                                rhovec(idxcpre(i))
                  sumrho_cim(k,iv) = sumrho_cim(k,iv) + factcp(i)* &
                                                rhovec(idxcpim(i))
               enddo
            enddo
         enddo
         call mbe_Doppler
         rhsprre(:,0) = hzennglp*sumrho_pre(:,0)
         rhsprim(:,0) = hzennglp*sumrho_pim(:,0)
         rhscpre(:,0) = hzennglc*sumrho_cre(:,0)
         rhscpim(:,0) = hzennglc*sumrho_cim(:,0)
!
!  Adams-Moulton step
!
!  4th order
         Field_p_re_vec = Fd_p_re_prev - ( 9.0_kd*rhsprim(:,0)  + &
                                         19.0_kd*rhsprim(:,-1) - &
                                          5.0_kd*rhsprim(:,-2) + &
                                               rhsprim(:,-3) )/24.0_kd !&
         Field_p_im_vec = Fd_p_im_prev + ( 9.0_kd*rhsprre(:,0)  + &
                                         19.0_kd*rhsprre(:,-1) - &
                                          5.0_kd*rhsprre(:,-2) + &
                                               rhsprre(:,-3) )/24.0_kd !& 
         Field_c_re_vec = Fd_c_re_prev - ( 9.0_kd*rhscpim(:,0)  + &
                                         19.0_kd*rhscpim(:,-1) - &
                                          5.0_kd*rhscpim(:,-2) + &
                                               rhscpim(:,-3) )/24.0_kd !&
         Field_c_im_vec = Fd_c_im_prev + ( 9.0_kd*rhscpre(:,0)  + &
                                         19.0_kd*rhscpre(:,-1) - &
                                          5.0_kd*rhscpre(:,-2) + &
                                               rhscpre(:,-3) )/24.0_kd !&
!  3rd order
!        Field_p_re_vec = Fd_p_re_prev - ( 5.0_kd*rhsprim(:,0)  + &
!                                         8.0_kd*rhsprim(:,-1) - &
!                                              rhsprim(:,-2) )/12.0_kd  
!        Field_p_im_vec = Fd_p_im_prev + ( 5.0_kd*rhsprre(:,0)  + &
!                                         8.0_kd*rhsprre(:,-1) - &
!                                              rhsprre(:,-2) )/12.0_kd  
!        Field_c_re_vec = Fd_c_re_prev - ( 5.0_kd*rhscpim(:,0)  + &
!                                         8.0_kd*rhscpim(:,-1) - &
!                                              rhscpim(:,-2) )/12.0_kd  
!        Field_c_im_vec = Fd_c_im_prev + ( 5.0_kd*rhscpre(:,0)  + &
!                                         8.0_kd*rhscpre(:,-1) - &
!                                              rhscpre(:,-2) )/12.0_kd  
!
      else
!
         print*,'Logic error in mbe_propagate_2 or mbe_calc_Fd_4 ',j
         call mbe_error
!
      endif
!
      end subroutine mbe_calc_Fd_3

!!!!!!!!!!!!!!!!!!!

      subroutine mbe_Doppler
!
      implicit none
      real(kd), dimension(:,:), allocatable :: sum_mod,sum_arg
      real(kd), dimension(:), allocatable :: sum_mod_v,sum_arg_v
      real(kd), dimension(:), allocatable :: s_pre,s_pim,s_cre,s_cim
      real(kd), dimension(:), allocatable :: vtest
      real(kd) :: v
      real(kd) :: tot1,tot2,ynt
      integer :: istat,j,numb
!
      if(.not.fl_Doppler)then
         sumrho_pre(:,0) = sumrho_pre(:,1)
         sumrho_pim(:,0) = sumrho_pim(:,1)
         sumrho_cre(:,0) = sumrho_cre(:,1)
         sumrho_cim(:,0) = sumrho_cim(:,1)
         return
      endif
!
!  Calculation with Doppler averaging.
!
      sumrho_pre(:,0) = 0.0_kd
      sumrho_pim(:,0) = 0.0_kd
      sumrho_cre(:,0) = 0.0_kd
      sumrho_cim(:,0) = 0.0_kd
      do i = 1,ivmax
         sumrho_pre(:,0) = sumrho_pre(:,0)+sumrho_pre(:,i)*fMvweight(i)
         sumrho_pim(:,0) = sumrho_pim(:,0)+sumrho_pim(:,i)*fMvweight(i)
         sumrho_cre(:,0) = sumrho_cre(:,0)+sumrho_cre(:,i)*fMvweight(i)
         sumrho_cim(:,0) = sumrho_cim(:,0)+sumrho_cim(:,i)*fMvweight(i)
      enddo
!
      end subroutine mbe_Doppler

!!!!!!!!!!!!!!!!!!!
 
      end subroutine mbe_propagate_2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine mbe_calc_coher(rhs_mat,rho_vec,iinit)
!
!  Essentially the same as ldbl_calc_coher.
!
!  Given the populations, calculate the coherences in the rate equation
!  approximation. The populations are passed to mbe_calc_coher through
!  the respective elements of rho_vec, and the populations and coherences
!  are passed back to the calling program through that array. The initial 
!  values of the elements of rho_vec corresponding to the coherences are
!  unused and overwritten.
!
!  This subroutine makes use of the arrays ipiv_r and arhs which it
!  inherits from the calling program (these arrays are global within the
!  module) and of the array rhs_mat. The latter is passed to mbe_calc_coher
!  through the list of arguments. This array is allocatable. It is meant to
!  have been properly allocated prior the call to mbe_calc_coher; no check
!  is made on whether this happened.
!
!  The calculation is simply prepared (i.e., the matrix representing the system
!  is factorized) if iinit.eq.1 at entry. In this case, rho_vec is not used
!  at all. Clearly, mbe_calc_coher must be called with iinit set to 1 prior
!  coherences are calculated (which happens when iinit.ne.1 at entry) and/or
!  each time rhs_mat is modified.
!  
      implicit none
      integer, parameter :: neqs = nst*nst
      real(kd), dimension (:,:), intent(in) :: rhs_mat
      real(kd), dimension(neqs), intent(inout) :: rho_vec
      integer :: iinit
      real(kd), dimension(nst) :: popul
      integer :: i,ifail,info,lda,ldb,m,numb,nrhs
      logical, save :: fl_init = .true.
      external dgesv,sgesv
!
!  arhs is calculated and factorized if iinit.eq.1; otherwise rho_vec is
!  calculated.
!
      if(iinit.eq.1)then
!
!  The matrix representing the system is arhs.
!
         arhs = rhs_mat
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
            if(all(rhs_mat(i,:) .eq. 0.0_kd))then
               arhs(i,i) = 1.0_kd     
            endif
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
!  Factorizes arhs. This is done by the subroutine s/dgetrf of lapack.
!    numb: the number of equations in the linear system to be solved
!    lda: leading dimension of the array arhs
!    info: an integer variable whose value is set to zero by s/dgesv in case
!       of successful completion or to another value in case of problem
!       (see online information about s/dgesv for details)
!
         numb = neqs
         lda = neqs
!
         if(kd.eq.kind(1.0))then
            call sgetrf(numb,numb,arhs,lda,ipiv_r,info)
         elseif(kd.eq.kind(1.d0))then
            call dgetrf(numb,numb,arhs,lda,ipiv_r,info)
         else
            print*,'Logic error detected in mbe_calc_coher'
            call mbe_error
         endif
!
!  info.ne.0 : error condition
!  info.eq.0 : normal return from s/dgetrf
!
         if(info.ne.0)then
            print*,'mbe_calc_coher: s/dgetrf returns with info = ',info
            call mbe_error
         endif
!
!  That's all for an initialization call. Reset fl_init.
         fl_init = .false.
!
      else
!
!  Calculate the coherences.
!
!  First check that the value fl_init.eqv.false
!
         if(fl_init)then
            print*,'mbe_calc_coher was called with iinit.ne.1 before being'
            print*,'first called with iinit.eq.1...'
            call mbe_error
         endif
!
!  Set all the elements of rho_vec to zero, except the populations.
         do i = 1,nst
            call ldbl_pop_index(i,m)
            popul(i) = rho_vec(m)
         enddo
         rho_vec = 0.0_kd
         do i = 1,nst
            call ldbl_pop_index(i,m)
            rho_vec(m) = popul(i)
         enddo
!
!  Calculate the rhs of the system of linear equations.
!
         rho_vec = -matmul(rhs_mat,rho_vec)
!
!  Solve the system of linear equations with rho_vec as their rhs.
!  Variables passed to s/dgetrs (solution of system of linear equations,
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
!  Programming note: numb and lda are set here, too, as they are not saved.
!
         numb = neqs
         nrhs = 1
         lda = neqs
         ldb = neqs
!
         if(kd.eq.kind(1.0))then
            call sgetrs('N',numb,nrhs,arhs,lda,ipiv_r,rho_vec,ldb,info)
         elseif(kd.eq.kind(1.d0))then
            call dgetrs('N',numb,nrhs,arhs,lda,ipiv_r,rho_vec,ldb,info)
         else
            print*,'Logic error detected in mbe_calc_coher'
            call mbe_error
         endif
!
!  info.ne.0 : error condition
!  info.eq.0 : normal return from s/dgetrs.
!
         if(info.ne.0)then
            print*,'mbe_calc_coher: s/dgetrs returns with info = ',info
            call mbe_error
         endif
!
!  At this stage, rho_vec contains the correct coherences (within the rate
!  equation approximation) but not the correct populations. Replace the
!  (wrong) populations by those passed to mbe_calc_coher (they were saved
!  in the array popul):
            do i = 1,nst
               call ldbl_pop_index(i,m)
               rho_vec(m) = popul(i)
            enddo
!
!  That's all...
!
      endif
!
      end subroutine mbe_calc_coher
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
      subroutine mbe_set_Doppler
!
!  Set the global variables urms, n_v_values, vmesh, fMvweight and
!  fl_Doppler_notinit through a call to obe_get_Dopplerpar. An error
!  condition is flagged if this call has not been preceded by a call to
!  obe_set_Doppler.  A successful call to obe_set_Doppler has the
!  effect of (re)-allocating the arrays vmesh and fMvweight, and
!  overwriting all these global variables with new information.
!
!  urms is the r.m.s. speed in the light propagation direction, in m/s.
!  n_v_values is the number of abscissas in the integration over velocity
!  classes, vmesh are the abscissas and fMvweight are the corresponding
!  weights multiplied by the Maxwell distribution.
!
      implicit none
!
      call obe_get_Dopplerpar(fl_Doppler_notinit,urms,n_v_values, &
                              vmesh,fMvweight)
!
      if(fl_Doppler_notinit)then
         print*,'mbe_set_Doppler: The parameters of the integration'
         print*,'over the velocity classes should have been set'
         print*,'by a prior call to obe_set_Doppler.'
         call mbe_error
      endif
      if(.not.allocated(vmesh) .or. .not.allocated(fMvweight))then
         print*,'mbe_set_Doppler: Logic error...'
         print*,allocated(vmesh),allocated(fMvweight)
         call mbe_error
      endif
      if(size(vmesh).ne.n_v_values .or. size(fMvweight).ne.n_v_values)then
         print*,'mbe_set_Doppler: Inconsistent arrays sizes.'
         print*,n_v_values,size(vmesh),size(fMvweight)
         call mbe_error
      endif
!
      end subroutine mbe_set_Doppler

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine mbe_set_envlp(field_name,pulse_type,tw,t0,t1,force0)
!
!  This subroutine is used to set the global variables pulse_type_p,
!  tw_p, t0_p, t1_p, force0_p, pulse_type_c, tw_c, t0_c, t1_c and
!  force0_c defining the envelopes of the probe and coupling fields.
!  See mbe_calc_envelope for the possible types of pulse and the meaning
!  of these parameters.
!
!  The arguments tw, t0, t1 and force0 are optional. tw, t0 and force0 must
!  be specified for Gaussian and sech pulses, t0 for square-step pulses,
!  and t0 and t1 for flat-top pulses. None of tw, t0, t1 or force0 need
!  to be specified for cw fields. Any unnecessary argument is not used, even if
!  specified.
!
      implicit none
      double precision, intent(in), optional :: tw,t0,t1
      character*5, intent(in) :: field_name
      character*2, intent(in) :: pulse_type
      logical, intent(in), optional :: force0
!
      if(pulse_type.eq.'Gs' .and. (.not.present(tw) .or.   &
            .not.present(t0) .or. .not.present(force0)))then
         print*,'mbe_set_envlp: tw, t0 and force0 should be specified for a'
         print*,'Gaussian pulse.'
         call mbe_error
      endif
      if(pulse_type.eq.'hs' .and. (.not.present(tw) .or.   &
            .not.present(t0) .or. .not.present(force0)))then
         print*,'mbe_set_envlp: tw, t0 and force0 should be specified for a'
         print*,'sech pulse.'
         call mbe_error
      endif
      if(pulse_type.eq.'sq' .and. .not.present(t0))then
         print*,'mbe_set_envlp: t0 should be specified for a'
         print*,'square step pulse.'
      endif
      if(pulse_type.eq.'ft' .and.                        &
               (.not.present(t0) .or. .not.present(t1)))then
         print*,'mbe_set_envlp: t0 and t1 should be specified for a'
         print*,'flat-top pulse.'
      endif
!
      if(field_name.eq.'probe')then
         pulse_type_p = pulse_type
         if(present(tw))tw_p = tw
         if(present(t0))t0_p = t0
         if(present(t1))t1_p = t1
         if(present(force0))force0_p = force0
         fl_p_env_notset = .false.
      elseif(field_name.eq.'coupl')then
         pulse_type_c = pulse_type
         if(present(tw))tw_c = tw
         if(present(t0))t0_c = t0
         if(present(t1))t1_c = t1
         if(present(force0))force0_c = force0
         fl_c_env_notset = .false.
      else
         print*,'mbe_set_envlp: field name not recognized'
         print*,field_name
         call mbe_error
      endif
!
      end subroutine mbe_set_envlp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine mbe_calc_applfld(t,tw,pulse_type,t0,t1,cfield0,cfieldt)
!
!  Calculates an applied field at time t. This applied field could be
!  a Gaussian pulse centered at t = t0 (if pulse_type = 'Gs'), or
!  a sech pulse centered at t = t0 (if pulse_type = 'hs'), or a square
!  step pulse turned on at t = t0 (if pulse_type = 'sq'), or a flat-top
!  pulse smoothly turned on with a cos^2-profile between t0 and t1 (if
!  pulse_type = 'ft'). A cw field is also possible (pulse_type =
!  'cw'). The value of t1 is relevant only for flat-top pulses.
!
      implicit none
      double precision, intent(in) :: t,tw,t0,t1
      character*2, intent(in) :: pulse_type
      double complex, intent(in) :: cfield0
      double complex, intent(out) :: cfieldt
      double precision, parameter :: fourlog2 = 4.d0*dlog(2.d0)
      double precision, parameter :: pi = 4.d0*datan(1.d0)
      double precision :: ts
!
      if(pulse_type.eq.'cw')then
         cfieldt = (1.0d0,0.0d0)*cfield0
      elseif(pulse_type.eq.'Gs')then
         if(tw.eq.0.d0)then
            print*,'mbe_calc_applfld: tw is zero...'
            call mbe_error
         endif
         cfieldt = (1.0d0,0.0d0)*cfield0*exp(-fourlog2*((t-t0)/tw)**2)
      elseif(pulse_type.eq.'hs')then
         if(tw.eq.0.d0)then
            print*,'mbe_calc_applfld: tw is zero...'
            call mbe_error
         endif
         ts = tw/2.633916d0
         cfieldt = (1.0d0,0.0d0)*cfield0/dcosh((t-t0)/ts)
      elseif(pulse_type.eq.'sq')then
         if(t.lt.t0)then
            cfieldt = (0.0d0,0.0d0)
         else
            cfieldt = (1.0d0,0.0d0)*cfield0
         endif
      elseif(pulse_type.eq.'ft')then
         if(t1.le.t0)then
            print*,'mbe_calc_applfld: t1 must be larger than t0 for a'
            print*,'flat-top pulse smoothly turned on.', t0,t1
            call mbe_error
         else
            if(t.lt.t0)then
               cfieldt = (0.0d0,0.0d0)
            elseif(t.le.t1)then
               cfieldt = (1.0d0,0.0d0)*cfield0*cos(pi*(t1-t)/(2.d0*(t1-t0)))**2
            else
               cfieldt = (1.0d0,0.0d0)*cfield0
            endif
         endif
      else
         print*,'mbe_calc_applfld: Illegal pulse type. ', pulse_type
         call mbe_error
      endif
!
      end subroutine mbe_calc_applfld

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine mbe_set_tdfields_A(tmin,tmax,n_time_steps_in,nflds,  &
                                    iinterp,cfield_p,cfield_c)
!
!  Calculate the temporal envelope of the fields, or prepare their
!  calculation, given the information previously passed to mbe through
!  call or calls to mbe_set_envlp.
!
!  n_time_steps_in is the number of time steps. The total number of
!  values of t at which the density matrix must be obtained by the mbe routines,
!  including t = tmin and t = tmax, is n_time_steps + 1. This last
!  number is also the number of times at which the applied field(s) is (are)
!  tabulated by mbe_set_tdfields_A, if required. The variable nflds
!  should be 1 for single-field calculations or 2 for two-field
!  calculations. The variables cfield_p and (optionally) cfield_c are the
!  complex amplitudes of the respective fields at the maximum of their
!  envelope.
!
!  The variable iinterp should be either 0 or 1. If 1, the applied
!  fields are tabulated by mbe_set_tdfields_A, and the fields are then
!  interpolated as necessary in the course of the time integration. If
!  0, the time mesh is set by mbe_set_tdfields_A and the applied fields are
!  still tabulated, but they are calculated directly rather than by
!  interpolation at intermediate values of t.
!  
      implicit none
      double complex, intent(in) :: cfield_p
      double complex, intent(in), optional :: cfield_c
      double precision, intent(in) :: tmin,tmax
      integer, intent(in) :: iinterp,n_time_steps_in,nflds
      double precision :: time_step
      integer :: istat,k
!
!  Sanity checks...
      if(n_time_steps_in.lt.1 .or. nflds.gt.2 .or. nflds.le.0)then
         print*,'mbe_set_tdfields_A was called with an incorrect value'
         print*,'of n_t_values or nflds... ',n_time_steps_in,nflds
         call mbe_error
      endif
      if(fl_p_env_notset)then
         print*,'The details of the envelope of the probe field should'
         print*,'have been set before the first call to mbe_set_tdfields_A.'
         call mbe_error
      elseif(fl_c_env_notset .and. nflds.eq.2)then
         print*,'The details of the envelope of the coupling field should'
         print*,'have been set before the first call to mbe_set_tdfields_A'
         print*,'for a 2-field calculation.'
         call mbe_error
      endif
      if(nflds.eq.2 .and. .not.present(cfield_c))then
         print*,'mbe_set_tdfields_A: the optional argument cfield_c is'
         print*,'missing although a 2-field calculation is required.'
         call mbe_error
      endif
      if(tmin.ge.tmax)then
         print*,'mbe_set_tdfields_A: tmin is not smaller than tmax.'
         print*,tmin,tmax
         call mbe_error
      endif
      if(iinterp.ne.0 .and. iinterp.ne.1)then
         print*,'mbe_set_tdfields_A: invalid value of iinterp: ',iinterp
         call mbe_error
      endif
!
!  cfield0_p and cfield0_c are global variables within mbe.
      cfield0_p = cfield_p
      if(nflds.eq.2)cfield0_c = cfield_c
!
!  n_time_steps and t_mesh are global variables within mbe.
      n_time_steps = n_time_steps_in
      if(allocated(t_mesh))deallocate(t_mesh)
      allocate(t_mesh(0:n_time_steps), stat = istat)
      if(istat.ne.0)then
         print*,'mbe_set_tdfields_A: error at the allocation of t_mesh.'
         call mbe_error
      endif
!
      time_step = (tmax - tmin)/dble(n_time_steps)
      do k = 0, n_time_steps-1
         t_mesh(k) = tmin + dble(k)*time_step
      enddo
      t_mesh(n_time_steps) = tmax
!
      if(iinterp.eq.0)then
!  The time dependence of the applied field will be directly calculated 
!  later on:
         fl_spline = .false.
      else
!  The time dependence of the applied field will be calculated by
!  interpolation between tabulation points.
         fl_spline = .true.
      endif
!
!  The array Field_p_vec is allocated by mbe_calc_envelope. The same
!  also applies to Field_c_vec, unless nflds.eq.1 in which case this
!  array is allocated below.
      if(allocated(Field_p_vec))deallocate(Field_p_vec)
      if(allocated(Field_c_vec))deallocate(Field_c_vec)
      call mbe_calc_envelope(tw_p,pulse_type_p,t0_p,t1_p,force0_p,  &
                             cfield0_p,Field_p_vec)
      if(nflds.eq.2)then
         call mbe_calc_envelope(tw_c,pulse_type_c,t0_c,t1_c,force0_c, &
                                cfield0_c,Field_c_vec)
      else
         allocate(Field_c_vec(0:n_time_steps), stat = istat)
         if(istat.ne.0)then
            print*,'mbe_set_tdfields_A: error at the allocation '
            print*,'of Field_c_vec.'
            call mbe_error
         endif
         Field_c_vec = (0.d0,0.d0)
      endif
!
      fl_notinit_fields = .false.
!
      contains

!!!!!!!!!!!!!!!!!!!

      subroutine mbe_calc_envelope(tw,pulse_type,t0,t1,force0,  &
                                   cfield0,field_vec)
!
!  Returns a table of values of an applied field calculated at the times
!  specified in the array t_mesh.
!
!  For Gaussian or sech pulses, the field is forced to be zero at
!  t=t_mesh(0) (which is tmin) if force0 .eqv. .true..
!
      implicit none
      double complex, intent(in) :: cfield0
      double precision, intent(in) :: tw,t0,t1
      character*2, intent(in) :: pulse_type
      logical, intent(in) :: force0
      double complex, dimension(:), allocatable, intent(out) :: field_vec
      double precision, parameter :: fourlog2 = 4.d0*dlog(2.d0)
      double precision, parameter :: pi = 4.d0*datan(1.d0)
      double precision :: t,ts
      integer :: istat
!
      if(allocated(field_vec))deallocate(field_vec)
      allocate(field_vec(0:n_time_steps), stat = istat)
      if(istat.ne.0)then
         print*,'mbe_calc_envelope: error at the allocation of the array.'
         call mbe_error
      endif
!
      do k = 0, n_time_steps
         t = t_mesh(k)
         call mbe_calc_applfld(t,tw,pulse_type,t0,t1,cfield0,field_vec(k))
      enddo
!
      if(iinterp.eq.0 .and. force0)then
         print*,'mbe_calc_envelope: iinterp = 0 is incompatible with'
         print*,'forcing the applied to be zero at the start.'
         call mbe_error
      endif
      if(force0 .and. (pulse_type.eq.'Gs' .or. pulse_type.eq.'hs'))then
         field_vec(0) = (0.d0,0.d0)
      endif
!
      end subroutine mbe_calc_envelope

!!!!!!!!!!!!!!!!!!!
 
      end subroutine mbe_set_tdfields_A

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      subroutine mbe_set_tdfields_B(n_time_steps_in,t_mesh_in,probef,couplf)
!
!  Pass the complex field amplitudes of the time-dependent fields to
!  mbe and copy those into global variables accessible throughout the module.
!
!  The argument couplf (the time-dependent complex amplitude of the
!  coupling field) is optional. If not present, mbe assumes that the
!  coupling field has zero amplitude at all times.
!
!  n_time_steps_in + 1 is the number of values of t, probef and couplf
!  passed on to mbe. The arrays t_mesh_in, probef and couplf
!  must have a dimension of at least n_time_steps + 1.
!
      implicit none
      double complex, dimension(:), intent(in) :: probef
      double complex, dimension(:), optional, intent(in) :: couplf
      double precision, dimension(:), intent(in) :: t_mesh_in
      integer, intent(in) :: n_time_steps_in
      integer :: imin,imax,istat,k
!
      imin = lbound(t_mesh_in,1)
      imax = ubound(t_mesh_in,1)
      if(lbound(probef,1).ne.imin .or. ubound(probef,1).ne.imax .or.  &
               (present(couplf) .and. lbound(couplf,1).ne.imin) .or.  &
               (present(couplf) .and. ubound(couplf,1).ne.imax))then
         print*,'mbe_set_tdfields_B detected an incorrect array bound.'
         call mbe_error
      endif
!
      if(n_time_steps_in.gt.(imax-imin))then
         print*,'mbe_set_tdfields_B: The value of n_time_steps is inconsistent'
         print*,'with the dimensions of the arrays.'
         print*,n_time_steps_in,imin,imax
         call mbe_error
      elseif(n_time_steps_in.le.0)then
         print*,'mbe_set_tdfields_B: The value of n_time_steps is nonsensical.'
         print*,n_time_steps_in
         call mbe_error
      endif
!
      n_time_steps = n_time_steps_in
!
      if(allocated(Field_p_vec))deallocate(Field_p_vec)
      if(allocated(Field_c_vec))deallocate(Field_c_vec)
      if(allocated(t_mesh))deallocate(t_mesh)
      allocate(Field_p_vec(0:n_time_steps),                 &
               Field_c_vec(0:n_time_steps),                 &
               t_mesh(0:n_time_steps), stat = istat)
      if(istat.ne.0)then
         print*,'mbe_set_tdfields_B: error at the allocation of the arrays.'
         call mbe_error
      endif
!
      do k = 0, n_time_steps
         Field_p_vec(k) = probef(imin+k)
         if(present(couplf))then
            Field_c_vec(k) = couplf(imin+k)
         else
            Field_c_vec(k) = (0.d0,0.d0)
         endif
         t_mesh(k) = real(t_mesh_in(imin+k),kd)
      enddo
!
!  The time dependence of the applied field will be calculated through
!  spline interpolation later on (fl_spline is a global variable within
!  mbe).
      fl_spline = .true.
!
      fl_notinit_fields = .false.
!
      end subroutine mbe_set_tdfields_B

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine mbe_calc_susc(Field_re,Field_im,sumrho_re,sumrho_im, &
                               ak,density,refr_index,alpha)
!
!  Calculate the complex susceptibility and returns the (real) refractive
!  index and the absorption coefficient (a.k.a. the attenuation coefficient).
!
      use obe_constants
      implicit none
      real(kd), intent(in) :: Field_re,Field_im,sumrho_re,sumrho_im,ak,density
      real(kd), intent(out) :: alpha,refr_index
      complex(kd) :: chi,Field,sumrho
!
      Field = cmplx(Field_re,Field_im,kd)
!  sumrho is meant to be the sum of the relevant complex coherences.
      sumrho = cmplx(sumrho_re,sumrho_im,kd)
!
!  Calculate the refractive index and the absorption coefficient. These
!  two quantities are set to arbitrary values if the electric field is zero.
!  ak is the vacuum wave number.
      if(Field_re**2+Field_im**2 .gt. 0.0_kd)then
!  No minus sign here as the minus sign is already in sumrho.
         chi =  2.0_kd * density * sumrho / (epsilon0 * Field)
         refr_index = real(sqrt((1.0_kd,0.0_kd)+chi),kd)
         alpha = 2.0_kd * ak * aimag(sqrt((1.0_kd,0.0_kd)+chi))
      else
         refr_index = -1.0e3_kd   ! Arbitrary value
         alpha = 0.0_kd           ! Arbitrary value, convenient for the
                                  ! propagation with izrule.eq.0.
      endif
!
      end subroutine mbe_calc_susc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine mbe_tdint_2(istart,nsubsteps,imethod,iDoppler,    &
               rhovec_av,iunit,iunit_ampl,popinit,icont,nsb_var)
!
!  istart: as for mbe_propagate_2, with the additional possibility of
!  the option of -1 for istart, which means that the initial density matrix,
!  for all the velocity classes (if Doppler averaging), is taken to be the 
!  content of rhovec_av at entry. In this case, the argument popinit
!  does not need to be present (the content of popinit is ignored if
!  this optional argument is present).
!
!  icont (optional): If .eq. 1, simply re-use the existing ommats arrays.
!  If .ne. 1, reset fl_init_ommats to .true. to force the re-calculation of
!  these arrays
!
!  iDoppler.eq.1 is the signal that the calculation involved an averaging
!  over velocity classes. The calculation then uses the abcissas and weights
!  (multiplied by the relevant Maxwell frunction) contained in the arrays
!  vmesh and fMvweight of the obe module. The variable iDoppler should be
!  set to 0 if no Doppler averaging is intended.
!
!  iunit: as for ldbl_int. In short, iunit is the unit number of the file used
!  to write out the results if iunit > 0, results are written on special
!  files if iunit.eq.-1, and nothing is written if iunit.eq.0. The density
!  matrix (or Doppler-averaged density matrix) at the last time step is
!  returned to the calling program through rhovec_av at the end of the
!  calculation.
!
!  The subroutine mbe_set_tdfields_A (or mbe_set_tdfields_B), which sets
!  the time-dependent fields as well as n_time_steps and the time mesh,
!  must have been called prior to a call to mbe_tdint_2. The time mesh
!  (in the global array t_mesh) is expected to be specified in microseconds.
!
!  nsubsteps and nsb_var: as for mbe_propagate_2.
!
!  Note: Calculations in the rate approximation are not yet supported by this
!  subroutine.
!
      implicit none
      integer, parameter :: neqs = nst*nst
      double precision, dimension (nst), optional, intent(in) :: popinit
      integer, intent(in), optional :: icont
      integer, dimension (:), intent(in), optional :: nsb_var
      integer, intent(in) :: iDoppler,imethod,istart,iunit,iunit_ampl,nsubsteps
      real(kd), dimension (neqs), intent(inout) :: rhovec_av
      real(kd), dimension (:,:), allocatable :: rhovec
      real(kd), dimension (:), allocatable :: rhovec_0
      real(kd), dimension (:,:), allocatable :: rhs_mat
      integer, dimension (:), allocatable :: nsbsteps
      integer :: i,istat,iv,ivmax,j,jmin,k,kmin,m,mrl,mim,nvsize
      character(len=11) :: file_format
      character(len=24) :: fmt
      logical :: fl
      logical, save :: fl_Doppler_old = .false.
!
!  Initialize variables controlling the integration, if dop853 is
!  to be used. rtol_dop853, atol_dop853 and fl_set_tol_dop853 are
!  global variables.
      if(imethod.eq.4)then
         call obe_get_tol_dop853(rtol_dop853,atol_dop853,fl_set_tol_dop853)
      endif
!
!  Set fl_twofields to .true. as this subroutine is for a pair of fields.
      fl_twofields = .true.
!
!  Chek that nsubsteps and nsb_var (if present) make sense.
      if(nsubsteps.eq.-1)then
         if(present(nsb_var))then
            if(size(nsb_var).lt.n_time_steps)then
               print*,'mbe_tdint_2: the size of nsb_var is inconsistent'
               print*,'with the value of n_time_steps passed to'
               print*,'mbe_set_tdfields.'
               call mbe_error
            else
               kmin = lbound(nsb_var,1)
            endif
         else
            print*,'mbe_tdint_2: nsubsteps = -1 but nsb_var is not present.'
            call mbe_error
         endif
      elseif(nsubsteps.lt.1)then
         print*,'mbe_tdint_2: nsubsteps is nonsensical. ',nsubsteps
         call mbe_error
      endif
      if(allocated(nsbsteps))deallocate(nsbsteps)
      allocate(nsbsteps(n_time_steps),stat=istat)
      if(istat.gt.0)then
         print*,'mbe_tdint_2: Error at the allocation of nsbsteps.'
         call mbe_error
      endif
      do k = 1,n_time_steps
         if(nsubsteps.gt.0)then
            nsbsteps(k) = nsubsteps
         else
            nsbsteps(k) = nsb_var(kmin + k - 1)
            if(nsbsteps(k).lt.1)then
               print*,'mbe_tdint_2: Nonsensical value in nsb_var.'
               print*,kmin,k,kmin+k-1,nsb_var(kmin+k-1)
               call mbe_error
            endif
         endif
      enddo
!
!  Check that popinit is present if this variable will be used.
      if(.not.present(popinit) .and. istart.ge.0)then
         print*,'mbe_tdint_2: popinit should be present unless istart.eq.-1.'
         call mbe_error
      endif
!
!  Check that the values of iunit and iunit_ampl make sense, and find out
!  whether the corresponding files (if any) are meant to be formatted or
!  unformatted. The field receiving the details of the applied fields is
!  assumed always to be formatted.
      if(iunit.gt.0)then
         if(iunit.eq.5 .or. iunit.eq.6)then
            print*,'Wrong unit number in mbe_tdint_2 ',iunit
            call mbe_error
         endif
         inquire(unit=iunit,exist=fl)
         if(.not.fl)then
            print*,'mbe_tdint_2: the unit iunit does not exist ',iunit
            call mbe_error
         endif
         inquire(unit=iunit,opened=fl)
         if(.not.fl)then
            print*,'mbe_tdint_2: the unit iunit is not connected',iunit
            call mbe_error
         endif
         inquire(unit=iunit,form=file_format)
         if(nst*nst+2.le.999999)then
            write(fmt,'("(1x,",i6,"(1pe18.11,1x))")')nst*nst+2
         else
            print*,'mbe_tdint_2: The number of states is larger than'
            print*,'what this subroutine can cope with in regards to'
            print*,'preparing the format fmt. The code'
            print*,'should be updated in the unlikely event that'
            print*,'nst*nst + 2 > 999,999.'
            call mbe_error
         endif
      elseif(iunit.eq.-1)then
         call mbe_set_outputfiles
      elseif(iunit.lt.-1)then
         print*,'mbe_tdint_2: illegal value of iunit. ',iunit
         call mbe_error
      endif
!
      if(iunit_ampl.gt.0)then
         if(iunit_ampl.eq.5 .or. iunit_ampl.eq.6)then
            print*,'Wrong unit number in mbe_tdint_2 for iunit_ampl',iunit_ampl
            call mbe_error
         endif
         inquire(unit=iunit_ampl,exist=fl)
         if(.not.fl)then
            print*,'mbe_tdint_2: the unit iunit_ampl does not exist ', &
                   iunit_ampl
            call mbe_error
         endif
         inquire(unit=iunit_ampl,opened=fl)
         if(.not.fl)then
            print*,'mbe_tdint_2: the unit iunit_ampl is not connected', &
                   iunit_ampl
            call mbe_error
         endif
         if(iunit.eq.iunit_ampl .and. iunit.ne.0 .and. iunit.ne.-1)then
            print*,'mbe_tdint_2: Conflict between iunit and iunit_ampl ', &
                   iunit,iunit_ampl
            call mbe_error
         endif
         if(fl_outputfiles)then
            do j = 1,nst
               do i = 1,nst
                  if(iunits_arr(i,j).eq.iunit_ampl)then
                     print*,'mbe_tdint_2: Conflict between iunit_ampl and'
                     print*,'iunits_arr for i,j = ',i,j
                     call mbe_error
                  endif
               enddo
            enddo
         endif
      endif
!
!  Check that the numerical integration over the velocity classes
!  has been initialised though a prior call to obe_set_Doppler, if the
!  calculation involves Doppler averaging. The logical variable
!  fl_Doppler_notinit is imported from the obe module.
      if(iDoppler.eq.1)then
         call mbe_set_Doppler
         if(fl_Doppler_notinit)then
            print*,'mbe_tdint_2: attempt at a calculation involving'
            print*,'a Doppler averaging before a call to obe_set_Doppler.'
            call mbe_error
         else
            fl_Doppler = .true.
            ivmax = n_v_values
         endif
      elseif(iDoppler.eq.0)then
         fl_Doppler = .false.
         ivmax = 1
      else
         print*,'mbe_tdint_2: illegal value of iDoppler ',iDoppler
         call mbe_error
      endif
!
!  Initialise rhsmat_0
      call mbe_set_rhsmat_0
!
!  Signal to re-initialize the arrays set by mbe_set_rhsmat (the logical
!  variable fl_init_ommats is global within the module), unless we continue
!  with the arrays initialized at a previous run.
      fl_init_ommats = .true.
      if(present(icont))then
         if(icont.eq.1)then
            if(fl_Doppler.neqv.fl_Doppler_old)then
               print*,'mbe_tdint_2: a continuation run with a different'
               print*,'Doppler option is not allowed.'
               call mbe_error
            endif
            fl_init_ommats = .false.
         endif
      endif
      fl_Doppler_old = fl_Doppler
!
!  Check that the complex field amplitudes have been provided (they must
!  be provided through a call to mbe_set_tdfields) and if they have been
!  transcribe them into the relevant arrays.
      if(fl_notinit_fields)then
         print*,'mbe_tdint_2 was called before a call to either'
         print*,'mbe_set_tdfields_A or mbe_set_tdfields_B.'
         call mbe_error
      endif
!  Write the input fields into the relevant arrays.
      if(allocated(Field_p_re_vec))deallocate(Field_p_re_vec)
      if(allocated(Field_p_im_vec))deallocate(Field_p_im_vec)
      if(allocated(Field_c_re_vec))deallocate(Field_c_re_vec)
      if(allocated(Field_c_im_vec))deallocate(Field_c_im_vec)
      allocate(Field_p_re_vec(0:n_time_steps),                 &
               Field_p_im_vec(0:n_time_steps),                 &
               Field_c_re_vec(0:n_time_steps),                 &
               Field_c_im_vec(0:n_time_steps), stat = istat)
      if(istat.ne.0)then
         print*,'mbe_tdint_2: Error at the allocation of the arrays (1).'
         call mbe_error
      endif
      do k = 0, n_time_steps
         Field_p_re_vec(k) = real(dreal(Field_p_vec(k)),kd)
         Field_p_im_vec(k) = real(dimag(Field_p_vec(k)),kd)
         Field_c_re_vec(k) = real(dreal(Field_c_vec(k)),kd)
         Field_c_im_vec(k) = real(dimag(Field_c_vec(k)),kd)
      enddo
!
!  Allocate the rest of the storage:
!
      allocate(rhovec(neqs,ivmax),       &
               rhs_mat(neqs,neqs),       &
               stat=istat)
      if(istat.gt.0)then
         print*,'mbe_tdint_2: Error at the allocation of the arrays (2).'
         call mbe_error
      endif
!
!  Prepare the spline interpolation used in the calculation of the density
!  matrix, if necessary
      if(fl_spline)then
         allocate(Fd_pre2(0:n_time_steps),  &
                  Fd_pim2(0:n_time_steps),  &
                  Fd_cre2(0:n_time_steps),  &
                  Fd_cim2(0:n_time_steps),  &
                  stat=istat)
         if(istat.gt.0)then
            print*,'mbe_tdint_2: Error at the allocation of the arrays (3).'
            call mbe_error
         endif
         call mbe_spline_prep(Field_p_re_vec,Fd_pre2)
         call mbe_spline_prep(Field_p_im_vec,Fd_pim2)
         call mbe_spline_prep(Field_c_re_vec,Fd_cre2)
         call mbe_spline_prep(Field_c_im_vec,Fd_cim2)
      endif
!
!  Allocation of storage used in case of calculations using the
!  rate equation approximation for the initial time. The arrays arhs and
!  ipiv_r are global within this module.
!
      if(istart.eq.1)then
         if(allocated(arhs))deallocate(arhs)
         if(allocated(ipiv_r))deallocate(ipiv_r)
         allocate(arhs(neqs,neqs),ipiv_r(neqs),stat=istat)
         if(istat.gt.0)then
            print*,'mbe_tdint_2: Error at the allocation of arhs or ipiv_r.'
            call mbe_error
         endif
      endif
!
!  Propagate the density matrix in time.
!  The initial condition is governed by the variable istart, which is passed to
!  mbe_tdint_2 and must be set to 0, 1, 2 or -1 by the calling program. If
!  istart is 2, the atoms are assumed to be initially in the steady-state they
!  evolve into if starting in the state with a non-zero initial population.
!  Otherwise, the atoms are assumed to be initially in the mixed state defined
!  by the population array popinit: if istart is 0, the coherences are zero,
!  whereas if istart is 1 the coherences are calculated in the rate equations
!  approximation. The initial density matrix is the initial content of
!  rhovec_av if istart = -1, and in this case popinit is not used if present.
!  Note that istart = -1 might not make any sense for a calculation with
!  Doppler averaging, as setting istart = -1 means that the initial density
!  matrix is the same for all the velocity classes.
      if(istart.eq.-1)then
         if(allocated(rhovec_0))deallocate(rhovec_0)
         allocate(rhovec_0(neqs),stat=istat)
         if(istat.gt.0)then
            print*,'mbe_tdint_2: Error at the allocation of rhovec_0.'
            call mbe_error
         endif
         rhovec_0 = rhovec_av
      endif
!
      do k = 0,n_time_steps
!
!  Loop over velocity (runs only if the calculation involves a Doppler average)
!
         if(fl_Doppler)rhovec_av = 0.0_kd
         do iv = 1,ivmax
            if(fl_Doppler)then
               current_vel = vmesh(iv)
            else
               current_vel = 0.0_kd ! The value of current_vel is relevant 
                                    ! only when fl_Doppler is .true.
            endif
!
            if(k.eq.0)then
               Field_p_re = Field_p_re_vec(0)
               Field_p_im = Field_p_im_vec(0)
               Field_c_re = Field_c_re_vec(0)
               Field_c_im = Field_c_im_vec(0)
               call mbe_set_rhsmat2(rhs_mat,0)
               if(istart.eq.-1)then
                  rhovec(:,iv) = rhovec_0
               else
                  rhovec(:,iv) = 0.0_kd
                  do i = 1,nst
                     call ldbl_pop_index(i,m)
                     rhovec(m,iv) = popinit(i)
                  enddo
! Unless istart is 1 or 2, the coherences are initially zero.
                  if(istart.eq.2)then
! Steady state
                     call ldbl_steadystate(rhovec(:,iv),2,rhs_mat)
                  elseif(istart.eq.1)then
! Rate approximation
                     call mbe_calc_coher(rhs_mat,rhovec(:,iv),1)
                     call mbe_calc_coher(rhs_mat,rhovec(:,iv),0)
                  elseif(istart.ne.0)then
                     print*,'mbe_tdint_2: illegal value of istart ',istart
                     call mbe_error
                  endif
               endif
!
            else
!
               ktime = k - 1
               call mbe_intobe(imethod,t_mesh(k-1),t_mesh(k),         &
                                           nsbsteps(k),rhovec(:,iv))
!
            endif
!
!  Averaging over velocity classes if needed.
            if(fl_Doppler)then
               rhovec_av = rhovec_av + fMvweight(iv) * rhovec(:,iv)
            else
               rhovec_av = rhovec(:,iv)
            endif
!
         enddo
!
         if(iunit.gt.0 .or. iunit.eq.-1 .or. iunit_ampl.gt.0)then
            call mbe_output_td(k,t_mesh(k))
         endif
!
      enddo
!
!  Reset fl_twofields to its default value.
      fl_twofields = .true.
!
      deallocate(Field_p_re_vec,Field_p_im_vec,Field_c_re_vec,Field_c_im_vec)
      if(fl_spline)deallocate(Fd_pre2,Fd_pim2,Fd_cre2,Fd_cim2)
      deallocate(rhovec,rhs_mat)
!
      contains

!!!!!!!!!!!!!!!!!!!

      subroutine mbe_output_td(k,t)
!
!  Save the results on file. A value of 0 for the
!  input variable k is a signal to rewind the unit iunit (if used) before
!  writing.
!
      implicit none
      real(kd), intent(in) :: t
      integer, intent(in) :: k
      integer :: i,iost,j,m,mre,mim
!
      if(iunit.gt.0)then
!
         if(k.eq.0)rewind iunit
         if(file_format.eq.'UNFORMATTED')then
            write(unit=iunit,iostat=iost)t,(rhovec_av(j),j=1,neqs)
         else
            write(iunit,fmt,iostat=iost)t,(rhovec_av(j),j=1,neqs)
         endif
         if(iost.ne.0)then
            print*,'mbe_output_td: write returns with iostat = ',iost
            call mbe_error
         endif
!
      elseif(iunit.eq.-1)then
!
         if(fl_outputfiles)then
            do j = 1,nst
               do i = 1,nst
                  if(iunits_arr(i,j).gt.0)then
                     if(i.eq.j)then
                        call ldbl_pop_index(i,m)
                        write(iunits_arr(i,j),1501)t,rhovec_av(m)
                     elseif(i.lt.j)then
                        call ldbl_coher_index(i,j,mre,mim)
                        write(iunits_arr(i,j),1502)t,rhovec_av(mre), &
                                                     rhovec_av(mim)
                     else
                        call ldbl_coher_index(j,i,mre,mim)
                        write(iunits_arr(i,j),1502)t,rhovec_av(mre), &
                                                    -rhovec_av(mim)
                     endif
                  endif
               enddo
            enddo
         else
            print*,'mbe_output_td: The value of iunit is -1 but the unit'
            print*,'numbers of the output files have not been provided.'
            call mbe_error
         endif
      endif
!
      if(iunit_ampl.gt.0)then
!
         call mbe_calc_flds(t)
         write(iunit_ampl,1602)t,Field_p_re,Field_p_im,Field_c_re,Field_c_im
!
      endif
!
 1501 format(1x,2(1pe12.5,2x))
 1502 format(1x,3(1pe12.5,2x))
 1602 format(1x,5(1pe12.5,2x))
!
      end subroutine mbe_output_td

!!!!!!!!!!!!!!!!!!!

      end subroutine mbe_tdint_2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine mbe_propagate_1(istart,popinit,imethod,nsubsteps,         &
               izrule,zmax_in,n_z_steps,fieldpr,                           &
               density_in,iDoppler,nz_writeout,nt_writeout,output_1_pr,    &
               icont,nsb_var)
!
!  The same as mbe_propagate_2 but for the propagation of the probe field
!  only. The only difference to note with respect to mbe_propagate_2 is
!  the different structure of the calling list of output_pr.
!
!  See mbe_propagate_2 for further information.
!
      use obe_constants
      implicit none
      integer, parameter :: neqs = nst*nst
      double precision, dimension (nst), intent(in) :: popinit
      double precision, intent(in) :: density_in,zmax_in
      integer, dimension (:), intent(in), optional :: nsb_var
      integer, intent(in) :: imethod,istart,izrule,nsubsteps,       &
         nt_writeout,nz_writeout, n_z_steps, iDoppler
      integer :: npr
      integer, intent(in), optional :: icont
      type(obecfield), intent(in) :: fieldpr
      real(kd), dimension(:), allocatable :: z_mesh
      real(kd) :: density,zmax
      real(kd) :: atten_p,t,time_step,zmin
      real(kd), dimension (:), allocatable :: Fd_p_re_prev
      real(kd), dimension (:), allocatable :: Fd_p_im_prev
      real(kd), dimension (:), allocatable :: alpha_p,refr_index_p
      real(kd), dimension (:,:), allocatable :: rhs_mat
      real(kd), dimension (:,:), allocatable :: rhsprre,rhsprim
      real(kd), dimension (:,:), allocatable :: sumrho_pre,sumrho_pim
      real(kd), dimension(:), allocatable :: factpr
      integer, dimension (:), allocatable :: idxprre,idxprim
      real(kd), dimension (neqs) :: rhovec
      real(kd) :: akpr,hzennglp,unit_l,pi,z_step
      integer, dimension (:), allocatable :: nsbsteps
      integer :: i,istat,iv,ivmax,j,jmin,k,kmin,m,mrl,mim,nvsize
      logical, save :: fl_Doppler_old = .false.
      external :: output_1_pr
!
      if(.not.fl_spline)then
         print*,'mbe_propagate_1 requires interpolation...'
         call mbe_error
      endif
!
!  Initialize variables controlling the integration, if dop853 is
!  to be used. rtol_dop853, atol_dop853 and fl_set_tol_dop853 are
!  global variables.
      if(imethod.eq.4)then
         call obe_get_tol_dop853(rtol_dop853,atol_dop853,fl_set_tol_dop853)
      endif
!
!  Set fl_twofields to .false. as this subroutine is for a single field.
      fl_twofields = .false.
!
      if(istart.lt.0 .or. istart.gt.2)then
         print*,'mbe_propagate_1: illegal value of istart. ',istart
         call mbe_error
      endif
!
!  Chek that nsubsteps and nsb_var (if present) make sense.
      if(nsubsteps.eq.-1)then
         if(present(nsb_var))then
            if(size(nsb_var).lt.n_time_steps)then
               print*,'mbe_propagate_1: the size of nsb_var is inconsistent'
               print*,'with the value of n_time_steps passed to'
               print*,'mbe_set_tdfields.'
               call mbe_error
            else
               kmin = lbound(nsb_var,1)
            endif
         else
            print*,'mbe_propagate_1: nsubsteps = -1 but nsb_var is not present.'
            call mbe_error
         endif
      elseif(nsubsteps.lt.1)then
         print*,'mbe_propagate_1: nsubsteps is nonsensical. ',nsubsteps
         call mbe_error
      endif
      if(allocated(nsbsteps))deallocate(nsbsteps)
      allocate(nsbsteps(n_time_steps),stat=istat)
      if(istat.gt.0)then
         print*,'mbe_propagate_1: Error at the allocation of nsbsteps.'
         call mbe_error
      endif
      do k = 1,n_time_steps
         if(nsubsteps.gt.0)then
            nsbsteps(k) = nsubsteps
         else
            nsbsteps(k) = nsb_var(kmin + k - 1)
            if(nsbsteps(k).lt.1)then
               print*,'mbe_propagate_1: Nonsensical value in nsb_var.'
               print*,kmin,k,kmin+k-1,nsb_var(kmin+k-1)
               call mbe_error
            endif
         endif
      enddo
!
!  Convert to kd precision
      density = real(density_in,kd)
      zmax = real(zmax_in,kd)
!
!  Check that the numerical integration over the velocity classes
!  has been initialised though a prior call to obe_set_Doppler.
      if(iDoppler.eq.1)then
         call mbe_set_Doppler
         if(fl_Doppler_notinit)then
            print*,'mbe_propagate_1: attempt at a calculation involving'
            print*,'a Doppler averaging before a call to mbe_set_Doppler.'
            call mbe_error
         else
            fl_Doppler = .true.
            ivmax = n_v_values  ! from obe
         endif
      elseif(iDoppler.eq.0)then
         fl_Doppler = .false.
         ivmax = 1
      else
         print*,'mbe_propagate_1: illegal value of iDoppler ',iDoppler
         call mbe_error
      endif
!
!  Initialise rhsmat_0
      call mbe_set_rhsmat_0
!
!  Define the factors used in the calculation
!
!  First for the probe field.
!
!  Start by finding the number of transitions which need to be taken into
!  account and preparing the idx arrays. The content of fieldpr%dip_mom
!  will already have been checked by obe_setfields, which avoids the
!  possibility that a same transition would be counted twice (once as
!  i,j, once as j,i), but for safety we check again.
      npr = 0
      do i = 1,nst
         do j = 1,nst
            if(fieldpr%dip_mom(j,i).ne.(0.d0,0.d0))then
               if(fieldpr%dip_mom(i,j).ne.(0.d0,0.d0))then
                  print*,'mbe_propagate_1: error for the probe field,'
                  print*,'the array dip_mom is incorrectly formatted.'
                  print*,j,i,fieldpr%dip_mom(j,i),fieldpr%dip_mom(i,j)
                  call mbe_error
               endif
               npr = npr + 1
            endif
         enddo
      enddo
!
      if(allocated(idxprre))deallocate(idxprre)
      if(allocated(idxprim))deallocate(idxprim)
      allocate(idxprre(npr),idxprim(npr),factpr(npr),stat=istat)
      if(istat.gt.0)then
         print*,'mbe_propagate_1: Error at the allocation of the'
         print*,'idxppre, idxprim and factpr arrays.'
         call mbe_error
      endif
!
      npr = 0
      do i = 1,nst
         do j = 1,nst
            if(fieldpr%dip_mom(j,i).ne.(0.d0,0.d0))then
               npr = npr + 1
               if(i.ge.j)then
                  print*,'mbe_propagate_1: error in dip_mom for the probe'
                  print*,'field, or error in the labelling of the states.'
                  print*,'The subroutine will work only if dip_mom(j,i) is'
                  print*,'non-zero only if state j is higher in energy than'
                  print*,'state i and j > i.', i, j
                  call mbe_error
               endif
               call ldbl_coher_index(i,j,idxprre(npr),idxprim(npr))
               if(real(dimag(fieldpr%dip_mom(j,i)),kd).ne.0.0_kd)then
                  print*,'mbe_propagate_1 cannot handle complex transition dipole'
                  print*,'matrix elements.'
                  print*,j,i,fieldpr%dip_mom(j,i),' (probe field)'
                  call mbe_error
               else
                  factpr(npr) = real(fieldpr%dip_mom(j,i),kd)
               endif
            endif
         enddo
      enddo
!
!  Wavenumber (wavelengths are meant to be specified in nm):
      pi = 4.0_kd*atan(1.0_kd)
      if(fieldpr%wavelength.eq.0.0_kd)then
         print*,'mbe_propagate_1: The wavelength of the probe field is zero.'
         call mbe_error
      endif
      akpr = 2.0_kd * pi / (fieldpr%wavelength*1.0e-9_kd)
!     
!  Signal to re-initialize the arrays set by mbe_set_rhsmat (the logical
!  variable fl_init_ommats is global within the module), unless we continue
!  with the arrays initialized at a previous run.
      fl_init_ommats = .true.
      if(present(icont))then
         if(icont.eq.1)then
            if(fl_Doppler.neqv.fl_Doppler_old)then
               print*,'mbe_propagate_1: a continuation run with a different'
               print*,'Doppler option is not allowed.'
               call mbe_error
            endif
            fl_init_ommats = .false.
         endif
      endif
      fl_Doppler_old = fl_Doppler
!
      if(nt_writeout.le.0 .or. nz_writeout.le.0)then
         print*,'mbe_propagate_1: the value of nz_writeout or nt_writeout is'
         print*,'nonsensical. ',nz_writeout,nt_writeout
         call mbe_error
      endif
!
!  Check that the complex field amplitudes have been provided (they must
!  be provided through a call to mbe_set_tdfields) and if they have been
!  transcribe them into the relevant arrays.
      if(fl_notinit_fields)then
         print*,'mbe_propagate_1 was called before a call to mbe_set_tdfields.'
         call mbe_error
      endif
!  Write the input fields into the relevant arrays.
      if(allocated(Field_p_re_vec))deallocate(Field_p_re_vec)
      if(allocated(Field_p_im_vec))deallocate(Field_p_im_vec)
      allocate(Field_p_re_vec(0:n_time_steps),                 &
               Field_p_im_vec(0:n_time_steps),                 &
               stat = istat)
      if(istat.ne.0)then
         print*,'mbe_propagate_1: Error at the allocation of the arrays (1).'
         call mbe_error
      endif
      do k = 0, n_time_steps
         Field_p_re_vec(k) = real(dreal(Field_p_vec(k)),kd)
         Field_p_im_vec(k) = real(dimag(Field_p_vec(k)),kd)
      enddo
!    
!  Allocate the rest of the storage:
!
      allocate(Fd_p_re_prev(0:n_time_steps),  &
               Fd_p_im_prev(0:n_time_steps),  &
               Fd_pre2(0:n_time_steps),  &
               Fd_pim2(0:n_time_steps),  &
               rhs_mat(neqs,neqs),            &
               stat=istat)
      if(istat.gt.0)then
         print*,'mbe_propagate_1: Error at the allocation of the arrays (2).'
         call mbe_error
      endif
      if(izrule.eq.0)then
         allocate(z_mesh(-1:n_z_steps),           &
                  alpha_p(0:n_time_steps),        &
                  refr_index_p(0:n_time_steps),   &
                  stat=istat)
      elseif(izrule.eq.2)then
         allocate(z_mesh(-1:n_z_steps),           &
                  rhsprre(0:n_time_steps,-2:-1),  &
                  rhsprim(0:n_time_steps,-2:-1),  &
                  stat=istat)
      elseif(izrule.eq.3)then
         allocate(z_mesh(-6:n_z_steps),          &
                  rhsprre(0:n_time_steps,-5:0),  &
                  rhsprim(0:n_time_steps,-5:0),  &
                  stat=istat)
      else
         print*,'mbe_propagate_1: Illegal value of izrule. ',izrule
         call mbe_error
      endif
      if(istat.gt.0)then
         print*,'mbe_propagate_1: Error at the allocation of the arrays (3).'
         call mbe_error
      endif
      allocate(sumrho_pre(0:n_time_steps,0:ivmax),  &
               sumrho_pim(0:n_time_steps,0:ivmax),  &
               stat=istat)
      if(istat.gt.0)then
         print*,'mbe_propagate_1: Error at the allocation of the arrays (4).'
         call mbe_error
      endif
!
!  Allocation of storage used in case of calculations using the
!  rate equation approximation for the initial time. The arrays arhs and
!  ipiv_r are global within this module.
!
      if(istart.eq.1)then
         if(allocated(arhs))deallocate(arhs)
         if(allocated(ipiv_r))deallocate(ipiv_r)
         allocate(arhs(neqs,neqs),ipiv_r(neqs),stat=istat)
         if(istat.gt.0)then
            print*,'mbe_propagate_1: Error at the allocation of arhs or ipiv_r.'
            call mbe_error
         endif
      endif
!
!  Space mesh. The unit of length is set to 1 mum for propagation.
!
      unit_l = 1.0e-6_kd
!
!  We start the loop over space at j = -1 or j = -6 to accommodate the
!  extra initial steps required by either the mid-point formula (j = -1)
!  or the 4th-order Runge-Kutta formula (j=-6) used between z = 0 and
!  either z = z_step or z = 2z_step.
!
      zmin = 0.0_kd
      z_step = (zmax - zmin)/real(n_z_steps,kd)
      if(izrule.eq.0 .or. izrule.eq.2)then
         jmin = -1
         z_mesh(-1) = zmin
         z_mesh(0) = zmin + z_step/2.0_kd
         do j = 1,n_z_steps
            z_mesh(j) = zmin+real(j,kd)*z_step
         enddo
      else
         if(n_z_steps.lt.2)then
            print*,'Bad value of n_z_steps in mbe_propagate_1'
            call mbe_error
         endif
         jmin = -6
         z_mesh(-6) = zmin
         z_mesh(-5) = zmin + z_step/2.0_kd
         z_mesh(-4) = zmin + z_step/2.0_kd
         z_mesh(-3) = zmin + z_step
         z_mesh(-2) = zmin + z_step
         z_mesh(-1) = zmin + z_step + z_step/2.0_kd
         z_mesh( 0) = zmin + z_step + z_step/2.0_kd
         z_mesh( 1) = zmin + z_step + z_step
         do j = 2,n_z_steps
            z_mesh(j) = zmin+real(j,kd)*z_step
         enddo
      endif
!
      if(izrule.eq.0)then
!
!  Initialization of the calculation without propagation: the absorption
!  coefficient and refractive index at the very front of the medium
!  are set to their vacuum values.
!
         alpha_p = 0.0_kd
         refr_index_p = 1.0_kd
!
      else
!
!  Propagation parameter (the density is meant to be specified in
!  number of atoms / m^3).
!  
      hzennglp = z_step * (akpr * density / epsilon0) * unit_l
!
!  Minus sign to compensate for an erroneous minus sign in the
!  internal subprograms mbe_calc_Fd_2 and mbe_calc_Fd_3.
      hzennglp = -hzennglp
!
      endif
!
!
!  Start the loop over the space steps:
!
      do j = jmin,n_z_steps
!
!  Get the complex field amplitudes at that point in the medium. The
!  algorithm will depend on the value of izrule and whether we do a real
!  propagation calculation or only use Beer's Law (the latter is used when
!  izrule.eq.0).
!
         if(izrule.eq.0)then
!
!  This part of the code needs revisiting...
            print*,'mbe_propagate_2: calculations with izrule.eq.0 might be'
            print*,'incorrect in regards to the variation of the phase'
            print*,'of the field, as only the modulus of the amplitude'
            print*,'is changed, which amounts to positing that Re chi = 0.'
            print*,'Also, the calculation of the complex susceptibility might'
            print*,'need correction and re-alignment on obe_susceptiblity.'
            print*,'This part of the code is thus unsafe pending re-analysis.'
            call mbe_error
!  Beer's Law (no change for j.eq.jmin, which is not really a propagation step)
            if(j.gt.jmin)then
               do k = 0,n_time_steps
                  atten_p = exp(-alpha_p(k) * refr_index_p(k) *        &
                                (z_mesh(j)-z_mesh(j-1))*unit_l/2.0_kd)
                  Field_p_re_vec(k) = Field_p_re_vec(k) * atten_p
                  Field_p_im_vec(k) = Field_p_im_vec(k) * atten_p
               enddo
            endif
         elseif(izrule.eq.2)then
            call mbe_calc_Fd_2
         else
            call mbe_calc_Fd_3
         endif
!
!  Prepare the spline interpolation used in the calculation of the density
!  matrix:
         call mbe_spline_prep(Field_p_re_vec,Fd_pre2)
         call mbe_spline_prep(Field_p_im_vec,Fd_pim2)
!
!  Loop over velocity (runs only if the calculation involves a Doppler average)
!
         do iv = 1,ivmax
            if(fl_Doppler)then
               current_vel = vmesh(iv)
            else
               current_vel = 0.0_kd ! The value of current_vel is relevant 
                                    ! only when fl_Doppler is .true.
            endif
!
!  Propagate the density matrix in time at the current point in the medium.
!  The initial condition is governed by the variable istart, which is passed to
!  mbe_propagate_2 and must be set to 0, 1 or 2 by the calling program. If istart
!  is 2, the atoms are assumed to be initially in the steady-state they evolve
!  into if starting in the state with a non-zero initial population. Otherwise,
!  the atoms are assumed to be initially in the mixed state defined by the
!  population array popinit: if istart is 0, the coherences are zero, whereas
!  if istart is 1 the coherences are calculated in the rate equations 
!  approximation.
!
            do k = 0,n_time_steps
!
               if(k.eq.0)then
                  rhovec(:) = 0.0_kd
                  do i = 1,nst
                     call ldbl_pop_index(i,m)
                     rhovec(m) = real(popinit(i),kd)
                  enddo
                  Field_p_re = Field_p_re_vec(0)
                  Field_p_im = Field_p_im_vec(0)
                  call mbe_set_rhsmat1(rhs_mat,0)
! Unless istart is 1 or 2, the coherences are initially zero.
                  if(istart.eq.2)then
! Steady state
                     call ldbl_steadystate(rhovec,2,rhs_mat)
                  elseif(istart.eq.1)then
! Rate approximation
                     call mbe_calc_coher(rhs_mat,rhovec,1)
                     call mbe_calc_coher(rhs_mat,rhovec,0)
                  endif
               else
!
                  ktime = k - 1
                  call mbe_intobe(imethod,t_mesh(k-1),t_mesh(k),   &
                                               nsbsteps(k),rhovec)
!
               endif
!
!  Sum the relevant coherences (for the propagation calculation)
               sumrho_pre(k,iv) = 0.0_kd
               sumrho_pim(k,iv) = 0.0_kd
               do i = 1,npr
                  sumrho_pre(k,iv) = sumrho_pre(k,iv) - factpr(i)* &
                                               rhovec(idxprre(i))
                  sumrho_pim(k,iv) = sumrho_pim(k,iv) + factpr(i)* &
                                               rhovec(idxprim(i))
               enddo
!
!  Write out the density matrix and the fields (except for the intermediate
!  step between zmin and zmin + z_step or zmin + 2 z_step required by the
!  rule used to start the spatial integration):
               if(      j.eq.jmin                                            &
                  .or. (izrule.eq.2 .and. j.eq. 1 .and. nz_writeout.eq.1)    &
                  .or. (izrule.eq.3 .and. j.eq.-2 .and. nz_writeout.eq.1)    &
                  .or. (j.gt.1 .and. mod(j,nz_writeout).eq.0) )then
                  if(mod(k,nt_writeout).eq.0)then
                     call output_1_pr(j,k,iv,z_mesh(j),t_mesh(k),rhovec,     &
                                    Field_p_re_vec(k),Field_p_im_vec(k))
                  endif
               endif
!
            enddo
!
         enddo
!
         call mbe_Doppler
!
!  If "no propagation", calculate the absorption coefficient at this
!  point in the medium and skip the rest. Otherwise, save the
!  sum / integrated coherences for use in the next space step.
!
         if(izrule.eq.0 .and. j.lt.n_z_steps)then
            do k = 0,n_time_steps
               call mbe_calc_susc(Field_p_re_vec(k),Field_p_im_vec(k),       &
                                  sumrho_pre(k,0),sumrho_pim(k,0),           &
                                  akpr,density,refr_index_p(k),alpha_p(k))
            enddo
         elseif(izrule.eq.2)then
            if(j.ne.-1 .and. j.ne.1)then
               rhsprre(:,-2) = rhsprre(:,-1)
               rhsprim(:,-2) = rhsprim(:,-1)
            endif
            rhsprre(:,-1) = hzennglp*sumrho_pre(:,0)
            rhsprim(:,-1) = hzennglp*sumrho_pim(:,0)
         elseif(izrule.eq.3)then
            if(j.eq.-6)then
               rhsprre(:,-5) = hzennglp*sumrho_pre(:,0)
               rhsprim(:,-5) = hzennglp*sumrho_pim(:,0)
               rhsprre(:,-4) = rhsprre(:,-5)
               rhsprim(:,-4) = rhsprim(:,-5) 
            elseif(j.eq.-5 .or. j.eq.-1)then
               rhsprre(:,-3) = hzennglp*sumrho_pre(:,0)
               rhsprim(:,-3) = hzennglp*sumrho_pim(:,0)
            elseif(j.eq.-4 .or. j.eq.0)then
               rhsprre(:,-2) = hzennglp*sumrho_pre(:,0)
               rhsprim(:,-2) = hzennglp*sumrho_pim(:,0)
            elseif(j.eq.-3 .or. j.eq.1)then
               rhsprre(:,-1) = hzennglp*sumrho_pre(:,0)
               rhsprim(:,-1) = hzennglp*sumrho_pim(:,0)
            elseif(j.eq.-2)then
               rhsprre(:,-4) = hzennglp*sumrho_pre(:,0)
               rhsprim(:,-4) = hzennglp*sumrho_pim(:,0)
            elseif(j.eq.2)then
               rhsprre(:,-3) = rhsprre(:,-5)
               rhsprim(:,-3) = rhsprim(:,-5) 
               rhsprre(:,-2) = rhsprre(:,-4)
               rhsprim(:,-2) = rhsprim(:,-4) 
               rhsprre(:,-1) = hzennglp*sumrho_pre(:,0)
               rhsprim(:,-1) = hzennglp*sumrho_pim(:,0)
            elseif(j.gt.2)then
               rhsprre(:,-3) = rhsprre(:,-2)
               rhsprim(:,-3) = rhsprim(:,-2)
               rhsprre(:,-2) = rhsprre(:,-1)
               rhsprim(:,-2) = rhsprim(:,-1)
               rhsprre(:,-1) = hzennglp*sumrho_pre(:,0)
               rhsprim(:,-1) = hzennglp*sumrho_pim(:,0)
            endif
         endif
!
      enddo
!
      deallocate(Field_p_re_vec,Field_p_im_vec)
      deallocate(Fd_p_re_prev,Fd_p_im_prev,                                   &
         Fd_pre2,Fd_pim2,rhs_mat,z_mesh,                                      &
         sumrho_pre,sumrho_pim)
      if(izrule.eq.0)then
         deallocate(alpha_p,refr_index_p)
      else
         deallocate(rhsprre,rhsprim)
      endif
!
      contains
 
!!!!!!!!!!!!!!!!!!!
 
      subroutine mbe_calc_Fd_2
!
!  Warning: the right-hand sides are multiplied by an erroneous -1 factor.
!  This is compensated by a change of sign of hzennglp.
!
!  Calculate the new field amplitude. Use the mid point formula for the 
!  first step and the 2nd order Adams-Bashforth for the rest.
!
!  Entrance of the medium (j = jmin): The field is the applied fields.
!                                     mbe_propagate_1 expects to find the
!                                     corresponding complex field
!                                     already in Field_p_re_vec etc.
!  First step into the medium: Use the mid-point formula to calculate
!                              the field.
!  Subsequent steps: Use the Adams-Bashforth method.
!
      if(j.eq.-1)then
!  Use the existing contents of Field_p_re_vec etc. Save these values
!  for later use in the midpoint formula.
         Fd_p_re_prev = Field_p_re_vec
         Fd_p_im_prev = Field_p_im_vec
      elseif(j.eq.0)then
!  Midpoint formula (first half step)
         Field_p_re_vec = Fd_p_re_prev - rhsprim(:,-1)*0.5_kd
         Field_p_im_vec = Fd_p_im_prev + rhsprre(:,-1)*0.5_kd
      elseif(j.eq.1)then
!  Midpoint formula (second half step)
         Field_p_re_vec = Fd_p_re_prev - rhsprim(:,-1)
         Field_p_im_vec = Fd_p_im_prev + rhsprre(:,-1)
      elseif(j.ge.2)then
!  Adams-Bashforth step
         Field_p_re_vec = Field_p_re_vec - rhsprim(:,-1)*1.5_kd &
                                         + rhsprim(:,-2)*0.5_kd
         Field_p_im_vec = Field_p_im_vec + rhsprre(:,-1)*1.5_kd &
                                         - rhsprre(:,-2)*0.5_kd
      else
!
         print*,'Logic error in mbe_propagate_1 or mbe_calc_Fd_2 ',j
         call mbe_error
!
      endif
!
      end subroutine mbe_calc_Fd_2

!!!!!!!!!!!!!!!!!!!
 
      subroutine mbe_calc_Fd_3
!
!  Warning: the right-hand sides are multiplied by an erroneous -1 factor.
!  This is compensated by a change of sign of hzennglp.
!
!  Calculate the new field amplitude. Use RK4 for the first two steps and
!  the 3rd order Adams-Bashforth for the rest.
!
!  Entrance of the medium (j = jmin): The field is the applied field.
!                                     mbe_propagate_1 expects to find the
!                                     corresponding complex field
!                                     already in Field_p_re_vec etc.
!  First two steps into the medium: Use RK4 to calculate the field.
!  Subsequent steps: Use the Adams-Bashforth method.
!
      if(j.eq.-6)then
!
         continue  ! Use the existing contents of Field_p_re_vec etc
!
      elseif(j.eq.-5 .or. j.eq.-1)then
!
         Fd_p_re_prev = Field_p_re_vec
         Fd_p_im_prev = Field_p_im_vec
!  RK4 (first step)
         Field_p_re_vec = Fd_p_re_prev - rhsprim(:,-4)
         Field_p_im_vec = Fd_p_im_prev + rhsprre(:,-4)
      elseif(j.eq.-4 .or. j.eq.0)then
!  RK4 (second step)
         Field_p_re_vec = Fd_p_re_prev - rhsprim(:,-3)*0.5_kd 
         Field_p_im_vec = Fd_p_im_prev + rhsprre(:,-3)*0.5_kd 
      elseif(j.eq.-3 .or. j.eq.1)then
!  RK4 (third step)
         Field_p_re_vec = Fd_p_re_prev - rhsprim(:,-2)*0.5_kd 
         Field_p_im_vec = Fd_p_im_prev + rhsprre(:,-2)*0.5_kd 
      elseif(j.eq.-2 .or. j.eq.2)then
!  RK4 (fourth step)
         Field_p_re_vec = Fd_p_re_prev - (       rhsprim(:,-4)+        &
                                          2.0_kd*rhsprim(:,-3)+        &
                                          2.0_kd*rhsprim(:,-2)+        &
                                                 rhsprim(:,-1)  )/6.0_kd
         Field_p_im_vec = Fd_p_im_prev + (       rhsprre(:,-4)+        &
                                          2.0_kd*rhsprre(:,-3)+        &
                                          2.0_kd*rhsprre(:,-2)+        &
                                                 rhsprre(:,-1)  )/6.0_kd
!
      elseif(j.ge.3)then
!
         Fd_p_re_prev = Field_p_re_vec
         Fd_p_im_prev = Field_p_im_vec
!
!  Adams-Bashforth step
!  3rd order
         Field_p_re_vec = Fd_p_re_prev - (23.0_kd*rhsprim(:,-1) - &
                                          16.0_kd*rhsprim(:,-2) + &
                                           5.0_kd*rhsprim(:,-3))/12.0_kd 
         Field_p_im_vec = Fd_p_im_prev + (23.0_kd*rhsprre(:,-1) - &
                                          16.0_kd*rhsprre(:,-2) + &
                                           5.0_kd*rhsprre(:,-3))/12.0_kd 
!  2nd order
!        Field_p_re_vec = Fd_p_re_prev - rhsprim(:,-1)*1.50_kd &
!                                      + rhsprim(:,-2)*0.50_kd
!        Field_p_im_vec = Fd_p_im_prev + rhsprre(:,-1)*1.50_kd &
!                                      - rhsprre(:,-2)*0.50_kd
!  1st order
!        Field_p_re_vec = Fd_p_re_prev - rhsprim(:,-1)
!        Field_p_im_vec = Fd_p_im_prev + rhsprre(:,-1)

!  Calculate the density matrix with the new field amplitudes:
!
!  Prepare the spline interpolation used in the calculation of the density
!  matrix:
         call mbe_spline_prep(Field_p_re_vec,Fd_pre2)
         call mbe_spline_prep(Field_p_im_vec,Fd_pim2)
!  Loop over velocity 
         do iv = 1,ivmax
            if(fl_Doppler)then
               current_vel = vmesh(iv)
            else
               current_vel = 0.0_kd
            endif
!  Propagate the density matrix in time at the current point in the medium.
            do k = 0,n_time_steps
               if(k.eq.0)then
                  rhovec(:) = 0.0_kd
                  do i = 1,nst
                     call ldbl_pop_index(i,m)
                     rhovec(m) = real(popinit(i),kd)
                  enddo
                  Field_p_re = Field_p_re_vec(0)
                  Field_p_im = Field_p_im_vec(0)
                  call mbe_set_rhsmat1(rhs_mat,0)
! Unless istart is 1 or 2, the coherences are initially zero.
                  if(istart.eq.2)then
! Steady state
                     call ldbl_steadystate(rhovec,2,rhs_mat)
                  elseif(istart.eq.1)then
! Rate approximation
                     call mbe_calc_coher(rhs_mat,rhovec,1)
                     call mbe_calc_coher(rhs_mat,rhovec,0)
                  endif
               else
                  ktime = k - 1
                  call mbe_intobe(imethod,t_mesh(k-1),t_mesh(k),   &
                                               nsbsteps(k),rhovec)
               endif
!  Sum the relevant coherences (for the propagation calculation)
               sumrho_pre(k,iv) = 0.0_kd
               sumrho_pim(k,iv) = 0.0_kd
               do i = 1,npr
                  sumrho_pre(k,iv) = sumrho_pre(k,iv) - factpr(i)* &
                                                rhovec(idxprre(i))
                  sumrho_pim(k,iv) = sumrho_pim(k,iv) + factpr(i)* &
                                                rhovec(idxprim(i))
               enddo
            enddo
         enddo
         call mbe_Doppler
         rhsprre(:,0) = hzennglp*sumrho_pre(:,0)
         rhsprim(:,0) = hzennglp*sumrho_pim(:,0)
!
!  Adams-Moulton step
!
!  4th order
         Field_p_re_vec = Fd_p_re_prev - ( 9.0_kd*rhsprim(:,0)  + &
                                         19.0_kd*rhsprim(:,-1) - &
                                          5.0_kd*rhsprim(:,-2) + &
                                               rhsprim(:,-3) )/24.0_kd !&
         Field_p_im_vec = Fd_p_im_prev + ( 9.0_kd*rhsprre(:,0)  + &
                                         19.0_kd*rhsprre(:,-1) - &
                                          5.0_kd*rhsprre(:,-2) + &
                                               rhsprre(:,-3) )/24.0_kd !& 
!  3rd order
!        Field_p_re_vec = Fd_p_re_prev - ( 5.0_kd*rhsprim(:,0)  + &
!                                         8.0_kd*rhsprim(:,-1) - &
!                                              rhsprim(:,-2) )/12.0_kd  
!        Field_p_im_vec = Fd_p_im_prev + ( 5.0_kd*rhsprre(:,0)  + &
!                                         8.0_kd*rhsprre(:,-1) - &
!                                              rhsprre(:,-2) )/12.0_kd  
!
      else
!
         print*,'Logic error in mbe_propagate_1 or mbe_calc_Fd_4 ',j
         call mbe_error
!
      endif
!
      end subroutine mbe_calc_Fd_3

!!!!!!!!!!!!!!!!!!!

      subroutine mbe_Doppler
!
      implicit none
      real(kd), dimension(:,:), allocatable :: sum_mod,sum_arg
      real(kd), dimension(:), allocatable :: sum_mod_v,sum_arg_v
      real(kd), dimension(:), allocatable :: s_pre,s_pim
      real(kd), dimension(:), allocatable :: vtest
      real(kd) :: v
      real(kd) :: tot1,tot2,ynt
      integer :: istat,j,numb
!
      if(.not.fl_Doppler)then
         sumrho_pre(:,0) = sumrho_pre(:,1)
         sumrho_pim(:,0) = sumrho_pim(:,1)
         return
      endif
!
!  Calculation with Doppler averaging.
!
      sumrho_pre(:,0) = 0.0_kd
      sumrho_pim(:,0) = 0.0_kd
      do i = 1,ivmax
         sumrho_pre(:,0) = sumrho_pre(:,0)+sumrho_pre(:,i)*fMvweight(i)
         sumrho_pim(:,0) = sumrho_pim(:,0)+sumrho_pim(:,i)*fMvweight(i)
      enddo
!
      end subroutine mbe_Doppler

!!!!!!!!!!!!!!!!!!!
 
      end subroutine mbe_propagate_1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine mbe_tdint_1(istart,nsubsteps,imethod,iDoppler,    &
               rhovec_av,iunit,iunit_ampl,popinit,icont,nsb_var)
!
!  The same as mbe_tdint_2 for a single field (the "probe field").
!
!  See mbe_tdint_2 for further information.
!
      implicit none
      integer, parameter :: neqs = nst*nst
      double precision, dimension (nst), optional, intent(in) :: popinit
      integer, intent(in), optional :: icont
      integer, dimension (:), intent(in), optional :: nsb_var
      integer, intent(in) :: iDoppler,imethod,istart,iunit,iunit_ampl,nsubsteps
      real(kd), dimension (neqs), intent(inout) :: rhovec_av
      real(kd), dimension (:,:), allocatable :: rhovec
      real(kd), dimension (:), allocatable :: rhovec_0
      real(kd), dimension (:,:), allocatable :: rhs_mat
      integer, dimension (:), allocatable :: nsbsteps
      integer :: i,istat,iv,ivmax,j,jmin,k,kmin,m,mrl,mim,nvsize
      character(len=11) :: file_format
      character(len=24) :: fmt
      logical :: fl
      logical, save :: fl_Doppler_old = .false.
!
!  Initialize variables controlling the integration, if dop853 is
!  to be used. rtol_dop853, atol_dop853 and fl_set_tol_dop853 are
!  global variables.
      if(imethod.eq.4)then
         call obe_get_tol_dop853(rtol_dop853,atol_dop853,fl_set_tol_dop853)
      endif
!
!  Set fl_twofields to .false. as this subroutine is for a single field.
      fl_twofields = .false.
!
!  Chek that nsubsteps and nsb_var (if present) make sense.
      if(nsubsteps.eq.-1)then
         if(present(nsb_var))then
            if(size(nsb_var).lt.n_time_steps)then
               print*,'mbe_tdint_1: the size of nsb_var is inconsistent'
               print*,'with the value of n_time_steps passed to'
               print*,'mbe_set_tdfields.'
               call mbe_error
            else
               kmin = lbound(nsb_var,1)
            endif
         else
            print*,'mbe_tdint_1: nsubsteps = -1 but nsb_var is not present.'
            call mbe_error
         endif
      elseif(nsubsteps.lt.1)then
         print*,'mbe_tdint_1: nsubsteps is nonsensical. ',nsubsteps
         call mbe_error
      endif
      if(allocated(nsbsteps))deallocate(nsbsteps)
      allocate(nsbsteps(n_time_steps),stat=istat)
      if(istat.gt.0)then
         print*,'mbe_tdint_1: Error at the allocation of nsbsteps.'
         call mbe_error
      endif
      do k = 1,n_time_steps
         if(nsubsteps.gt.0)then
            nsbsteps(k) = nsubsteps
         else
            nsbsteps(k) = nsb_var(kmin + k - 1)
            if(nsbsteps(k).lt.1)then
               print*,'mbe_tdint_1: Nonsensical value in nsb_var.'
               print*,kmin,k,kmin+k-1,nsb_var(kmin+k-1)
               call mbe_error
            endif
         endif
      enddo
!
!  Check that popinit is present if this variable will be used.
      if(.not.present(popinit) .and. istart.ge.0)then
         print*,'mbe_tdint_1: popinit should be present unless istart.eq.-1.'
         call mbe_error
      endif
!
!  Check that the value of iunit makes sense, and find out whether the
!  corresponding file (if any) is meant to be formatted or unformatted.
      if(iunit.gt.0)then
         if(iunit.eq.5 .or. iunit.eq.6)then
            print*,'Wrong unit number in mbe_tdint_1 ',iunit
            call mbe_error
         endif
         inquire(unit=iunit,exist=fl)
         if(.not.fl)then
            print*,'mbe_tdint_1: the unit iunit does not exist ',iunit
            call mbe_error
         endif
         inquire(unit=iunit,opened=fl)
         if(.not.fl)then
            print*,'mbe_tdint_1: the unit iunit is not connected',iunit
            call mbe_error
         endif
         inquire(unit=iunit,form=file_format)
         if(nst*nst+2.le.999999)then
            write(fmt,'("(1x,",i6,"(1pe18.11,1x))")')nst*nst+2
         else
            print*,'mbe_tdint_1: The number of states is larger than'
            print*,'what this subroutine can cope with in regards to'
            print*,'preparing the format fmt. The code'
            print*,'should be updated in the unlikely event that'
            print*,'nst*nst + 2 > 999,999.'
            call mbe_error
         endif
      elseif(iunit.eq.-1)then
         call mbe_set_outputfiles
      elseif(iunit.lt.-1)then
         print*,'mbe_tdint_1: illegal value of iunit. ',iunit
         call mbe_error
      endif
!
      if(iunit_ampl.ge.0)then
         if(iunit_ampl.eq.5 .or. iunit_ampl.eq.6)then
            print*,'Wrong unit number in mbe_tdint_1 for iunit_ampl'
            print*,iunit_ampl
            call mbe_error
         endif
         inquire(unit=iunit_ampl,exist=fl)
         if(.not.fl)then
            print*,'mbe_tdint_1: the unit iunit_ampl does not exist ', &
                   iunit_ampl
            call mbe_error
         endif
         inquire(unit=iunit_ampl,opened=fl)
         if(.not.fl)then
            print*,'mbe_tdint_1: the unit iunit_ampl is not connected', &
                   iunit_ampl
            call mbe_error
         endif
         if(iunit.eq.iunit_ampl .and. iunit.ne.0 .and. iunit.ne.-1)then
            print*,'mbe_tdint_1: Conflict between iunit and iunit_ampl ', &
                   iunit,iunit_ampl
            call mbe_error
         endif
         if(fl_outputfiles)then
            do j = 1,nst
               do i = 1,nst
                  if(iunits_arr(i,j).eq.iunit_ampl)then
                     print*,'mbe_tdint_1: Conflict between iunit_ampl and'
                     print*,'iunits_arr for i,j = ',i,j
                     call mbe_error
                  endif
               enddo
            enddo
         endif
      endif
!
!  Check that the numerical integration over the velocity classes
!  has been initialised though a prior call to obe_set_Doppler, if the
!  calculation involves Doppler averaging. The logical variable
!  fl_Doppler_notinit is imported from the obe module.
      if(iDoppler.eq.1)then
         call mbe_set_Doppler
         if(fl_Doppler_notinit)then
            print*,'mbe_tdint_1: attempt at a calculation involving'
            print*,'a Doppler averaging before a call to obe_set_Doppler.'
            call mbe_error
         else
            fl_Doppler = .true.
            ivmax = n_v_values
         endif
      elseif(iDoppler.eq.0)then
         fl_Doppler = .false.
         ivmax = 1
      else
         print*,'mbe_tdint_1: illegal value of iDoppler ',iDoppler
         call mbe_error
      endif
!
!  Initialise rhsmat_0
      call mbe_set_rhsmat_0
!
!  Signal to re-initialize the arrays set by mbe_set_rhsmat (the logical
!  variable fl_init_ommats is global within the module), unless we continue
!  with the arrays initialized at a previous run.
      fl_init_ommats = .true.
      if(present(icont))then
         if(icont.eq.1)then
            if(fl_Doppler.neqv.fl_Doppler_old)then
               print*,'mbe_tdint_1: a continuation run with a different'
               print*,'Doppler option is not allowed.'
               call mbe_error
            endif
            fl_init_ommats = .false.
         endif
      endif
      fl_Doppler_old = fl_Doppler
!
!  Check that the complex field amplitudes have been provided (they must
!  be provided through a call to mbe_set_tdfields) and if they have been
!  transcribe them into the relevant arrays.
      if(fl_notinit_fields)then
         print*,'mbe_tdint_1 was called before a call to mbe_set_tdfields.'
         call mbe_error
      endif
!  Write the input fields into the relevant arrays.
      if(allocated(Field_p_re_vec))deallocate(Field_p_re_vec)
      if(allocated(Field_p_im_vec))deallocate(Field_p_im_vec)
      allocate(Field_p_re_vec(0:n_time_steps),                 &
               Field_p_im_vec(0:n_time_steps),                 &
               stat = istat)
      if(istat.ne.0)then
         print*,'mbe_tdint_1: Error at the allocation of the arrays (1).'
         call mbe_error
      endif
      do k = 0, n_time_steps
         Field_p_re_vec(k) = real(dreal(Field_p_vec(k)),kd)
         Field_p_im_vec(k) = real(dimag(Field_p_vec(k)),kd)
      enddo
!
!  Allocate the rest of the storage:
!
      allocate(rhovec(neqs,ivmax),       &
               rhs_mat(neqs,neqs),       &
               stat=istat)
      if(istat.gt.0)then
         print*,'mbe_tdint_1: Error at the allocation of the arrays (2).'
         call mbe_error
      endif
!
!  Prepare the spline interpolation used in the calculation of the density
!  matrix if necessary:
      if(fl_spline)then
         allocate(Fd_pre2(0:n_time_steps),  &
                  Fd_pim2(0:n_time_steps),  &
                  stat=istat)
         if(istat.gt.0)then
            print*,'mbe_tdint_1: Error at the allocation of the arrays (3).'
            call mbe_error
         endif
         call mbe_spline_prep(Field_p_re_vec,Fd_pre2)
         call mbe_spline_prep(Field_p_im_vec,Fd_pim2)
      endif
!
!  Allocation of storage used in case of calculations using the
!  rate equation approximation for the initial time. The arrays arhs and
!  ipiv_r are global within this module.
!
      if(istart.eq.1)then
         if(allocated(arhs))deallocate(arhs)
         if(allocated(ipiv_r))deallocate(ipiv_r)
         allocate(arhs(neqs,neqs),ipiv_r(neqs),stat=istat)
         if(istat.gt.0)then
            print*,'mbe_tdint_1: Error at the allocation of arhs or ipiv_r.'
            call mbe_error
         endif
      endif
!
!  Propagate the density matrix in time.
!  The initial condition is governed by the variable istart, which is passed to
!  mbe_propagate_1 and must be set to 0, 1, 2 or -1 by the calling program. If
!  istart is 2, the atoms are assumed to be initially in the steady-state they
!  evolve into if starting in the state with a non-zero initial population.
!  Otherwise, the atoms are assumed to be initially in the mixed state defined
!  by the population array popinit: if istart is 0, the coherences are zero,
!  whereas if istart is 1 the coherences are calculated in the rate equations
!  approximation. The initial density matrix is the initial content of
!  rhovec_av if istart = -1, and in this case popinit is not used if present.
!  Note that istart = -1 might not make any sense for a calculation with
!  Doppler averaging, as setting istart = -1 means that the initial density
!  matrix is the same for all the velocity classes.
      if(istart.eq.-1)then
         if(allocated(rhovec_0))deallocate(rhovec_0)
         allocate(rhovec_0(neqs),stat=istat)
         if(istat.gt.0)then
            print*,'mbe_tdint_1: Error at the allocation of rhovec_0.'
            call mbe_error
         endif
         rhovec_0 = rhovec_av
      endif
!
      do k = 0,n_time_steps
!
!  Loop over velocity (runs only if the calculation involves a Doppler average)
!
         if(fl_Doppler)rhovec_av = 0.0_kd
         do iv = 1,ivmax
            if(fl_Doppler)then
               current_vel = vmesh(iv)
            else
               current_vel = 0.0_kd ! The value of current_vel is relevant 
                                    ! only when fl_Doppler is .true.
            endif
!
            if(k.eq.0)then
               Field_p_re = Field_p_re_vec(0)
               Field_p_im = Field_p_im_vec(0)
               call mbe_set_rhsmat1(rhs_mat,0)
               if(istart.eq.-1)then
                  rhovec(:,iv) = rhovec_0
               else
                  rhovec(:,iv) = 0.0_kd
                  do i = 1,nst
                     call ldbl_pop_index(i,m)
                     rhovec(m,iv) = popinit(i)
                  enddo
!  Unless istart is 1 or 2, the coherences are initially zero.
                  if(istart.eq.2)then
!  Steady state
                     call ldbl_steadystate(rhovec(:,iv),2,rhs_mat)
                  elseif(istart.eq.1)then
!  Rate approximation
                     call mbe_calc_coher(rhs_mat,rhovec(:,iv),1)
                     call mbe_calc_coher(rhs_mat,rhovec(:,iv),0)
                  elseif(istart.ne.0)then
                     print*,'mbe_tdint_1: illegal value of istart ',istart
                     call mbe_error
                  endif
               endif
!
            else
!
               ktime = k - 1
               call mbe_intobe(imethod,t_mesh(k-1),t_mesh(k),         &
                                           nsbsteps(k),rhovec(:,iv))
!
            endif
!
!  Averaging over velocity classes if needed.
            if(fl_Doppler)then
               rhovec_av = rhovec_av + fMvweight(iv) * rhovec(:,iv)
            else
               rhovec_av = rhovec(:,iv)
            endif
!
         enddo
!
         if(iunit.gt.0 .or. iunit.eq.-1 .or. iunit_ampl.gt.0)then
            call mbe_output_td(k,t_mesh(k))
         endif
!
      enddo
!
!  Reset fl_twofields to its default value.
      fl_twofields = .true.
!
      deallocate(Field_p_re_vec,Field_p_im_vec)
      if(fl_spline)deallocate(Fd_pre2,Fd_pim2)
      deallocate(rhovec,rhs_mat)
!
      contains

!!!!!!!!!!!!!!!!!!!

      subroutine mbe_output_td(k,t)
!
!  Save the results on file. A value of 0 for the
!  input variable k is a signal to rewind the unit iunit before writing.
!
      implicit none
      real(kd), intent(in) :: t
      integer, intent(in) :: k
      integer :: i,iost,j,m,mre,mim
!
      if(iunit.gt.0)then
!
         if(k.eq.0)rewind iunit
         if(file_format.eq.'UNFORMATTED')then
            write(unit=iunit,iostat=iost)t,(rhovec_av(j),j=1,neqs)
         else
            write(iunit,fmt,iostat=iost)t,(rhovec_av(j),j=1,neqs)
         endif
         if(iost.ne.0)then
            print*,'mbe_output_td: write returns with iostat = ',iost
            call mbe_error
         endif
 
      elseif(iunit.eq.-1)then
 
         if(fl_outputfiles)then
            do j = 1,nst
               do i = 1,nst
                  if(iunits_arr(i,j).gt.0)then
                     if(i.eq.j)then
                        call ldbl_pop_index(i,m)
                        write(iunits_arr(i,j),1501)t,rhovec_av(m)
                     elseif(i.lt.j)then
                        call ldbl_coher_index(i,j,mre,mim)
                        write(iunits_arr(i,j),1502)t,rhovec_av(mre), &
                                                     rhovec_av(mim)
                     else
                        call ldbl_coher_index(j,i,mre,mim)
                        write(iunits_arr(i,j),1502)t,rhovec_av(mre), &
                                                    -rhovec_av(mim)
                     endif
                  endif
               enddo
            enddo
         else
            print*,'mbe_output_td: The value of iunit is -1 but the unit'
            print*,'numbers of the output files have not been provided.'
            call mbe_error
         endif
      endif
!
      if(iunit_ampl.gt.0)then
!
         call mbe_calc_flds(t)
         write(iunit_ampl,1601)t,Field_p_re,Field_p_im
!
      endif
!
 1501 format(1x,2(1pe12.5,2x))
 1502 format(1x,3(1pe12.5,2x))
 1601 format(1x,3(1pe12.5,2x))
!
      end subroutine mbe_output_td

!!!!!!!!!!!!!!!!!!!

      end subroutine mbe_tdint_1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end module mbe
