      module ldblstore
!
!  This module is used by ldbl_tdint to store the density matrix
!  calculated at each time step.
!
      use general_settings
!
      implicit none
      real(kd), dimension(:,:), allocatable, public  :: rho_allt
      real(kd), public :: factor
!
      end module ldblstore
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      module ldbl
!
!  The program units contained in this module integrate the Lindblad master
!  equation.
!
!  The parameters defining the system (components of the Hamiltonian matrix,
!  decay rates and T2 dephasing rates) must be passed to ldbl through a
!  subroutine... 
!
      use general_settings
!
      implicit none
! 
!  The following statement makes the contents of ldbl "invisible" to other
!  parts of the programs, with the exceptions of the variables and subprograms
!  declared below as public.
!
      private   
!
!------------------------------------------------------------------------------
!
!  The following subprograms contained in ldbl may be called from outside the
!  module:
!
      public ldbl_tdint, ldbl_steadystate, ldbl_mat2vec, ldbl_vec2mat, &
             ldbl_set_tol_dop853, ldbl_pop_index, ldbl_coher_index, &
             ldbl_set_rhsmat, ldbl_steadystate_prep, ldbl_steadystate_calc
!
!  All the other subprograms contained in ldbl can be accessed only from ldbl.
!
!------------------------------------------------------------------------------
!
!  The following variables will be "known" by all the subprograms
!  contained in this module but are not directly accessible from the outside.
!
      integer, parameter :: ndm = nst ! <-- number of dimensions
!
!
!  Logical variable used to check or control the job flow
      logical :: fl_indexes_not_init = .true.
!
!  The following variables are used to convert between storage formats
      integer, dimension (ndm,ndm), save :: kindex
      integer, dimension (ndm,ndm), save :: mindex
!
!  The following variables are used by ldbl_steadystate_calc and _prep.
      integer, save :: j_state = -1
      integer, save :: m_stdy = -1
!
!  Variables used to communicate with the dop853 ODE solver
      double precision, save :: rtol_dop853,atol_dop853
      logical, save :: fl_set_tol_dop853 = .false.
!
!
!------------------------------------------------------------------------------

      contains
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
      subroutine ldbl_error
!
!  Control may be passed to this subprogram from elsewhere within ldbl
!  following detection of an error.
!
!  At the moment ldbl_error simply stops the execution without saving files.
!  The user may want to replace this by a more sophisticated error handling
!  procedure.
!
      stop 'Error in the ldbl module...'
!
      end subroutine ldbl_error
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
      subroutine ldbl_set_indexes
!
!  Set the index arrays kindex and mindex used by ldbl_reformat_rhs_cmat,
!  ldbl_mat2vec, ldbl_vec2mat and possibly also elsewhere in the program.
!  These indexes are global variables, as is the flag fl_indexes_not_init.
!  (The latter, if .true., indicates that these index arrays are not
!  initialized; it is initially .true. and is reset to .false. by
!  ldbl_set_indexes.)
!
!  The contents of these index arrays depend only on ndm, the size of the space,
!  and not on any other parameters of the problem. These arrays need to be
!  set only once in any run.
!
      implicit none
      integer :: i,j,k,m
!
!  Set kindex: the complex density matrix can be represented either as
!  a ndm x ndm square matrix or as a vector of ndm**2 components. The index
!  kindex(i,j) is the index of the component of the vector corresponding to
!  the (i,j) component of the matrix.
      k = 0
      do j=1,ndm
         do i=1,ndm
            k = k + 1
            kindex(i,j) = k
         enddo
      enddo
!
!  Set mindex: The diagonal components of the density matrix (which are real)
!  and the real and imaginary parts of the components forming the rest of
!  its upper triangle can be arranged into a vector of ndm*ndm components. 
!  The m-th component of that vector is either rho(i,i) or real[rho(i,j)] with
!  i < j or imag[rho(i,j)] with i < j. For i = j, mindex(i,j) is the index of
!  the component of this vector which corresponds to rho(i,i). For i < j,
!  mindex(i,j) is the index of the component which corresponds to
!  real[rho(i,j)] and mindex(j,i) the index of the component which corresponds
!  to imag[rho(i,j)] (note the swap between real and imaginary parts if i > j).
!
      m = 0
      do j=1,ndm
         do i=1,j
            m = m + 1
            mindex(i,j) = m
            if(i.lt.j)then
               m = m + 1
               mindex(j,i) = m
            endif
         enddo
      enddo
!
!  Signal that the index arrays have been initialized:
      fl_indexes_not_init = .false.
!
      end subroutine ldbl_set_indexes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
      subroutine ldbl_set_rhsmat(rhs_mat,rhs_cmat)
!
!  Calculates the matrix which gives rhodot_vec when acting on rho_vec,
!  where rho_vec is the density matrix represented as a vector of ndm**2
!  components, and rhodot_vec is the vector containing the time derivatives
!  of the components of the density matrix.
!
!  The calculation makes use of the array of indexes kindex and mindex.
!  These arrays are initialized by ldbl_set_rhsmat if they aren't yet.
!
      implicit none
      interface
         subroutine ext_setsys(h_mat,decrates_mat,t2dephas_mat,ncoll)
         use general_settings
         complex(kd), dimension (nst,nst) :: h_mat
         real(kd), dimension (:,:,:) :: decrates_mat
         complex(kd), dimension (nst,nst) :: t2dephas_mat
         integer, intent(out) :: ncoll
         end subroutine ext_setsys
      end interface
      integer, parameter :: neqs = ndm*ndm
      real(kd), dimension(neqs,neqs), intent(out), optional :: rhs_mat
      complex(kd), dimension(neqs,neqs), intent(out), optional :: rhs_cmat
      complex(kd), dimension(ndm,ndm) :: h_mat
      real(kd), dimension(:,:,:), allocatable :: decrates_mat
      complex(kd), dimension(ndm,ndm) :: t2dephas_mat
      complex(kd), allocatable, dimension(:,:) :: rhs_cmat_local
      complex(kd), dimension(neqs) :: crho_vec,crhodot_vec
      real(kd), dimension(neqs) :: rho_vec
      complex(kd), dimension(ndm,ndm) :: rho,rhodot
      complex(kd), dimension(ndm,ndm) :: superop_mat
      integer, dimension(neqs) :: indx
      integer :: i,j,k,ncoll,istat
      logical :: flreal,flcomplex
!
!  Check the arguments are correct.
!
      if(present(rhs_mat))then
         flreal = .true.
      else
         flreal = .false.
      endif
      if(present(rhs_cmat))then
         flcomplex = .true.
      else
         flcomplex = .false.
      endif
      if(.not.flreal .and. .not.flcomplex)then
         print*,'ldbl_set_rhsmat: Missing argument(s).'
         call ldbl_error
      endif
!
!  To start, check whether the index arrays kindex and mindex have already
!  been set, and call the subprogram which sets them if thay haven't.
!  These arrays are global variables. Also, allocate the (possibly very
!  large) array rhs_cmat. This array is deallocated before exit from this
!  subroutine.
      if(fl_indexes_not_init) call ldbl_set_indexes
      if(.not.flcomplex)then
         allocate(rhs_cmat_local(neqs,neqs),stat=istat)
         if(istat.gt.0)then
            print*,'ldbl_set_rhsmat: Error at the allocation of '
            print*,'rhs_cmat_local.'
            call ldbl_error
         endif
      endif
!
!  The calculation proceeds in two steps: (1) Obtain the matrix representing
!  the right-hand sides in the description using complex coherences. This
!  gives a complex matrix, rhs_cmat. (2) Transform this matrix into the
!  matrix representing the right-hand in the description using real and
!  imaginary parts of the coherences (i,j) for j > i. The real part of this
!  second matrix is the matrix rhs_mat returned by this subroutine. This
!  second part of the calculation is done only if the optional argument
!  rhs_mat is specified in the calling list. 
!
!  First step. We calculate rhs_cmat column by column, each column
!  being obtained by calculating the derivative of rho for a rho which
!  has only one non-zero component. We loop over all the columns of rhs_cmat
!  in this way. Note how the do-loops are nested. This ensures that the order
!  in which the columns are stored is consistent with the order of the
!  components in the vector crhodot_vec.
!
      call ext_setsys(h_mat,decrates_mat,t2dephas_mat,ncoll)
!
      do j=1,ndm
         do i=1,ndm
            k = kindex(i,j)
!
!  Set the appropriate component of rho to 1, and the others to 0:
            rho = (0.0_kd,0.0_kd)
            rho(i,j) = (1.0_kd,0.0_kd)
!
!  Obtain the corresponding right-hand side of the optical Bloch equations
!  (the subroutine ldbl_lindblad is contained within ldbl_set_rhsmat and puts
!  the result into the square matrix rhodot):
            call ldbl_lindblad
!
!  Repackage the square array rhodot into a vector
            crhodot_vec = reshape(rhodot,shape=(/neqs/))
!
!  Store this as the k-th column of rhs_cmat:
            if(flcomplex)then
               rhs_cmat(:,k) = crhodot_vec
            else
               rhs_cmat_local(:,k) = crhodot_vec
            endif
!
!  End of the double loop
         enddo
      enddo
!
!  Skip the next step if it is not required
!
      if(.not.flreal)go to 100
!
!  If required, transform rhs_cmat into the right-hand matrix for calculations
!  on the real and imaginary parts separately. There are two methods, one
!  using the internal subroutine ldbl_reformat_rhs_cmat and one not using
!  this subroutine. The latter method has a simpler logic and makes it
!  possible to treat models in which the diagonal elements of the
!  Hamiltonian might have a non-zero imaginary part.
!
!  First method
!!!   if(flcomplex)then
!!!      call ldbl_reformat_rhs_cmat(rhs_cmat)
!!!      rhs_mat = real(rhs_cmat,kd)
!!!   else
!!!      call ldbl_reformat_rhs_cmat(rhs_cmat_local)
!!!      rhs_mat = real(rhs_cmat_local,kd)
!!!   endif
!
!  Second method
!
      do k = 1,neqs
         rho_vec = 0.0_kd
         rho_vec(k) = 1.0_kd
         call ldbl_vec2mat(rho_vec,rho)
         crho_vec = reshape(rho,shape=(/neqs/))
         if(flcomplex)then
            crhodot_vec = matmul(rhs_cmat,crho_vec)
         else
            crhodot_vec = matmul(rhs_cmat_local,crho_vec)
         endif
         rhodot = reshape(crhodot_vec,shape=(/ndm,ndm/))
         call ldbl_mat2vec(rhodot,rhs_mat(:,k))
      enddo
!
 100  continue
!
      if(allocated(rhs_cmat_local))deallocate(rhs_cmat_local)
!
      contains

!!!!!!!!!!!!!!!!!!

      subroutine ldbl_reformat_rhs_cmat(rhs_cmat)
!
!  Calculation for upper triangle (still complex).
!
!  This subroutine takes the matrix representing the right-hand sides of
!  the optical Bloch equations, when written in terms of a full density
!  matrix of real populations and complex coherences, and transforms it
!  into the matrix representing the right-hand sides when the equations are
!  written in terms of a set of real populations, real parts of the coherences
!  (i,j) for j > i, and imaginary parts of the coherences (i,j) for j > i.
!
!  In short, in both representations the ndm x ndm density matrix is written
!  in the form of a vector of ndm**2 components. If M is the matrix of
!  right-hand sides in the initial representation and T is the matrix
!  transforming the vector from the new representation to the initial
!  representation, then the transformed right-hand-sides matrix is
!  M' = T**(-1) M T, where T**(-1) is the inverse of T.  The calculation is
!  programmed as M' = transpose[transpose(M T) transpose(T**(-1))], as 
!  tranpose(T**(-1)) is almost identical to T and organising it in this
!  way avoids manipulating array rows (it is hoped that the tranpose function
!  is hightly optimized).
!
!  The calculation proceeds in four steps: Calculate Mstore =  M T;
!  transpose Mstore; repeat the first step changing for the transpose Mstore
!  and the transpose inverse T, and finally transpose the end result to get M'.
!
      implicit none
      integer, parameter :: neqs = ndm*ndm
      complex(kd), dimension(neqs) :: vec_storel,vec_storeu
      complex(kd), dimension(:,:), intent(inout) :: rhs_cmat
      complex(kd), dimension(neqs,neqs) :: rhs_store
      integer :: k,kl,ku,m,mrl,mim
      integer :: i,j
!
      do j=1,ndm
         do i=1,j
            if(i.eq.j)then
!  For i.eq.j, simply copy the k(i,i)-th column of the original matrix into
!  the m(i,i)-th column of the storage matrix 
               m = mindex(i,i)
               k = kindex(i,i)
               rhs_store(:,m) = rhs_cmat(:,k)
            else
!  For i.ne.j, put the sum of the k(i,j)-th and k(j,i)-th columns into the 
!  the m(i,j)-th column of the storage matrix. Also, put the difference of these
!  two columns, multiplied by the imaginary i, into the m(j,i)-th column of
!  the storage matrix.
               mrl = mindex(i,j)
               mim = mindex(j,i)
               ku = kindex(i,j)
               kl = kindex(j,i)
               vec_storeu = rhs_cmat(:,ku)
               vec_storel = rhs_cmat(:,kl)
               rhs_store(:,mrl) = vec_storeu + vec_storel
               rhs_store(:,mim) = (0.0_kd,1.0_kd)*(vec_storeu - vec_storel)
            endif
         enddo
      enddo
!
      rhs_store = transpose(rhs_store)
!
      do j=1,ndm
         do i=1,j
            if(i.eq.j)then
!  For i.eq.j, simply copy the k(i,i)-th column of the storage matrix into
!  the m(i,i)-th column of the new matrix
               m = mindex(i,i)
               k = kindex(i,i)
               rhs_cmat(:,m) = rhs_store(:,k)
            else
!  For i.ne.j, put the sum of the k(i,j)-th and k(j,i)-th columns (divided by 2)
!  into the the m(i,j)-th column of the new matrix. Also, put the difference
!  of these two columns, multiplied by -i/2, into the m(j,i)-th column of
!  the new matrix.
               mrl = mindex(i,j)
               mim = mindex(j,i)
               ku = kindex(i,j)
               kl = kindex(j,i)
               vec_storeu = rhs_store(:,ku)
               vec_storel = rhs_store(:,kl)
               rhs_cmat(:,mrl) = 0.5_kd*(vec_storeu + vec_storel)
               rhs_cmat(:,mim) = (0.0_kd,-0.5_kd)*(vec_storeu - vec_storel)
            endif
         enddo
      enddo
!
      rhs_cmat = transpose(rhs_cmat)
!
      end subroutine ldbl_reformat_rhs_cmat

!!!!!!!!!!!!!!!!!!
 
      subroutine ldbl_lindblad
! 
!  Given the matrices h_mat, decrates_mat and t2dephas_mat and given the
!  density matrix rho, calculate the time derivative of the latter, rhodot,
!  as given by the Lindblad master equation.
!
      implicit none
      integer :: ii,jj,kk
      complex(kd), dimension(ndm,ndm) :: coll_mat
!
!  Calculate -i [H,rho] (no 1/hbar factor here as H is expressed in terms of
!  angular frequencies, not energies).
!  We take into account the fact that at this point of the execution, all the
!  elements of the matrix rho are zero except the element (i,j) which is 1.
!  hence, the matrix H rho is a matrix whose j-th column is the same as the
!  i-th column of H and all other columns are zero (and the other way round
!  for the matrix rho H).
      rhodot = (0.0_kd,0.0_kd)
      rhodot(:,j) = rhodot(:,j) - (0.0_kd,1.0_kd)*h_mat(:,i)
      rhodot(i,:) = rhodot(i,:) + (0.0_kd,1.0_kd)*h_mat(j,:)
!
!  Contribution from the spontaneous decay:
!  Sum over all the decay channels, and for each of them add the
!  corresponding contribution to rhodot.
!  ii is the final state, jj is the initial state.
!
      if(ncoll.ge.1)then
        do kk=1,ncoll
           coll_mat = (1.0_kd,0.0_kd)*decrates_mat(:,:,kk)
           call ldbl_set_superopmat_coll(coll_mat,superop_mat,i,j)
           rhodot = rhodot + superop_mat
        enddo
      elseif(ncoll.eq.-1)then
         do ii=1,ndm
            do jj=1,ndm
               if(decrates_mat(ii,jj,1).gt.0.0_kd)then
                  call ldbl_set_superopmat_fast(decrates_mat(ii,jj,1), &
                     ii,jj,superop_mat,i,j)
                  rhodot = rhodot + superop_mat
               elseif(decrates_mat(ii,jj,1).lt.0.0_kd)then
                  print*,'ldbl_lindblad: Invalid decay rate'
                  print*,ii,jj,decrates_mat(ii,jj,1)
                  call ldbl_error
               endif
            enddo
         enddo
      else
         print*,'ldbl_lindblad: Illegal value of ncoll ',ncoll
         call ldbl_error
      endif
!
!  Add the contribution from T2 dephasing processes. Note that the * operation
!  in the following statement amounts to multiplying each element of rho by the
!  corresponding element of t2dephas_mat; this is not a matrix multiplication
!  in the usual mathematical sense. The contribution is subtracted in this
!  version of the module, as t2dphas_mat is contains the little gammas, not
!  minus the little gammas.
      rhodot = rhodot - t2dephas_mat*rho
!
      end subroutine ldbl_lindblad

!!!!!!!!!!!!!!!!!!

      end subroutine ldbl_set_rhsmat
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
      subroutine ldbl_set_superopmat_coll(coll_mat,superop_mat,i,j)
!
!  Construct the matrix representing the Lindblad superoperator arising
!  from spontaneous decay, given a collapse operator provided by the user.
!  In matrix form, this operator is meant to be sqrt(decrate) times coll_mat.
!  We use the fact that at this stage of the program the density operator, rho,
!  is here a matrix whose elements are all zero except the element (i,j) which
!  is 1. Hence, with C denoting the coll_mat matrix,
!  (i)   C rho Cdagger is the matrix obtained by taking the outer product of
!        the column vector C_(:,i) and the row vector Cdagger(j,:);
!  (ii)  rho Cdagger C is the matrix whose i-th row is the product of
!        the row vector Cdagger(j,:) and the matrix C (and all the other rows
!        of rho Cdagger C are zero);
!  (iii) Cdagger C rho is the matrix whose j-th column is the product of
!        the matrix Cdagger with the column vector C(:,i) (and all the other
!        columns of Cdagger C rho are zero).
!
      implicit none
      complex(kd), dimension (ndm,ndm), intent(out) :: superop_mat
      integer, intent(in) :: i,j
      complex(kd), dimension (ndm,ndm) :: coll_mat
      complex(kd), dimension (ndm,ndm) :: cdagger
      complex(kd), dimension (ndm,1) :: column
      complex(kd), dimension (1,ndm) :: row
!
      cdagger = conjg(transpose(coll_mat))
!
      column(:,1) = coll_mat(:,i)
      row(1,:) = cdagger(j,:)
      superop_mat = matmul(column,row)
      superop_mat(i,:) = superop_mat(i,:) - matmul(cdagger(j,:),coll_mat)/2
      superop_mat(:,j) = superop_mat(:,j) - matmul(cdagger,coll_mat(:,i))/2
!
      end subroutine ldbl_set_superopmat_coll
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
      subroutine ldbl_set_superopmat_fast(decrate,ii,jj,superop_mat,i,j)
!
!  Construct the matrix representing the Lindblad superoperator arising
!  from spontaneous decay, given the decay rates provided by the user.
!
!  Only for transitions between non-degenerate states.
!
      implicit none
      complex(kd), dimension (ndm,ndm), intent(out) :: superop_mat
      integer, intent(in) :: i,ii,j,jj
      real(kd), intent(in) :: decrate
!
      superop_mat = (0.0_kd,0.0_kd)
      if(j.eq.jj .and. i.eq.j)superop_mat(ii,ii) =     &
            superop_mat(ii,ii) + decrate
      if(i.eq.jj)superop_mat(i,j) = superop_mat(i,j) - decrate/2.0_kd
      if(j.eq.jj)superop_mat(i,j) = superop_mat(i,j) - decrate/2.0_kd
!
      end subroutine ldbl_set_superopmat_fast
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
      subroutine ldbl_tdint(irate,t1,t2,n_time_steps,imethod,rho_vec,  &
         istore,iunit,rhsmatin,iunits_arr,isteady)
!
!  Integrate the Lindblad equations over time from t1 to t2,
!  starting from the initial values passed to the subprogram through
!  the vector rho_vec, and return the result to the calling program
!  through the same vector (which is thus overwritten by ldbl_tdint).
!  This vector contains the populations and the real and imaginary parts
!  of the (i,j) coherences for j > i.
!
!  The variable irate must be 0 or 1 at entry. The rate equation approximation
!  is made if irate.eq.1 and is not made if irate.eq.0.
!
!  The steady state density matrix is calculated if the optional
!  argument isteady is present and equal to 1. The normal course of
!  action, i.e., integrating the optical Bloch equations for a finite
!  time interval, is followed if isteady = 0 or if this optional
!  argument is not present. The only legal values of isteady are 0 and 1.
!  If isteady = 1, then it is necessary that n_time_steps = 1 and
!  imethod = 3.
!
!  In the normal course of action, the interval [t1,t2] is divided into
!  n_time_steps steps of duration tstep, and the solution is propagated over
!  each step. The variable imethod determines the integration method. Its legal
!  values currently are:
!    imethod = 1:  use the fourth-order Runge-Kutta formula
!    imethod = 2:  use a user-supplied ODE solver
!    imethod = 3:  use an eigenvector-based solver (see ldbl_tdeig).
!    imethod = 4:  use Hairer, Norsett and Wanner's DOP853 ODE solver.
!    imethod = 5:  use a fifth-order Runge-Kutta formula
!
!  The variables istore and iunit control what ldbl_tdint writes out and
!  where and how this is written:
!
!  The variable istore can take three different values:
!     istore = 1:  ldbl_tdin allocates or re-allocates the array rho_allt
!                  of the module ldblstore and at each time step (including
!                  t1 and t2) writes the whole density matrix in the
!                  corresponding column of this array. This array has
!                  ndm**2 + 1 rows and n_time_steps + 1 columns. The first
!                  row (index 0) is used to store the value of t for the
!                  corresponding column. The first column (index 0) corresponds
!                  to t = t1, the last column (index n_time_steps) to t = t2.
!     istore = 0:  No use of rho_allt is made.
!     istore = -1: Same as istore = 1, but the density matrix calculated
!                  by ldbl_tdint at each time step, multiplied by the variable
!                  factor of the module ldblstore, is added to the appropriate
!                  column of the array rho_allt. This array is not
!                  allocated or de-allocated by ldbl_tdint when istore.eq.-1.
! 
!  If iunit .gt. 0, the whole density matrix is written out on file at
!  each integration step, sequentially, using iunit as the unit number.
!  (The corresponding file is first positioned at its beginning by ldbl_tdint.
!  This operation will thus overwrite whatever information this file initially
!  contained. Note that ldbl_tdint does *not* open or close that file.) The
!  value of iunit cannot be 5 or 6 (these unit numbers are reserved for the
!  standard input and standard output).
!
!  If iunit .eq. -1, the elements of the density matrix for which a
!  positive unit number is specified in the optional argument iunits_arr
!  are written in the corresponding files. These files are assumed to be
!  ready for being written on, and the unit numbers are assumed to be
!  suitable (this is not verified).
!
!  If iunit has any other value, no results are written out within this
!  subroutine (the density matrix at t = t2 is simply passed to the calling
!  program through the array rho_vec).
!  
!  The matrix representing the right-hand sides can be passed to ldbl_tdint
!  as its last argument. This argument is optional: if this matrix is not 
!  passed to ldbl_tdint, then ldbl_tdint will obtain it directly from
!  ldbl_set_rhsmat.
!
!  Programming note: As this subprogram makes use of optional arguments, it
!  will require an interface block in each calling program if taken out of the
!  module.
!  
      use ldblstore
      implicit none
      integer, parameter :: neqs = ndm*ndm
      real(kd), dimension (neqs), intent(inout) :: rho_vec
      real(kd), dimension (neqs,neqs), optional, intent(in) :: rhsmatin
      integer, dimension (ndm,ndm), optional, intent(in) :: iunits_arr
      real(kd), intent(in) :: t1,t2
      integer, intent(in), optional :: isteady
      integer, intent(in) :: imethod,irate,istore,iunit,n_time_steps
      real(kd), dimension (neqs,neqs) :: rhs_mat
      real(kd), dimension (:,:), allocatable, save :: arhs
      real(kd), dimension (:,:), allocatable, save :: rhs_mat_rate
      real(kd), dimension (:), allocatable :: pop_vec
      real(kd), dimension (:), allocatable :: rhovec_temp,rhovecdot_temp
      real(kd), dimension(:), allocatable :: rhodot,rhok,rhokdot1,rhokdot2
      integer, dimension (:), allocatable, save :: ipiv
      real(kd) :: t,tstep
      integer :: i,istat,m,n,jsteady
      logical :: fl_cont,fl
      logical, save :: fl_init_rate = .true.
      character(len=11) :: file_format
!!!!  external ext_ode
!
!  Check that the values of the arguments are consistent if the optional
!  argument istead_in is present.
!
      if(present(isteady))then
         if(isteady.eq.1)then
            if(imethod.ne.3 .or. n_time_steps.ne.1)then  
               print*,'ldbl_tdint: imethod must be 3 and n_time_steps'
               print*,'must be 1 for a calculation with isteady = 1.'
               print*,imethod,n_time_steps
               call ldbl_error
            endif
            jsteady = 1
         elseif(isteady.eq.0)then
            jsteady = 0
         else
            print*,'ldbl_tdint: the optional argument isteady must be either'
            print*,'0 or 1.',isteady
            call ldbl_error
         endif
      else
         jsteady = 0
      endif
!
!  Prepare the rate equation calculation if irate.eq.1: allocate
!  the array used by ldbl_calc_coher (which is contained within
!  this subprogram).
!
      if(irate.eq.1)then
         if(allocated(arhs))deallocate(arhs)
         if(allocated(rhs_mat_rate))deallocate(rhs_mat_rate)
         if(allocated(pop_vec))deallocate(pop_vec)
         if(allocated(rhovec_temp))deallocate(rhovec_temp)
         if(allocated(rhovecdot_temp))deallocate(rhovecdot_temp)
         if(allocated(ipiv))deallocate(ipiv)
         allocate(arhs(neqs,neqs),rhs_mat_rate(ndm,ndm),pop_vec(ndm), &
            rhovec_temp(neqs),rhovecdot_temp(neqs),ipiv(neqs),stat=istat)
         if(istat.gt.0)then
            print*,'ldbl_tdint: Error at the allocation of arhs etc.'
            call ldbl_error
         endif
      elseif(irate.ne.0)then
         print*,'ldbl_tdint: Illegal value of irate ',irate
         call ldbl_error
      endif
!
!  Prepare the arrays used by rk4, dop853 or rk5 (if needed)
!
      if(imethod.eq.1 .or. imethod.eq.4 .or. imethod.eq.5)then
         if(irate.eq.1)then
            n = ndm
         else
            n = neqs
         endif
      endif
      if(imethod.eq.1)then
         call ldbl_rk4(rho_vec,t,tstep,n) ! n > 0 means initialization only
      elseif(imethod.eq.4)then
         call ldbl_dop853(rho_vec,t,tstep,n) ! n > 0 means initialization only
      elseif(imethod.eq.5)then
         call ldbl_rk5(rho_vec,t,tstep,n) ! n > 0 means initialization only
      endif
!
!  Check that the value of iunit makes sense, and find out whether the
!  corresponding file (if any) is meant to be formatted or unformatted.
!
      if(iunit.gt.0)then
         if(iunit.eq.5 .or. iunit.eq.6)then
            print*,'Wrong unit number in ldbl_tdint ',iunit
            call ldbl_error
         endif
         inquire(unit=iunit,exist=fl)
         if(.not.fl)then
            print*,'ldbl_tdint: the unit iunit does not exist ',iunit
            call ldbl_error
         endif
         inquire(unit=iunit,opened=fl)
         if(.not.fl)then
            print*,'ldbl_tdint: the unit iunit is not connected',iunit
            call ldbl_error
         endif
         inquire(unit=iunit,form=file_format)
      elseif(iunit.eq.-1)then
         if(.not.present(iunits_arr))then
            print*,'ldbl_tdint: The array containing the units numbers'
            print*,'of the output files is not present although the'
            print*,'subroutine has been called with iunit = -1.'
            call ldbl_error
         endif
      endif
!
!  If the module ldblstore is to be used, and depending on whether istore.eq.1
!  or -1, (re-)allocate the array rho_allt of that module (which destroys
!  its previous content) or check that this array is already allocated and
!  has the necessary size.
!
      if(istore.eq.1)then
        if(allocated(rho_allt))deallocate(rho_allt)
        allocate(rho_allt(0:neqs,0:n_time_steps),stat=istat)
        if(istat.gt.0)then
           print*,'ldbl_tdint: Error at the allocation of rho_allt.'
           call ldbl_error
        endif
      elseif(istore.eq.-1)then
         if(.not.allocated(rho_allt))then
            print*,'ldbl_tdint is called with istore.eq.-1 although the array'
            print*,'rho_allt is not allocated.'
            call ldbl_error
         elseif(lbound(rho_allt,1).ne.0 .or. &
                lbound(rho_allt,2).ne.0)then
            print*,'ldbl_tdint: Wrong lower bound for rho_allt'
            print*,lbound(rho_allt,1),lbound(rho_allt,2)
            call ldbl_error
         elseif(ubound(rho_allt,1).ne.neqs .or. &
                ubound(rho_allt,2).ne.n_time_steps)then
            print*,'ldbl_tdint: Wrong upper bound for rho_allt'
            print*,ubound(rho_allt,1),ubound(rho_allt,2)
            call ldbl_error
         endif
      elseif(istore.ne.0)then
         print*,'ldbl_tdint: illegal value for istore ',istore
         call ldbl_error
      endif
!
!  Obtain the matrix used to calculate the right-hand sides, unless this
!  matrix was passed to ldbl_tdint through its list of arguments.
      if(present(rhsmatin))then
         rhs_mat = rhsmatin
      else
         call ldbl_set_rhsmat(rhs_mat)
      endif
!
!  Prepare the relevant arrays if the calculation is to be done
!  within the rate equation approximation:
      if(irate.eq.1)call ldbl_rate_prep
!
!  Calculate the time step for the required number of steps: 
      if(n_time_steps.lt.1)then
         print*,'ldbl_tdint has been called with a non-sensical value'
         print*,'of n_time_steps. ',n_time_steps
         call ldbl_error
      endif
      tstep = (t2-t1)/real(n_time_steps,kd)
!
!  Now integrate the system, starting with the populations and coherences
!  passed to ldbl_tdint through the array rho_vec. At each time step, rho_vec
!  is overwritten with the propagated solution and t is updated to t + tstep.
!  The density matrix is written on file if iunit > 0. (The file is first
!  rewound and the initial density matrix is also written.) A value of .true.
!  for the logical variable fl_cont is a signal to ext_ode that the calculation
!  is merely a continuation of the previous call, while .false. signals that
!  initialization may be required.
!
!  In the rate equation approximation, only the populations are propagated
!  in time. However, the coherences are calculated at each time step,
!  so that the results written on file are for a complete density matrix.
      t = t1
      if(irate.eq.1)then
         do i = 1,ndm
            call ldbl_pop_index(i,m)
            pop_vec(i) = rho_vec(m)
         enddo
         call ldbl_calc_coher(pop_vec,rho_vec)
      endif
!
      call ldbl_tdint_save(0)
!
      do n = 1,n_time_steps
         if(n.eq.1)then
            fl_cont = .false.
         else
            fl_cont = .true.
         endif
         if(imethod.eq.1)then
            if(irate.eq.1)then
               call ldbl_rk4(pop_vec,t,tstep,0) ! Propagate the populations.
               call ldbl_calc_coher(pop_vec,rho_vec)
            else
               call ldbl_rk4(rho_vec,t,tstep,0) ! Propagate everything.
            endif
         elseif(imethod.eq.4)then
            if(irate.eq.1)then
               call ldbl_dop853(pop_vec,t,tstep,0) ! Propagate the populations.
               call ldbl_calc_coher(pop_vec,rho_vec)
            else
               call ldbl_dop853(rho_vec,t,tstep,0) ! Propagate everything.
            endif
         elseif(imethod.eq.5)then
            if(irate.eq.1)then
               call ldbl_rk5(pop_vec,t,tstep,0) ! Propagate the populations.
               call ldbl_calc_coher(pop_vec,rho_vec)
            else
               call ldbl_rk5(rho_vec,t,tstep,0) ! Propagate everything.
            endif
         elseif(imethod.eq.2)then
!!!!        call ext_ode(rhs_mat,rho_vec,t,tstep,fl_cont)
            print*,'ldbl_tdint: Option 2 for imethod is currently disabled.'
            call ldbl_error
         elseif(imethod.eq.3)then
            if(irate.eq.1)then
               ! Propagate the populations:
               call ldbl_tdeig_rate(rhs_mat_rate,pop_vec,t,tstep,fl_cont,    &
                                                                    jsteady)
               call ldbl_calc_coher(pop_vec,rho_vec)
            else
               ! Propagate everything:
               call ldbl_tdeig(rhs_mat,rho_vec,t,tstep,fl_cont,jsteady)
            endif
         else
            print*,'Illegal value of imethod in ldbl_tdint: ',imethod
            call ldbl_error
         endif
!
         call ldbl_tdint_save(n)
!
      enddo
!
      contains
 
!!!!!!!!!!!!!!!!!!!

      subroutine ldbl_tdint_save(n)
!
!  Save the results on file or in the array rho_allt (from ldblstore).
!  If writing the whole density matrix on file (iunit.gt.0 in ldbl_tdint),
!  a value of 0 for the input variable n is a signal to rewind the unit iunit
!  before writing.
!
!  Also possible is to write selected elements of the density matrix,
!  as determined by the array iunits_arr (an optional argument of
!  ldbl_tdint). This option is activated by a value of 0 of iunit.
!
      use ldblstore
      implicit none
      integer, intent(in) :: n
      integer :: i,iost,j,m,mre,mim
!
      if(iunit.gt.0)then
         if(n.eq.0)rewind iunit
         if(file_format.eq.'UNFORMATTED')then
            write(unit=iunit,iostat=iost)t,(rho_vec(j),j=1,neqs)
         else
            write(iunit,*,iostat=iost)t,(rho_vec(j),j=1,neqs)
         endif
         if(iost.ne.0)then
            print*,'ldbl_tdint_save: write returns with iostat = ',iost
            call ldbl_error
         endif
      elseif(iunit.eq.-1)then
         do j = 1,ndm
            do i = 1,ndm
               if(iunits_arr(i,j).gt.0)then
                  if(i.eq.j)then
                     call ldbl_pop_index(i,m)
                     write(iunits_arr(i,j),1501)t,rho_vec(m)
                  elseif(i.lt.j)then
                     call ldbl_coher_index(i,j,mre,mim)
                     write(iunits_arr(i,j),1502)t,rho_vec(mre),rho_vec(mim)
                  else
                     call ldbl_coher_index(j,i,mre,mim)
                     write(iunits_arr(i,j),1502)t,rho_vec(mre),-rho_vec(mim)
                  endif
               endif
            enddo
         enddo
      endif
      if(istore.eq.1)then
         rho_allt(0,n) = t 
         rho_allt(1:neqs,n) =  rho_vec(1:neqs)
      elseif(istore.eq.-1)then
         rho_allt(0,n) = t 
         rho_allt(1:neqs,n) =  &
            rho_allt(1:neqs,n) + factor * rho_vec(1:neqs)
      endif
!
 1501 format(1x,2(1pe12.5,2x))
 1502 format(1x,3(1pe12.5,2x))
!
      end subroutine ldbl_tdint_save
 
!!!!!!!!!!!!!!!!!!!

      subroutine ldbl_rhs(t,rho,rhodot)
!
!  Calculates the right-hand side of the master equation and returns
!  the result in rhodot.
!
!  At the moment, t (time) is not used by ldbl_rhs...
!
      implicit none
      real(kd), dimension (:), intent(in) :: rho
      real(kd), dimension (:), intent(out) :: rhodot
      real(kd), intent(in) :: t
      integer :: i,m
!
!  The matrix rhs_mat comes from the subroutine ldbl_tdint, in which
!  the subroutine ldbl_rhs is contained. The matrix rhs_mat_rate is
!  constructed by ldbl_rate_prep (which should be called prior to the
!  first call to ldbl_rhs).
!
      if(irate.eq.1)then
!  Calculation for the system of ndm rate equations.
!  First, check that ldbl_rate_prep has been called:
         if(fl_init_rate)then
            print*,'ldbl_rhs was called before ldbl_rate_prep...'
            call ldbl_error
         endif
         rhodot = matmul(rhs_mat_rate,rho)
      else
!  Calculation for the system of ndm**2 optical Bloch equations.
         rhodot = matmul(rhs_mat,rho)
      endif
!
      end subroutine ldbl_rhs

!!!!!!!!!!!!!!!!!!!

      subroutine ldbl_rate_prep
!
!  Prepare the calculation in the rate equation approximation. I.e., the
!  matrix representing the optical Bloch equation is factorized and
!  the smaller matrix representing the rate equations is calculated.
!
!  This subroutine makes use of the arrays rhs_mat, arhs and rhs_mat_rate,
!  which it inherits from ldbl_tdint.
!
!  This subroutine must be called prior any use of the arrays arhs (used
!  to calculate the coherences) and rhs_mat_rate (used to integrate the
!  rate equations). It must be called again each time rhs_mat is modified.
!  
      implicit none
      integer, parameter :: neqs = ndm*ndm
      real(kd), dimension(:), allocatable :: rho,rhodot
      integer :: i,info,istat,j,lda,m,numb,nrhs
      external dgetrf,sgetrf
!
!  We first calculate and factorize the matrix representing the system, arhs.
!
      arhs = rhs_mat
!
!  The populations should not be recalculated (which amounts to putting 1
!  on the corresponding diagonal elements of the matrix of the system, and
!  0 in the corresponding off-diagonal elements.
      do i = 1,ndm
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
         call sgetrf(numb,numb,arhs,lda,ipiv,info)
      elseif(kd.eq.kind(1.d0))then
         call dgetrf(numb,numb,arhs,lda,ipiv,info)
      else
         print*,'Logic error detected in ldbl_rate_prep'
         call ldbl_error
      endif
!
!  info.ne.0 : error condition
!  info.eq.0 : normal return from s/dgetrf
!
      if(info.ne.0)then
         print*,'ldbl_rate_prep: s/dgetrf returns with info = ',info
         call ldbl_error
      endif
!
!  Reset fl_init_rate (this is done at this stage, rather than at the
!  very end of ldbl_rate_prep, so that ldbl_calc_coher can be used a few
!  lines below).
      fl_init_rate = .false.
!
!  Construct the matrix representing the right-hand sides in the
!  rate equation approximation, column by column
!
      if(allocated(rho))deallocate(rho)
      if(allocated(rhodot))deallocate(rhodot)
      allocate(rho(ndm),rhodot(ndm),stat=istat)
      if(istat.gt.0)then
         print*,'ldbl_rate_prep: Error at the allocation of rho/rhodot'
         call ldbl_error
      endif
!
      do j = 1,ndm
         rho = 0.0_kd
         rho(j) = 1.0_kd
         call ldbl_calc_coher(rho,rhovec_temp)
         rhovecdot_temp = matmul(rhs_mat,rhovec_temp)
         do i = 1,ndm
            call ldbl_pop_index(i,m)
            rhodot(i) = rhovecdot_temp(m)
         enddo
         rhs_mat_rate(:,j) = rhodot
      enddo
!
      deallocate(rho,rhodot)
!
!  That's all...
!
      end subroutine ldbl_rate_prep

!!!!!!!!!!!!!!!!!!!

      subroutine ldbl_calc_coher(pop_vec,rho_vec)
!
!  Given the populations, calculate the coherences in the rate equation
!  approximation. The populations are passed to ldbl_calc_coher through
!  the array pop_vec, and the populations and coherences are passed back
!  to the calling program through the array rho_vec.
!
!  This subroutine makes use of the array arhs, which it inherits from
!  ldbl_tdint. This array is initialized by ldbl_rate_prep, which needs
!  to be called before the first call to ldbl_calc_coher.
!
      implicit none
      integer, parameter :: neqs = ndm*ndm
      real(kd), dimension(:), intent(in) :: pop_vec
      real(kd), dimension(neqs), intent(out) :: rho_vec
      integer :: i,info,lda,ldb,m,numb,nrhs
      external dgetrs,sgetrs
!
!  First check that the value fl_init_rate.eqv.false
!
      if(fl_init_rate)then
         print*,'ldbl_calc_coher was called before ldbl_rate_prep...'
         call ldbl_error
      endif
!
!  Set all the elements of rho_vec to zero, except the populations.
      rho_vec = 0.0_kd
      do i = 1,ndm
         call ldbl_pop_index(i,m)
         rho_vec(m) = pop_vec(i)
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
      numb = neqs
      nrhs = 1
      lda = neqs
      ldb = neqs
!
      if(kd.eq.kind(1.0))then
         call sgetrs('N',numb,nrhs,arhs,lda,ipiv,rho_vec,ldb,info)
      elseif(kd.eq.kind(1.d0))then
         call dgetrs('N',numb,nrhs,arhs,lda,ipiv,rho_vec,ldb,info)
      else
         print*,'Logic error detected in ldbl_calc_coher'
         call ldbl_error
      endif
!
!  info.ne.0 : error condition
!  info.eq.0 : normal return from s/dgetrs.
!
      if(info.ne.0)then
         print*,'ldbl_calc_coher: s/dgetrs returns with info = ',info
         call ldbl_error
      endif
!
!  At this stage, rho_vec contains the correct coherences (within the rate
!  equation approximation) but not the correct populations. Replace the
!  (wrong) populations by those passed to ldbl_calc_coher:
      do i = 1,ndm
         call ldbl_pop_index(i,m)
         rho_vec(m) = pop_vec(i)
      enddo
!
!  That's all...
!
      end subroutine ldbl_calc_coher

!!!!!!!!!!!!!!!!!!!
 
      subroutine ldbl_rk4(rho,tstart,tstep,n)
!
!  Propagates the array rho over one time step of duration tstep using
!  the fourth order Runge-Kutta method. The derivatives at time t are
!  calculated in the subroutine ldbl_rhs(t,rho,rhodot).
!
!  rho: At entry the array rho at t = tstart. At return, the array rho
!       at t = tstart + tstep
!  tstart: The initial time. At return, the value of tstart is reset to
!          the final time, tstart + tstep.
!  n: If n > 0, allocate the arrays using n as their dimension, and do
!     nothing else. Otherwise, integrates the system of equations.
!
      implicit none
      real(kd), dimension(:), intent(inout) :: rho
      real(kd), intent(inout) :: tstart
      real(kd), intent(in) :: tstep
      integer, intent(in) :: n
      real(kd) :: half_tstep,thalfstep
      real(kd), dimension(:), allocatable, save :: rhodot,rhok,rhokdot1, &
         rhokdot2
      integer :: istat
!
      if(n.gt.0)then
!
!  Only allocate the arrays used by rk4 (if needed)
!
         if(allocated(rhodot))deallocate(rhodot)
         if(allocated(rhok))deallocate(rhok)
         if(allocated(rhokdot1))deallocate(rhokdot1)
         if(allocated(rhokdot2))deallocate(rhokdot2)
         allocate(rhodot(n),rhok(n),rhokdot1(n),rhokdot2(n),stat=istat)
         if(istat.gt.0)then
            print*,'ldbl_rk4: Error at the allocation of rhodot etc'
            call ldbl_error
         endif
!
      else
!
!  Integrate the equations. Start by calculating the derivatives at the
!  initial value of t:
         call ldbl_rhs(tstart,rho,rhodot)
!
         half_tstep = 0.5_kd*tstep
         thalfstep = tstart + half_tstep
!
         rhok = rho + half_tstep*rhodot           ! rhok is now rho + k_1 / 2
         call ldbl_rhs(thalfstep,rhok,rhokdot1)
         rhok = rho + half_tstep*rhokdot1         ! rhok is now rho + k_2 / 2
         call ldbl_rhs(thalfstep,rhok,rhokdot2)
         rhodot = rhodot + 2.0_kd*(rhokdot1+rhokdot2)
         rhok = rho + tstep*rhokdot2              ! rhok is now rho + k_3
         call ldbl_rhs(tstart+tstep,rhok,rhokdot1)
         rhodot = rhodot + rhokdot1
!
         rho = rho + (tstep/6.0_kd)*rhodot
         tstart = tstart + tstep
!
      endif
!
      end subroutine ldbl_rk4

!!!!!!!!!!!!!!!!!!!

      subroutine ldbl_rk5(y,x,h,n)
!
!  Fifth order Runge Kutta ("Butcher's formula")
!
!  n: If n > 0, allocate the arrays using n as their dimension, and do
!     nothing else. Otherwise, integrates the system of equations.
!
      implicit none
      real(kd), dimension(:), intent(inout) :: y
      real(kd), intent(inout) :: x
      real(kd), intent(in) :: h
      integer, intent(in) :: n
      real(kd), dimension(:), allocatable, save :: yst,dy1,dy2,dy3,dy4,dy5,dy6
      real(kd) :: hh
      integer :: istat
!
      if(n.gt.0)then
!
!  Only allocate the arrays used by rk5 (if needed)
!
         if(allocated(yst))deallocate(yst)
         if(allocated(dy1))deallocate(dy1)
         if(allocated(dy2))deallocate(dy2)
         if(allocated(dy3))deallocate(dy3)
         if(allocated(dy4))deallocate(dy4)
         if(allocated(dy5))deallocate(dy5)
         if(allocated(dy6))deallocate(dy6)
         allocate(yst(n),dy1(n),dy2(n),dy3(n),dy4(n),dy5(n),dy6(n), &
            stat=istat)
         if(istat.gt.0)then
            print*,'ldbl_rk5: Error at the allocation of yst etc'
            call ldbl_error
         endif
!
      else
!
         call ldbl_rhs(x,y,dy1)
         yst = y + h*dy1/4.0_kd
         call ldbl_rhs(x+h/4.0_kd,yst,dy2)
         yst = y + h*(dy1+dy2)/8.0_kd
         call ldbl_rhs(x+h/4.0_kd,yst,dy3)
         yst = y + h*(-dy2/2.0_kd+dy3)
         call ldbl_rhs(x+h/2.0_kd,yst,dy4)
         yst = y + h*(3.0_kd*dy1+9.0_kd*dy4)/16.0_kd
         call ldbl_rhs(x+3.0_kd*h/4.0_kd,yst,dy5)
         yst = y + h*(-3.0_kd*dy1+2.0_kd*dy2+12.0_kd*dy3-12.0_kd*dy4+ &
                                                   8.0_kd*dy5)/7.0_kd
         call ldbl_rhs(x+h,yst,dy6)
         y = y + h*(7.0_kd*dy1+32.0_kd*dy3+12.0_kd*dy4+32.0_kd*dy5+ &
                                                   7.0_kd*dy6)/90.0_kd
         x = x + h
!
      endif
!
      end subroutine ldbl_rk5

!!!!!!!!!!!!!!!!!!!

      subroutine ldbl_dop853(rho,tstart,tstep,n)
!
!  Propagates the array rho over one time step of duration tstep using
!  the dop853 subroutines, which ldbl_dop853 is merely an interface
!  with.
!
!  The arguments of ldbl_dop853 are the same as those of ldbl_rk4. See
!  ldbl_rk4 for further information.
!
      implicit none
      real(kd), dimension(:), intent(inout) :: rho
      real(kd), intent(inout) :: tstart
      real(kd), intent(in) :: tstep
      integer, intent(in) :: n
      real(kd), dimension(:), allocatable :: work
      integer, dimension(:), allocatable :: iwork
      real(kd) , dimension(1) :: rtol,atol,rpar
      real(kd) :: uround
      integer, dimension(1) :: ipar
      integer :: itol,iout,idid,istat,lwork,liwork,nrdens
      integer :: neqs = -1
      external fcn_dummy,solout_dummy
      save
!
!  Check that the module is programmed in double precision, to avoid
!  a mismatch of variable type with dop853
      if(kd.ne.kind(1.d0))then
         print*,'ldbl_dop853 can be used only if the rest of the module'
         print*,'is programmed in double precision.'
         call ldbl_error
      endif
!
      if(n.gt.0)then
!
!  Only initialise the calculation and allocate the arrays used
!  by dop853.
!
         neqs = n
         nrdens = 0  ! No dense output required.
         lwork = 11*n + 8*nrdens + 21
         liwork = nrdens + 21
         if(allocated(work))deallocate(work)
         if(allocated(iwork))deallocate(iwork)
         allocate(work(lwork),iwork(liwork),stat=istat)
         if(istat.gt.0)then
            print*,'ldbl_dop853: Error at the allocation of the arrays.'
            call ldbl_error
         endif
!  fl_set_tol_dop853, rtol_top853 and atol_top853 are global parameters.
         if(.not.fl_set_tol_dop853)then
            print*,'ldbl_dop853: the tolerance parameters rtol and atol'
            print*,'should have been initialized by a call to'
            print*,'obe_set_tol_dop853 before the first call to'
            print*,'ldbl_dop853.'
            call ldbl_error
         endif
         rtol = rtol_dop853
         atol = atol_dop853
         itol = 0
         iout = 0
         uround = spacing(1.d0)
!
      else
!
         if(neqs.le.0)then
            print*,'ldbl_dop853: unexpected value of neqs.'
            print*,'Perhaps an error of logic.'
            call ldbl_error
         endif
         work = 0.d0
         work(1) = uround
         iwork = 0
!  The following arguments are not used by dop853 in the present
!  context, but need to be included in the calling list of consistency:
!      fcn_dummy, solout_dummy, rpar, ipar.
!  fcn_dummy and solout_dummy are the names of dummy external subroutines
!  which do not do anything. These two subroutines do not belong to
!  this module but have been appended at the end of this file.
!
!  Note that at the return from dop853, tstart will have been updated
!  to tstart+tstep, unless there would be an error.
!
         call dop853(neqs,fcn_dummy,tstart,rho,tstart+tstep,     &
             rtol,atol,itol,solout_dummy,iout,work,lwork,        &
             iwork,liwork,rpar,ipar,idid)
!
         if(idid.lt.0)then
            print*,'ldbl_dop853: dop853 returns with an error condition.'
            print*,idid
            call ldbl_error
         endif
! 
      endif
!
      end subroutine ldbl_dop853

!!!!!!!!!!!!!!!!!!!

! *********************************************************************
!  The following is a Fortran 90 translation of Hairer, Norsett and Wanner's
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
!! 
!!    double precision :: x,y,xend,rtol,atol,work,rpar
!!    double precision :: uround,safe,fac1,fac2,beta,hmax,h
      real(kd) :: x,y,xend,rtol,atol,work,rpar
      real(kd) :: uround,safe,fac1,fac2,beta,hmax,h
!!
      integer :: n,itol,iout,lwork,iwork,liwork,ipar,idid
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
      real(kd), parameter ::                 &
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
      real(kd), parameter ::                 &
     &  b1 =   5.42937341165687622380535766363D-2,   &
     &  b6 =   4.45031289275240888144113950566D0,    &
     &  b7 =   1.89151789931450038304281599044D0,    &
     &  b8 =  -5.8012039600105847814672114227D0,     &
     &  b9 =   3.1116436695781989440891606237D-1,    &
     &  b10 = -1.52160949662516078556178806805D-1,   &
     &  b11 =  2.01365400804030348374776537501D-1,   &
     &  b12 =  4.47106157277725905176885569043D-2
      real(kd), parameter ::                  &
     &  bhh1 = 0.244094488188976377952755905512D+00,  &
     &  bhh2 = 0.733846688281611857341361741547D+00,  &
     &  bhh3 = 0.220588235294117647058823529412D-01   
      real(kd), parameter ::                 &
     &  er1 =  0.1312004499419488073250102996D-01,   &
     &  er6 = -0.1225156446376204440720569753D+01,   &
     &  er7 = -0.4957589496572501915214079952D+00,   &
     &  er8 =  0.1664377182454986536961530415D+01,   &
     &  er9 = -0.3503288487499736816886487290D+00,   &
     &  er10 =  0.3341791187130174790297318841D+00,   &
     &  er11 =  0.8192320648511571246570742613D-01,   &
     &  er12 = -0.2235530786388629525884427845D-01
      real(kd), parameter ::                   &
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
      real(kd), parameter ::                   &
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
      real(kd), parameter ::                  &
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
      real(kd), parameter ::                  &
     &  a141 =  5.61675022830479523392909219681D-2,   &
     &  a147 =  2.53500210216624811088794765333D-1,   &
     &  a148 = -2.46239037470802489917441475441D-1,   &
     &  a149 = -1.24191423263816360469010140626D-1,   &
     &  a1410 =  1.5329179827876569731206322685D-1,   &
     &  a1411 =  8.20105229563468988491666602057D-3,  &
     &  a1412 =  7.56789766054569976138603589584D-3,  &
     &  a1413 = -8.298D-3
      real(kd), parameter ::                   &
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
      real(kd), parameter ::                   &
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
      real(kd), parameter ::                   &
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
      real(kd), parameter ::                   &
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
      real(kd), parameter ::                   &
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
      real(kd) :: x,xend,work,rpar
      integer :: n,itol,iout,lwork,iwork,liwork,ipar,idid,nrd
      real(kd) :: uround,safe,fac1,fac2,beta,hmax,h,cont
      integer :: nfcn,nstep,naccpt,nrejct,iprint,nmax,meth,  &
                 nstiff,nrdens,iek1,iek2,iek3,iek4,iek5,iek6,&
                 iek7,iek8,iek9,iek10,iey1,ieco,istore,icomp
      real(kd) :: Y(N),Y1(N),K1(N),K2(N),K3(N),K4(N),K5(N),K6(N)
      real(kd) :: K7(N),K8(N),K9(N),K10(N),ATOL(*),RTOL(*)     
      DIMENSION CONT(8*NRD),ICOMP(NRD),RPAR(*),IPAR(*)
      real(kd) :: facold,expo1,facc1,facc2,posneg,atoli, &
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
      POSNEG=SIGN(1.D0,dble(XEND-X))
! --- INITIAL PREPARATIONS   
      ATOLI=ATOL(1)
      RTOLI=RTOL(1)    
      LAST=.FALSE. 
      HLAMB=0.D0
      IASTI=0
      CALL ldbl_rhs(X,Y,K1) !!! CALL FCN(N,X,Y,K1,RPAR,IPAR)
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
!!! ldble module:
          print*,'DP86CO was called for a non-zero value of IOUT'
          call ldbl_error
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
         CALL ldbl_rhs(X,Y,K1) !!! CALL FCN(N,X,Y,K1,RPAR,IPAR)
      END IF
      DO 22 I=1,N 
      Y1(I)=Y(I)+H*A21*K1(I)  
  22  continue
      CALL ldbl_rhs(X+C2*H,Y1,K2) !!! CALL FCN(N,X+C2*H,Y1,K2,RPAR,IPAR)
      DO 23 I=1,N 
      Y1(I)=Y(I)+H*(A31*K1(I)+A32*K2(I))  
  23  continue
      CALL ldbl_rhs(X+C3*H,Y1,K3) !!! CALL FCN(N,X+C3*H,Y1,K3,RPAR,IPAR)
      DO 24 I=1,N 
      Y1(I)=Y(I)+H*(A41*K1(I)+A43*K3(I))  
  24  continue
      CALL ldbl_rhs(X+C4*H,Y1,K4) !!! CALL FCN(N,X+C4*H,Y1,K4,RPAR,IPAR)
      DO 25 I=1,N 
      Y1(I)=Y(I)+H*(A51*K1(I)+A53*K3(I)+A54*K4(I))
  25  continue
      CALL ldbl_rhs(X+C5*H,Y1,K5) !!! CALL FCN(N,X+C5*H,Y1,K5,RPAR,IPAR)
      DO 26 I=1,N 
      Y1(I)=Y(I)+H*(A61*K1(I)+A64*K4(I)+A65*K5(I))
  26  continue
      CALL ldbl_rhs(X+C6*H,Y1,K6) !!! CALL FCN(N,X+C6*H,Y1,K6,RPAR,IPAR)
      DO 27 I=1,N 
      Y1(I)=Y(I)+H*(A71*K1(I)+A74*K4(I)+A75*K5(I)+A76*K6(I))
  27  continue
      CALL ldbl_rhs(X+C7*H,Y1,K7) !!! CALL FCN(N,X+C7*H,Y1,K7,RPAR,IPAR)
      DO 28 I=1,N 
      Y1(I)=Y(I)+H*(A81*K1(I)+A84*K4(I)+A85*K5(I)+A86*K6(I)+A87*K7(I))  
  28  continue
      CALL ldbl_rhs(X+C8*H,Y1,K8) !!! CALL FCN(N,X+C8*H,Y1,K8,RPAR,IPAR)
      DO 29 I=1,N 
      Y1(I)=Y(I)+H*(A91*K1(I)+A94*K4(I)+A95*K5(I)+A96*K6(I)+A97*K7(I)  &
     &   +A98*K8(I))
  29  continue
      CALL ldbl_rhs(X+C9*H,Y1,K9) !!! CALL FCN(N,X+C9*H,Y1,K9,RPAR,IPAR)
      DO 30 I=1,N 
      Y1(I)=Y(I)+H*(A101*K1(I)+A104*K4(I)+A105*K5(I)+A106*K6(I)        &
     &   +A107*K7(I)+A108*K8(I)+A109*K9(I))
  30  continue
      CALL ldbl_rhs(X+C10*H,Y1,K10) !!! CALL FCN(N,X+C10*H,Y1,K10,RPAR,IPAR)
      DO 31 I=1,N 
      Y1(I)=Y(I)+H*(A111*K1(I)+A114*K4(I)+A115*K5(I)+A116*K6(I)        &
     &   +A117*K7(I)+A118*K8(I)+A119*K9(I)+A1110*K10(I))
  31  continue
      CALL ldbl_rhs(X+C11*H,Y1,K2) !!! CALL FCN(N,X+C11*H,Y1,K2,RPAR,IPAR)
      XPH=X+H
      DO 32 I=1,N 
      Y1(I)=Y(I)+H*(A121*K1(I)+A124*K4(I)+A125*K5(I)+A126*K6(I)        &
     &   +A127*K7(I)+A128*K8(I)+A129*K9(I)+A1210*K10(I)+A1211*K2(I))
  32  continue
      CALL ldbl_rhs(XPH,Y1,K3) !!! CALL FCN(N,XPH,Y1,K3,RPAR,IPAR)
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
         CALL ldbl_rhs(XPH,K5,K4) !!! CALL FCN(N,XPH,K5,K4,RPAR,IPAR)
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
            CALL ldbl_rhs(X+C14*H,Y1,K10)!!!CALL FCN(N,X+C14*H,Y1,K10,RPAR,IPAR)
            DO 52 I=1,N 
               Y1(I)=Y(I)+H*(A151*K1(I)+A156*K6(I)+A157*K7(I)         &
     &            +A158*K8(I)+A1511*K2(I)+A1512*K3(I)+A1513*K4(I)     &
     &            +A1514*K10(I))
  52        continue
            CALL ldbl_rhs(X+C15*H,Y1,K2)!!!CALL FCN(N,X+C15*H,Y1,K2,RPAR,IPAR)
            DO 53 I=1,N 
               Y1(I)=Y(I)+H*(A161*K1(I)+A166*K6(I)+A167*K7(I)         &
     &            +A168*K8(I)+A169*K9(I)+A1613*K4(I)+A1614*K10(I)     &
     &            +A1615*K2(I))                                        
  53        continue
            CALL ldbl_rhs(X+C16*H,Y1,K3)!!!CALL FCN(N,X+C16*H,Y1,K3,RPAR,IPAR)
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
      real(kd) :: x,y,xend,posneg,f0,f1,y1,hmax,atol,    &
                          rtol,rpar
      integer :: n,iord,itol,ipar
      DIMENSION Y(N),Y1(N),F0(N),F1(N),ATOL(*),RTOL(*)
      DIMENSION RPAR(*),IPAR(*)
      real(kd) :: dnf,dny,atoli,rtoli,sk,h,der2,h1,hinit,der12
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
      CALL ldbl_rhs(X+H,Y1,F1) !!! CALL FCN(N,X+H,Y1,F1,RPAR,IPAR) 
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
      real(kd) :: x,con,xold,h
      integer :: ii,icomp,nd
      DIMENSION CON(8*ND),ICOMP(ND)
      real(kd) :: s,s1,conpar,contd8
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

      end subroutine ldbl_tdint

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine ldbl_set_tol_dop853(rtol_in,atol_in,flag_in)
!
!  This subroutine is used by obe to communicate the values of the
!  tolerance parameters rtol and atol to ldbl, for use by the
!  dop853 ode solver. It is a simple mailbox.
!
      implicit none
      double precision, intent(in) :: rtol_in,atol_in
      logical, intent(in) :: flag_in
!
!  rtol_dop853, atol_dop853 and fl_set_tol_dop853 are global variables.
      rtol_dop853 = rtol_in
      atol_dop853 = atol_in
      fl_set_tol_dop853 = flag_in
!
!  Variables used to communicate with the dop853 ODE solver
      end subroutine ldbl_set_tol_dop853

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine ldbl_tdeig(rhs_mat,rho,t,tstep,fl_cont,isteady)
!
!  Propagate the vector rho over one time step of duration tstep by representing
!  this arrays as a linear combination of the right-hand eigenvectors of the
!  matrix rhs_mat, or obtain the steady state deriving from the vector rho
!  passed to ldbl_tdeig.
!
!  The calculation of the steady state is identical to the calculation of
!  the time-dependent density matrix, apart that for the steady state only
!  the eigenvectors belonging to a zero eigenvalue are used. For the purpose
!  of the calculation, a zero eigenvalue is defined as one whose real and
!  imaginary parts are both smaller than the parameter eps defined below,
!  in absolute value.
!
!  rho: At entry, the array rho at time t (used only if fl_cont.eqv..false.).
!       At return, the arrays rho at time = t + tstep, or the steady state
!       density matrix if isteady.eq.1.
!  t: The initial time. At return, the value of t is reset to the final time,
!     t + tstep, unless isteady.eq.1.
!  fl_cont: A .false. value triggers the calculation of the eigenvalues and
!           eigenvectors of rhs_mat (which is overwritten) before the array
!           rho is calculated. A .true. value means that rho is calculated
!           using the coefficient of the expansion calculated at the last
!           call to ldbl_tdeig with fl_cont .eqv. .false.
!  isteady: If isteady.eq.1, the subroutine returns the steady state density
!           matrix deriving from the input density matrix, and t is not
!           updated. This option is incompatible with fl_cont.eqv..true.
!           If isteady.eq.0, the calculation is done as normal. Any other
!           value of isteady is illegal. 

      implicit none
!
      real(kd) :: eps
!
      integer, parameter :: neqs=ndm*ndm
      real(kd), dimension(neqs,neqs), intent(inout) :: rhs_mat
      real(kd), dimension(neqs), intent(inout) :: rho
      real(kd), intent(inout) :: t
      real(kd), intent(in) :: tstep
      integer, intent(in) :: isteady
      logical, intent(in) :: fl_cont
!
      character*1 :: jobvl,jobvr
      integer :: info,lda,ldvl,ldvr,lwork,n
      real(kd), dimension(neqs) :: wr,wi
      real(kd), dimension(neqs,neqs), save :: vl,vr
      real(kd), dimension(4*neqs) :: work
!
      complex(kd), dimension(neqs,neqs), save :: vl_cplx,vr_cplx
      complex(kd), dimension(neqs), save :: cm_0,w
      complex(kd), dimension(neqs), save :: rho_cplx
      real(kd), save :: t_0
      real(kd) :: deltat
      integer :: i,j,m
      external dgeev,sgeev
!
!  Define the parameter eps, according to the precision
!
      if(kd.eq.kind(1.0))then
         eps = 1.e-3_kd
      else
         eps = 1.e-8_kd
      endif
!
!  Check that the input parameters are legal:
      if(isteady.eq.1)then
         if(fl_cont)then
            print*,'Wrong combination of isteady and fl_cont in ldbl_tdeig.'
            call ldbl_error
         endif
      elseif(isteady.ne.0)then
         print*,'Wrong value of isteady in ldbl_tdeig: ',isteady
         call ldbl_error
      endif
!
      if(.not.fl_cont)then
!
         jobvl = 'V'
         jobvr = 'V'
         n = neqs
         lda = neqs
         ldvl = neqs
         ldvr = neqs
         lwork = 4*neqs
         if(kd.eq.kind(1.0))then
            call sgeev(jobvl,jobvr,n,rhs_mat,lda,wr,wi,vl,ldvl,vr,ldvr,   &
                       work,lwork,info)
         elseif(kd.eq.kind(1.d0))then
            call dgeev(jobvl,jobvr,n,rhs_mat,lda,wr,wi,vl,ldvl,vr,ldvr,   &
                       work,lwork,info)
         else
            print*,'Logic error in ldbl_tdeig...',kd,kind(1.0),kind(1.d0)
            call ldbl_error
         endif
         if(info.ne.0)then
            print*,'s/dgeev (from Lapack) returned with info = ',info
            print*,'in ldbl_tdeig.'
            call ldbl_error
         endif

!  Calculate the coefficients of the expansion, and save the initial value of t.
         rho_cplx = rho
         m = 1
         do 
            if(m.gt.neqs)exit
            if(wi(m).eq.0.0_kd)then
               w(m) = cmplx(wr(m),0.0_kd,kd)
               vl_cplx(:,m) = cmplx(vl(:,m),0.0_kd,kd)
               vr_cplx(:,m) = cmplx(vr(:,m),0.0_kd,kd)
               if(abs(dot_product(vl_cplx(:,m),vr_cplx(:,m))).eq.0.0_kd)then
                  print*,'ldbl_tdeig: zero eigenvector...'
                  call ldbl_error
               endif
               cm_0(m) = dot_product(vl_cplx(:,m),rho_cplx)/  &
                         dot_product(vl_cplx(:,m),vr_cplx(:,m))
               m = m + 1
            else
               if(m.ge.neqs)then
                  print*,'Logic Error in ldbl_tdeig', m
                  call ldbl_error
               endif
               w(m) = cmplx(wr(m),wi(m),kd)
               vl_cplx(:,m) = cmplx(vl(:,m),vl(:,m+1),kd)
               vr_cplx(:,m) = cmplx(vr(:,m),vr(:,m+1),kd)
               if(abs(dot_product(vl_cplx(:,m),vr_cplx(:,m))).eq.0.0_kd)then
                  print*,'ldbl_tdeig: zero eigenvector...'
                  call ldbl_error
               endif
               cm_0(m) = dot_product(vl_cplx(:,m),rho_cplx)/  &
                         dot_product(vl_cplx(:,m),vr_cplx(:,m))
               m = m + 1
               w(m) = conjg(w(m-1))
               vl_cplx(:,m) = conjg(vl_cplx(:,m-1))
               vr_cplx(:,m) = conjg(vr_cplx(:,m-1))
               cm_0(m) = conjg(cm_0(m-1))
               m = m + 1
            endif
         enddo
         t_0 = t
!
      endif
!
!  Calculate rho at time t + tstep by making a linear combination of the
!  right-hand eigenvectors of rhs_mat. The coefficient of the m-th eigenvector
!  at time t + tstep is the coefficient at time t_0 (saved in cm_0(m))
!  multiplied by exp[w(m)*(t+tstep-t_0)], where w(m) is the
!  corresponding eigenvalue (this eigenvalue can be complex).
!
!  Note that this coefficient, the exponential and the eigenvector are all
!  real if the eigenvalue is real, while if it is complex then the imaginary
!  parts cancel between terms corresponding to eigenvalues that are complex
!  conjugate of each other, with the consequence that the sum is always real
!  even though the eigenvalues and eigenvectors may be complex.
!
!  First, check that the initial vector is suitably represented by
!  a linear combination of right eigenvectors
      rho = 0.0_kd
      do m = 1,neqs
         rho = rho + real(cm_0(m)*vr_cplx(:,m),kd)
      enddo
!  rho_cplx still has the initial vector...
      do m = 1,neqs
         if(abs(rho(m)-real(rho_cplx(m),kd)).gt.eps)then
            print*,'ldbl_tdeig: Warning of a possible issue with the basis of'
            print*,'right eigenvectors. ',m,abs(rho(m)-real(rho_cplx(m),kd))
         endif
      enddo
!
!  Normal calculation if isteady.eq.0, otherwise steady-state calculation
      if(isteady.eq.0)then
         deltat = t + tstep - t_0
         rho = 0.0_kd
         do m = 1,neqs
            rho = rho + real(exp(w(m)*deltat)*cm_0(m)*vr_cplx(:,m),kd)
         enddo
         t = t + tstep
      else
         rho = 0.0_kd
         do m = 1,neqs
            if(abs(real(w(m),kd)).lt.eps .and. abs(aimag(w(m))).lt.eps)then
               rho = rho + real(cm_0(m)*vr_cplx(:,m),kd)
            elseif(abs(real(w(m),kd)).lt.eps .and. abs(cm_0(m)).gt.eps)then
               print*,'ldbl_tdeig finds oscillating coherences.'
               print*,'The density matrix does not converge to'
               print*,'a constant matrix in the long time limit.'
               call ldbl_error
            endif
         enddo
      endif
!
      end subroutine ldbl_tdeig

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine ldbl_mat2vec(rho_mat,rho_vec)
!
!  Transform the matrix rho_mat into a vector of real and
!  imaginary parts of the components forming its upper triangle
!  (it is assumed that the diagonal elements are real).
!  There is no check that the lower triangle of rho_mat is consistent
!  with the upper triangle, i.e., that rho_mat is Hermitian.
!  The elements of rho_mat below the diagonal are ignored altogether.
!
      implicit none
      complex(kd), dimension (ndm,ndm), intent(in) :: rho_mat
      real(kd), dimension (ndm*ndm), intent(out) :: rho_vec
      integer :: i,j,m,mrl,mim
!
!  To start, check whether the index arrays kindex and mindex have already
!  been set, and call the subprogram which sets them if thay haven't.
!  These arrays are global variables.
      if(fl_indexes_not_init) call ldbl_set_indexes
!
      do j=1,ndm
         do i=1,j
            if(i.eq.j)then
               m = mindex(i,i)
               rho_vec(m) = real(rho_mat(i,i),kd)
            else
               mrl = mindex(i,j)
               mim = mindex(j,i)
               rho_vec(mrl) = real(rho_mat(i,j),kd)
               rho_vec(mim) = aimag(rho_mat(i,j))
            endif
         enddo
      enddo
!
      end subroutine ldbl_mat2vec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
      subroutine ldbl_vec2mat(rho_vec,rho_mat)
!
!  Transform the vector rho_vec containing the real and
!  imaginary parts of the components forming the upper triangle
!  of a matrix into this matrix itself, assuming that this matrix
!  is Hermitian.
!
      implicit none
      complex(kd), dimension (ndm,ndm), intent(out) :: rho_mat
      real(kd), dimension (ndm*ndm), intent(in) :: rho_vec
      integer :: i,j,m,mrl,mim
!
!  To start, check whether the index arrays kindex and mindex have already
!  been set, and call the subprogram which sets them if thay haven't.
!  These arrays are global variables.
      if(fl_indexes_not_init) call ldbl_set_indexes
!
      do j=1,ndm
         do i=1,j
            if(i.eq.j)then
               m = mindex(i,i)
               rho_mat(i,i) = cmplx(rho_vec(m),0.0_kd,kd)
            else
               mrl = mindex(i,j)
               mim = mindex(j,i)
               rho_mat(i,j) = cmplx(rho_vec(mrl),rho_vec(mim),kd)
               rho_mat(j,i) = conjg(rho_mat(i,j))
            endif
         enddo
      enddo 
!
      end subroutine ldbl_vec2mat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
      subroutine ldbl_pop_index(i,m)
!
!  Return the index of the component corresponding to the population rho(i,i)
!  in the 1D array representing the density matrix in vector form.
!
      integer, intent(in) :: i
      integer, intent(out) :: m
!
!  To start, check whether the index arrays kindex and mindex have already
!  been set, and call the subprogram which sets them if thay haven't.
!  These arrays are global variables.
      if(fl_indexes_not_init) call ldbl_set_indexes
!
!  Check that the value of i makes sense and exits if it doesn't.
      if(i.lt.1 .or. i.gt.ndm)then
         print*,'ldbl_pop_index was called with an illegal value of i: ',i
         call ldbl_error
      endif
!
      m = mindex(i,i)
!
      end subroutine ldbl_pop_index

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
      subroutine ldbl_coher_index(i,j,mrl,mim)
!
!  Return the indexes of the components corresponding to the real and the
!  imaginary parts of the coherence rho(i,j) (j > i) in the 1D array
!  representing the density matrix in vector form.
!
      integer, intent(in) :: i,j
      integer, intent(out) :: mrl,mim
!
!  To start, check whether the index arrays kindex and mindex have already
!  been set, and call the subprogram which sets them if thay haven't.
!  These arrays are global variables.
      if(fl_indexes_not_init) call ldbl_set_indexes
!
!  Check that the value of i makes sense and exits if it doesn't.
      if(i.lt.1 .or. i.gt.ndm .or. j.lt.1 .or. j.gt.ndm .or. j.le.i)then
         print*,'ldbl_coher_index was called with an illegal value of i'
         print*,'or of j: ',i,j
         call ldbl_error
      endif
!
      mrl = mindex(i,j)
      mim = mindex(j,i)
!
      end subroutine ldbl_coher_index

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine ldbl_steadystate(rho_st_vec,ioption,rhsmatin)
!
!  Calculate the steady-state density matrix.
!
!  rho_st_vec: If ioption.eq.0 or 2, this arrays should contain the
!  initial density matrix on entry. The content this array has on entry
!  is not used if iption.eq.1. On return, this array contains the steady-state
!  density matrix.
!
!  There is a choice between two methods of calculation. If ioption.eq.1, 
!  the (neqs x neqs) matrix representing the r.h.s. of the Lindblad 
!  equation is transformed into a (neqs-1 x neqs-1) square matrix a_mat and
!  a vector b of (neqs-1) components. The steady state is then found
!  by solving the linear system a_mat x = b. The latter is done by the
!  subroutine dgesv of the Lapack library, within the subroutine
!  ldbl_steadystate_calc. This method will be unsuccesful in the case of
!  multiple steady states.
!
!  In the second method (selected if ioption.eq.2), the initial density matrix
!  is projected on the zero-eigenvalue eigenvectors of the matrix representing
!  the right-hand sides of the Lindblad equations. This approach is valid
!  even in the case of the possibility of multiple steady states.
!
!  If called with ioption.eq.0, ldbl_steadtstate will first try the first
!  method, and then the second if the calculation cannot be done using the
!  first method.

!  The matrix representing the right-hand side of the system can be passed to
!  this subprogram as its last argument. This argument is optional: if this
!  array is not passed, then ldbl_steadystate will obtain it directly from
!  ldbl_set_rhsmat.
!
!  Programming note: As this subprogram makes use of an optional argument, it
!  will require an interface block in each calling program if taken out of the
!  module.
!
      implicit none
      integer, parameter :: neqs = ndm*ndm
      real(kd), dimension (neqs,neqs), optional, intent(in) :: rhsmatin
      real(kd), dimension (neqs), intent(inout) :: rho_st_vec
      integer, intent(in) :: ioption
      real(kd), dimension (neqs,neqs) :: rhs_mat
      real(kd), dimension (neqs-1,neqs-1) :: a_mat
      real(kd), dimension (neqs-1) :: b
      real(kd) :: t,tstep
      integer :: ifail,isteady
      logical :: fl_cont
!
!  Check the value of ioption:
      if(ioption.ne.0 .and. ioption.ne.1 .and. ioption.ne.2)then
         print*,'ldbl_steadystate: Illegal value of ioption at entry ',ioption
         call ldbl_error
      endif
!
!  Obtain the matrix used to calculate the right-hand sides, unless this
!  matrix was passed to ldbl_steadystate through its list of arguments.
      if(present(rhsmatin))then
         rhs_mat = rhsmatin
      else
         call ldbl_set_rhsmat(rhs_mat)
      endif
!
!  Set ifail to 1, to avoid the possibility of an error in case of a change
!  to the logic of this subprogram or ldbl_steadystate_calc (the
!  latter will set ifail to 0 in case of successful completion, and normally
!  sets ifail to 1 in case of unsuccessful completion).
      ifail = 1
!
!  First method.
!
      if(ioption.eq.0 .or. ioption.eq.1)then
!
!  Calculate the arrays a_mat and b:
         call ldbl_steadystate_prep(rhs_mat,a_mat,b)
!  Solve the system and rearrange the result:
         call ldbl_steadystate_calc(a_mat,b,rho_st_vec,ifail)
!  Error condition:
         if(ifail.eq.1)then
            if(ioption.eq.1)then
               print*,'ldbl_steadystate: Unsuccessful calculation.'
               call ldbl_error
            endif
            print*,'ldbl_steadystate: Method 1 failed. Method 2 is used.'
         endif
      endif
!
!  Second method. The calculation is done by ldbl_tdeig (option istead.eq.1,
!  the values of t and tstep do not matter).
!  New: We use ldbl_tdint instead of ldbl_tdeig...
!
      if((ioption.eq.0 .and. ifail.eq.1) .or. ioption.eq.2)then
         t = 0.0_kd
         tstep = 0.0_kd
         fl_cont = .false.
         isteady = 1
!! tdeig is now replaced by tdint.
!!       call ldbl_tdeig(rhs_mat,rho_st_vec,t,tstep,fl_cont,isteady)
         call ldbl_tdint(1,t,t+tstep,1,3,rho_st_vec,                   &
                                     0,0,rhsmatin=rhs_mat,isteady=1)
      endif
!
      end subroutine ldbl_steadystate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine ldbl_steadystate_calc(a_mat,b,rho_st_vec_r,ifail)
!
!  Calculate the steady-state density matrix, given the arrays a_mat and b
!  (which this subroutine modifies). The result is returned through
!  the array rho_st_vec_r. If the subroutine cannot find the steady state
!  density matrix, a warning message is issued on the standard output,
!  the variable ifail is set to 1 on exit and no change is made to the array
!  rho_st_vec_r. (Failure is normally due to the existence of several 
!  steady states, rather than just one, or to the absence of any steady state.
!  Either of these two cases results in the subroutine dgesv's returning
!  with an error condition.) The variable ifail is set to 0 on exit if the
!  calculation is succesful.
!
      implicit none
      integer, parameter :: neqs = ndm*ndm
      real(kd), dimension (neqs-1,neqs-1), intent(inout) :: a_mat 
      real(kd), dimension (neqs-1), intent(inout) :: b
      real(kd), dimension (neqs), intent(out) :: rho_st_vec_r
      integer, intent(out) :: ifail
      integer :: i,m_stdy_local
      integer :: numb,nrhs,lda,ldb,info
      integer, dimension (neqs-1) :: ipv
      external dgesv,sgesv
!
      if(j_state.le.0)then
         print*,'ldbl_steadystate_calc was called before the first call'
         print*,'to ldbl_steadystate_prep...'
         call ldbl_error
      endif
!
      call ldbl_pop_index(j_state,m_stdy_local)
      if(m_stdy_local.ne.m_stdy)then
         print*,'ldbl_steadystate_calc: the value of m_stdy_local is'
         print*,'inconsistent with the value set by ldbl_steadystate_prep.'
         print*,m_stdy_local,m_stdy
         call ldbl_error
      endif
!
!  Variables passed to s/dgesv (solution of system of linear equations,
!  from lapack):
!    numb: the number of equations in the linear system to be solved
!    nrhs: number of right-hand sides (here one only)
!    a_mat: the coefficients of the system (this array is overwritten
!       by s/dgesv)
!    lda: leading dimension of the array a_mat
!    ipv: an integer vector (overwritten by s/dgesv)
!    b: the column vector(s) forming the right-hand side(s) of the system
!       (here there is only one right-hand side; this array is overwritten
!       by s/dgesv)
!    ldb: the leading dimension of the array b 
!    info: an integer variable whose value is set to zero by is/dgesv in case
!       of successful completion or to another value in case of problem
!       (see online information about s/dgesv for details)
!
      numb = neqs - 1
      nrhs = 1
      lda = neqs - 1
      ldb = neqs - 1
!
      if(kd.eq.kind(1.0))then
         call sgesv(numb,nrhs,a_mat,lda,ipv,b,ldb,info)
      elseif(kd.eq.kind(1.d0))then
         call dgesv(numb,nrhs,a_mat,lda,ipv,b,ldb,info)
      else
         print*,'Logic error detected in ldbl_steadystate_calc'
         call ldbl_error
      endif
!
!  info.ne.0 : error condition
!  info.eq.0 : normal return from s/dgesv.
!
      if(info.ne.0)then
!
!  Keep going, but issue an error message and set ifail to 1.
         print*,'ldbl_steadystate_calc: s/dgesv returns with info = ',info
         ifail = 1
!
      else
!
!  Construct the complex steady-state density matrix from the elements of b,
!  ending by calculating the state 1 population. The variable ifail is then
!  set to 0 to signal a successful completion.
!
         if(m_stdy.eq.1)then
            rho_st_vec_r(2:neqs) = b(1:neqs-1) 
         elseif(m_stdy.eq.neqs)then
            rho_st_vec_r(1:neqs-1) = b(1:neqs-1) 
         else
            rho_st_vec_r(1:m_stdy-1) = b(1:m_stdy-1) 
            rho_st_vec_r(m_stdy+1:neqs) = b(m_stdy:neqs-1) 
         endif
         rho_st_vec_r(m_stdy) = 1.0_kd
         do i = 1,ndm
            if(i.ne.j_state)rho_st_vec_r(m_stdy) = rho_st_vec_r(m_stdy) - &
               rho_st_vec_r(mindex(i,i))
         enddo
         ifail = 0
!
      endif
!
      end subroutine ldbl_steadystate_calc
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
      subroutine ldbl_steadystate_prep(rhs_mat,a_mat,b,jst)
!
!  Calculate the arrays a_mat and b used in calculations of the steady state
!  density matrix.
!
!  The matrix representing the right-hand sides is passed to this subprogram
!  as its first argument.
!
!  The (1,1) population is set equal to 1, unless the optional jst
!  argument is specified in which case the (jst,jst) population is
!  set equal to 1.
!
      implicit none
      integer, parameter :: neqs = ndm*ndm
      real(kd), dimension (neqs,neqs), intent(in) :: rhs_mat
      real(kd), dimension (neqs-1,neqs-1), intent(out) :: a_mat
      real(kd), dimension (neqs-1), intent(out) :: b
      integer, intent(in), optional :: jst
      integer :: i,j,m
!
!  Transform rhs_mat into a (neqs-1 x neqs-1) square matrix a_mat and
!  a vector b of (neqs-1) components, so that the steady state is found
!  by solving the linear system Ax = b.
!
      if(present(jst))then
         if(jst.lt.1 .or. jst.gt.ndm)then
            print*,'ldbl_steadystate_prep: the optional argument jst'
            print*,'has a non-sensical value. ',jst
            call ldbl_error
         endif
         j_state = jst
      else
         j_state = 1
      endif
!
      call ldbl_pop_index(j_state,m_stdy)
!
!  We set rho_{j_state,j_state} = 1, hence the corresponding row and the
!  coressponding column of rhs_mat go and the column vector b is the negative
!  of the latter...
!
      if(m_stdy.eq.1)then
         b(1:neqs-1) = -rhs_mat(2:neqs,1)
         a_mat(1:neqs-1,1:neqs-1) = rhs_mat(2:neqs,2:neqs)
      elseif(m_stdy.eq.neqs)then
         b(1:neqs-1) = -rhs_mat(1:neqs-1,neqs)
         a_mat(1:neqs-1,1:neqs-1) = rhs_mat(1:neqs-1,1:neqs-1)
      else
         b(1:m_stdy-1) = -rhs_mat(1:m_stdy-1,m_stdy)
         b(m_stdy:neqs-1) = -rhs_mat(m_stdy+1:neqs,m_stdy)
         a_mat(1:m_stdy-1,1:m_stdy-1) =  &
                                      rhs_mat(1:m_stdy-1,1:m_stdy-1)
         a_mat(m_stdy:neqs-1,1:m_stdy-1) =  &
                                      rhs_mat(m_stdy+1:neqs,1:m_stdy-1)
         a_mat(1:m_stdy-1,m_stdy:neqs-1) =  &
                                      rhs_mat(1:m_stdy-1,m_stdy+1:neqs)
         a_mat(m_stdy:neqs-1,m_stdy:neqs-1) =  &
                                      rhs_mat(m_stdy+1:neqs,m_stdy+1:neqs)
      endif
!
!  Subtract the relevant quantity from the entries of a_mat corresponding
!  to diagonal elements of the density matrix (the arrays mindex is a
!  global variable, which will have been set by ldbl_set_rhsmat)
!
      do i=1,ndm
         if(i.ne.j_state)then
            m = mindex(i,i)
            if(m.lt.m_stdy)then
               a_mat(:,m) = a_mat(:,m) + b(:)
            elseif(m.gt.m_stdy)then
               a_mat(:,m-1) = a_mat(:,m-1) + b(:)
            else
               print*,'Logic error in ldbl_steadystate_prep ',i,m
               call ldbl_error
            endif
         endif
      enddo
!
      end subroutine ldbl_steadystate_prep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine ldbl_tdeig_rate(rhs_mat,rho,t,tstep,fl_cont,isteady)
!
!  The same as ldbl_tdeig, but here for calculations within the rate
!  equations approximation. Only the populations are propagated,
!  and rhs_mat is meant to be the rhs_mat matrix for the populations
!  only. The only difference with ldbl_tdeig is that the arrays are
!  of size ndm rather than ndm*ndm.
!
!  Propagate the vector rho over one time step of duration tstep by representing
!  this arrays as a linear combination of the right-hand eigenvectors of the
!  matrix rhs_mat, or obtain the steady state deriving from the vector rho
!  passed to ldbl_tdeig_rate.
!
!  The calculation of the steady state is identical to the calculation of
!  the time-dependent density matrix, apart that for the steady state only
!  the eigenvectors belonging to a zero eigenvalue are used. For the purpose
!  of the calculation, a zero eigenvalue is defined as one whose real and
!  imaginary parts are both smaller than the parameter eps defined below,
!  in absolute value.
!
!  rho: At entry, the array rho at time t (used only if fl_cont.eqv..false.).
!       At return, the arrays rho at time = t + tstep, or the steady state
!       density matrix if isteady.eq.1.
!  t: The initial time. At return, the value of t is reset to the final time,
!     t + tstep, unless isteady.eq.1.
!  fl_cont: A .false. value triggers the calculation of the eigenvalues and
!           eigenvectors of rhs_mat (which is overwritten) before the array
!           rho is calculated. A .true. value means that rho is calculated
!           using the coefficient of the expansion calculated at the last
!           call to ldbl_tdeig_rate with fl_cont .eqv. .false.
!  isteady: If isteady.eq.1, the subroutine returns the steady state density
!           matrix deriving from the input density matrix, and t is not
!           updated. This option is incompatible with fl_cont.eqv..true.
!           If isteady.eq.0, the calculation is done as normal. Any other
!           value of isteady is illegal. 

      implicit none
!
      real(kd) :: eps
!
      real(kd), dimension(ndm,ndm), intent(inout) :: rhs_mat
      real(kd), dimension(ndm), intent(inout) :: rho
      real(kd), intent(inout) :: t
      real(kd), intent(in) :: tstep
      integer, intent(in) :: isteady
      logical, intent(in) :: fl_cont
!
      character*1 :: jobvl,jobvr
      integer :: info,lda,ldvl,ldvr,lwork,n
      real(kd), dimension(ndm) :: wr,wi
      real(kd), dimension(ndm,ndm), save :: vl,vr
      real(kd), dimension(4*ndm) :: work
!
      complex(kd), dimension(ndm,ndm), save :: vl_cplx,vr_cplx
      complex(kd), dimension(ndm), save :: cm_0,w
      complex(kd), dimension(ndm) :: rho_cplx
      real(kd), save :: t_0
      real(kd) :: deltat
      integer :: i,j,m
      external dgeev,sgeev
!
!  Define the parameter eps, according to the precision
!
      if(kd.eq.kind(1.0))then
         eps = 1.e-3_kd
      else
         eps = 1.e-8_kd
      endif
!
!  Check that the input parameters are legal:
      if(isteady.eq.1)then
         if(fl_cont)then
            print*,'Wrong combination of isteady and fl_cont in ldbl_tdeig_rate.'
            call ldbl_error
         endif
      elseif(isteady.ne.0)then
         print*,'Wrong value of isteady in ldbl_tdeig_rate: ',isteady
         call ldbl_error
      endif
!
      if(.not.fl_cont)then
!
         jobvl = 'V'
         jobvr = 'V'
         n = ndm
         lda = ndm
         ldvl = ndm
         ldvr = ndm
         lwork = 4*ndm
         if(kd.eq.kind(1.0))then
            call sgeev(jobvl,jobvr,n,rhs_mat,lda,wr,wi,vl,ldvl,vr,ldvr,   &
                       work,lwork,info)
         elseif(kd.eq.kind(1.d0))then
            call dgeev(jobvl,jobvr,n,rhs_mat,lda,wr,wi,vl,ldvl,vr,ldvr,   &
                       work,lwork,info)
         else
            print*,'Logic error in ldbl_tdeig_rate...',kd,kind(1.0),kind(1.d0)
            call ldbl_error
         endif
         if(info.ne.0)then
            print*,'s/dgeev (from Lapack) returned with info = ',info
            print*,'in ldbl_tdeig_rate.'
            call ldbl_error
         endif

!  Calculate the coefficients of the expansion, and save the initial value of t.
         rho_cplx = rho
         m = 1
         do 
            if(m.gt.ndm)exit
            if(wi(m).eq.0.0_kd)then
               w(m) = cmplx(wr(m),0.0_kd,kd)
               vl_cplx(:,m) = cmplx(vl(:,m),0.0_kd,kd)
               vr_cplx(:,m) = cmplx(vr(:,m),0.0_kd,kd)
               if(abs(dot_product(vl_cplx(:,m),vr_cplx(:,m))).eq.0.0_kd)then
                  print*,'ldbl_tdeig_rate: zero eigenvector...'
                  call ldbl_error
               endif
               cm_0(m) = dot_product(vl_cplx(:,m),rho_cplx)/  &
                         dot_product(vl_cplx(:,m),vr_cplx(:,m))
               m = m + 1
            else
               if(m.ge.ndm)then
                  print*,'Logic Error in ldbl_tdeig_rate', m
                  call ldbl_error
               endif
               w(m) = cmplx(wr(m),wi(m),kd)
               vl_cplx(:,m) = cmplx(vl(:,m),vl(:,m+1),kd)
               vr_cplx(:,m) = cmplx(vr(:,m),vr(:,m+1),kd)
               if(abs(dot_product(vl_cplx(:,m),vr_cplx(:,m))).eq.0.0_kd)then
                  print*,'ldbl_tdeig_rate: zero eigenvector...'
                  call ldbl_error
               endif
               cm_0(m) = dot_product(vl_cplx(:,m),rho_cplx)/  &
                         dot_product(vl_cplx(:,m),vr_cplx(:,m))
               m = m + 1
               w(m) = conjg(w(m-1))
               vl_cplx(:,m) = conjg(vl_cplx(:,m-1))
               vr_cplx(:,m) = conjg(vr_cplx(:,m-1))
               cm_0(m) = conjg(cm_0(m-1))
               m = m + 1
            endif
         enddo
         t_0 = t
!
      endif
!
!  Calculate rho at time t + tstep by making a linear combination of the
!  right-hand eigenvectors of rhs_mat. The coefficient of the m-th eigenvector
!  at time t + tstep is the coefficient at time t_0 (saved in cm_0(m))
!  multiplied by exp[w(m)*(t+tstep-t_0)], where w(m) is the
!  corresponding eigenvalue (this eigenvalue can be complex).
!
!  Note that this coefficient, the exponential and the eigenvector are all
!  real if the eigenvalue is real, while if it is complex then the imaginary
!  parts cancel between terms corresponding to eigenvalues that are complex
!  conjugate of each other, with the consequence that the sum is always real
!  even though the eigenvalues and eigenvectors may be complex.
!
!  First, check that the initial vector is suitably represented by
!  a linear combination of right eigenvectors
      rho = 0.0_kd
      do m = 1,ndm
         rho = rho + real(cm_0(m)*vr_cplx(:,m),kd)
      enddo
!  rho_cplx still has the initial vector...
      do m = 1,ndm
         if(abs(rho(m)-real(rho_cplx(m),kd)).gt.eps)then
            print*,'ldbl_tdeig_rate: Warning of a possible issue with the basis'
            print*,'of right eigenvectors. ',m,abs(rho(m)-real(rho_cplx(m),kd))
         endif
      enddo
!
!  Normal calculation if isteady.eq.0, otherwise steady-state calculation
      if(isteady.eq.0)then
         deltat = t + tstep - t_0
         rho = 0.0_kd
         do m = 1,ndm
            rho = rho + real(exp(w(m)*deltat)*cm_0(m)*vr_cplx(:,m),kd)
         enddo
         t = t + tstep
      else
         rho = 0.0_kd
         do m = 1,ndm
            if(abs(real(w(m),kd)).lt.eps .and. abs(aimag(w(m))).lt.eps)then
               rho = rho + real(cm_0(m)*vr_cplx(:,m),kd)
            elseif(abs(real(w(m),kd)).lt.eps .and. abs(cm_0(m)).gt.eps)then
               print*,'ldbl_tdeig_rate finds oscillating populations.'
               print*,'The density matrix does not converge to'
               print*,'a constant matrix in the long time limit.'
               call ldbl_error
            endif
         enddo
      endif
!
      end subroutine ldbl_tdeig_rate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end module ldbl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine fcn_dummy
!
!  A dummy subprogram which does not do anything. It is specified
!  in the calling list of the subroutine dop853.
!
      end subroutine fcn_dummy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine solout_dummy
!
!  A dummy subprogram which does not do anything. It is specified
!  in the calling list of the subroutine dop853.
!
      end subroutine solout_dummy

