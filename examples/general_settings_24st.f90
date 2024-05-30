      module general_settings
!
      implicit none
!
      integer, parameter, public :: kd = kind(1.d0)  ! use double precision
!     integer, parameter, public :: kd = kind(1.0)   ! use single precision
!
      integer, parameter, public :: nst = 24         ! number of states
!
!  Lowest value of the index labelling the states outside the modules:
!
      integer, parameter, public :: nmn = 1      
!
      end module general_settings
