! -----------------------------------------------------------------------
subroutine Maximum(u,umax)
! -----------------------------------------------------------------------
!
!  Description: Subroutine to calculate  || u ||_inf  =  max |u| 
!
! -----------------------------------------------------------------------

  use header
  include "mpif.h"

  type(Vector), intent(inout) :: u
  real(kind=8), intent(out)   :: umax

  real(kind=8) :: mymax

  mymax = maxval(u%xx(u%ibeg:u%iend))
  
  call MPI_Reduce(mymax,umax,1,MPI_DOUBLE_PRECISION,MPI_MAX,0, &
                  MPI_COMM_WORLD,ierr)

end subroutine Maximum
