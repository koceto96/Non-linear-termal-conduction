!
!  Function to calculate the scalar product of two vectors, 
!  i.e to calculate d = a*b (parallel version).
!
!------------------------------------------------------------------------
function vec_dot(a,b) result (d)
  use header 
  implicit none
  include "mpif.h"

  type(Vector), intent(inout) :: a
  type(Vector), intent(inout) :: b
  real(kind=8) :: d
  real(kind=8) :: ddot
  integer :: ierror
  
  ! Compute local dot_product
  d = ddot(a%iend-a%ibeg+1, a%xx(a%ibeg), 1, b%xx(b%ibeg), 1)
  
  call MPI_Allreduce(MPI_IN_PLACE,d,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)

end function vec_dot
