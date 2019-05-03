subroutine Func(A,u,r,lambda,beta,negative,myid,nprocs)

  use header
  implicit none
  include "mpif.h"
 
  type(Matrix), intent(inout)  :: A
  type(Vector), intent(inout)  :: u
  type(Vector), intent(inout)  :: r
  real(kind=8), intent(in)  :: lambda,beta
  integer, intent(in) :: myid,nprocs
  logical, intent(out) :: negative

  real(kind=8) :: div
  integer :: i,ierror

  negative = .false.
  !--------------------------------------------------------------------------------
  ! Ser r = A_delta*U
  !--------------------------------------------------------------------------------
  call Mat_Mult(A,u,r,myid,nprocs)

  !--------------------------------------------------------------------------------
  ! Ser -r = -A_delta*U + lambda*e^(u/1+beta*u)
  !--------------------------------------------------------------------------------
  do i=r%ibeg,r%iend
    ! Check for a negative temperature
    negative = u%xx(i) .lt. 0.0_8
    
    div = 1.0_8 + beta*u%xx(i) 
    r%xx(i) = lambda*exp(u%xx(i) / div) - r%xx(i)
  end do
  
  !--------------------------------------------------------------------------------
  ! Check for errors on all processes
  !--------------------------------------------------------------------------------
  call MPI_Allreduce(MPI_IN_PLACE,negative,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,ierror)

end subroutine Func
