subroutine Jacobian(A,u,lambda,beta,J)

  use header
  implicit none 

  type(Matrix), intent(inout)  :: A
  type(Vector), intent(inout)  :: u
  type(Matrix), intent(inout)  :: J
  real(kind=8), intent(in)  :: lambda,beta

  real(kind=8)  :: div
  integer :: i
  !--------------------------------------------------------------------------------
  ! Compute J = A_delta - G'(U)
  !--------------------------------------------------------------------------------
  
  do i=u%ibeg,u%iend
    div = 1.0_8 + beta*u%xx(i)
     
    J%aa(J%ii(i)) = A%aa(J%ii(i)) - lambda*(1.0_8/(div**2))*exp((1.0_8/div)*u%xx(i))
  end do

end subroutine Jacobian
