! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! Subroutine that creates a file for visualise.m:
!               
!          Take a vector that corresponds to a 2D finite difference 
!          approximation on a uniform grid (e.g. the solution to the 
!          nonlinear thermal conduction problem) and write it to a
!          file in a format that can be used by the Matlab function 
!          postprocess.m to visualise it in a 3D surface plot using
!          the command
!                        visualise
!          in Matlab          
!               
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

subroutine save_solution(u,m)

  use header

  implicit none

  integer, intent(in)         :: m
  type(Vector), intent(inout) :: u

  real(kind=8), dimension(:,:), allocatable :: uu
  integer :: ierr,j,mdelta


  allocate(uu(m+1,m+1))

  do j=1,m
     uu(j,:) = u%xx((j-1)*(m+1)+1:j*(m+1))
  end do

  mdelta = m/2 + 1

  uu(m+1,1:mdelta)     = u%xx(m*(m+1)+1:m*(m+1)+mdelta)
  uu(m+1,mdelta+1:m+1) = 0
         
  open(3,file='solution.txt')

  write(3,*) m
  write(3,*) 
  write(3,*) uu(:,:)

  close(3)

  deallocate(uu)
  
end subroutine save_solution
