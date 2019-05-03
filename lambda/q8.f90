!
!  Description: Numerical solution of the nonlinear thermal conduction 
!               equation in a partially insulated slab (parallel version)
!
! -----------------------------------------------------------------------

program main

  use header
  implicit none
  include "mpif.h"

  real (kind=8) :: beta, lambda, tau, lambda_step
  integer :: kmax,m
  parameter (m = 64, beta = 0.25_8, tau = 1.0d-5, kmax = 20, lambda_step = 1d-2)
  
  type(Matrix)  :: A
  type(Vector)  :: u
  real (kind=8) :: umax
  integer       ::  n, its
  integer       :: myid, numprocs, nrows, ibeg, iend
  logical       :: ierr
  integer       :: i

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Beginning of program - Initialisation of MPI context
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, numprocs, ierr)

  if ( myid == 0 ) then
    open(unit=3,file='lambda.dat')
  end if

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Broadcast m to the other processes
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  n = (m+1)*m + m/2 + 1

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Calculate the start and end indices of the rows to be held locally
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  ibeg = int(myid*real(n)/numprocs)+1
  iend = int((myid+1)*real(n)/numprocs)
  nrows = iend-ibeg+1

  allocate(A%aa(5*nrows))
  allocate(A%jj(5*nrows))
  allocate(A%ii(n+1))

  A%n    = n
  A%ibeg = ibeg
  A%iend = iend

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     Construct the linear part of the Jacobian 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  call Laplace(A,m,ibeg,iend) 

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     Allocate space for u and set the initial guess to 0
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  allocate(u%xx(n))

  u%n    = n
  u%ibeg = ibeg
  u%iend = iend
  u%xx(u%ibeg:u%iend) = 0.0d0
  
  lambda = 0.0_8
  do i=0,int(floor(1.5_8/lambda_step))
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     Apply Newton's method to solve the system
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ierr = .false.
    call Newton(A,u,lambda,beta,tau,kmax,its,ierr,myid,numprocs)
    if (ierr) then
      if (myid == 0) then
        print '(AI4)', "FATAL: Negative temperature at iteration:", its
      end if
      goto 10
    end if
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     Calculate the maximum temperature
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    call Maximum(u,umax)
  
    if (myid == 0) then

      if (its > kmax) then
         print '(AI4A)','The Newton Method is diverging. Maximum number',kmax,' of iterations attained.'
       else
         print '(AI4AF17.12)','After',its,' Newton steps the maximum temperature umax =',umax
         write (3,*) lambda, umax
       end if

    end if
    lambda = lambda + lambda_step
  end do
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     Write the solution to a file for postprocessing in Matlab
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (myid .eq. 0 ) then
    close(3) 
  end if 

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     Deallocate memory
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  10 deallocate(A%aa)
  deallocate(A%jj)
  deallocate(A%ii)
  deallocate(u%xx)

  call MPI_Finalize(ierr)

end program main
