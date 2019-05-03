!
!  Description: Newton's method (parallel version)
!
! -----------------------------------------------------------------------
subroutine Newton(A,u,lambda,beta,tau,kmax,k,ierr,myid,nprocs) 
  
  use header
  implicit none
  include "mpif.h"

  type(Matrix), intent(in)  :: A
  type(Vector), intent(inout)  :: u
  real(kind=8), intent(in)  :: lambda,beta,tau
  integer, intent(in) :: myid,nprocs
  integer, intent(in) :: kmax
  integer, intent(out) :: k
  logical, intent(out) :: ierr
 
  type(Matrix) :: J
  type(Vector) :: r
  type(Vector) :: u_last

  real(kind=8) :: Vec_Dot
  real(kind=8) :: tol,last_norm,norm
  real(kind=8) :: t1,t2,t_total

  integer :: cg_kmax,cg_k,cg_total

  cg_kmax = 9999
  
  t_total = 0.0_8
  cg_total = 0
 
  allocate(r%xx(u%n))
  r%n = u%n
  r%ibeg = u%ibeg
  r%iend = u%iend
  
  allocate(u_last%xx(u%n))
  u_last%n = u%n
  u_last%ibeg = u%ibeg
  u_last%iend = u%iend

  ! Allocate the jacobian as a copy of A
  allocate(J%aa(A%nnz),J%jj(A%nnz),J%ii(A%n+1))
  J%n = A%n
  J%bw = A%bw
  J%nnz = A%nnz
  J%ibeg = A%ibeg
  J%iend = A%iend
  J%aa = A%aa
  J%jj = A%jj
  J%ii(J%ibeg:J%iend+1) = A%ii(J%ibeg:J%iend+1)
  
  !--------------------------------------------------------------------------------
  ! Begin Newton's method
  !--------------------------------------------------------------------------------
   
  if (myid .eq. 0) then
     print '(A)', "|---k---|------CG iterations------|---------||r||---------|"  
  end if
  !--------------------------------------------------------------------------------
  ! Set r = F(U)
  !--------------------------------------------------------------------------------
  call Func(A,u,r,lambda,beta,ierr,myid,nprocs) 
  if (ierr) goto 10
  
  tol = 0.1_8
  
  do k=0,kmax-1 
    !--------------------------------------------------------------------------------
    ! Calculate the norm of r
    !--------------------------------------------------------------------------------
    norm = sqrt(Vec_Dot(r,r))
    if (norm .le. tau) then
        exit
    end if 
     
    !--------------------------------------------------------------------------------
    ! Calculate F'(U)
    !--------------------------------------------------------------------------------
    call Jacobian(A,u,lambda,beta,J) 
    
    !--------------------------------------------------------------------------------
    ! Set the maximum tolerance according to the formula
    !--------------------------------------------------------------------------------
    if (k .ge. 1) then
       tol = min(0.1_8,0.9_8*(norm**2)/(last_norm)**2)
    end if
    
    !--------------------------------------------------------------------------------
    ! Calculate F'(U)s=-r and time the number of CG iterations
    !--------------------------------------------------------------------------------
    
    t1 = MPI_WTIME()
    call cg(J,u_last,r,tol,cg_kmax,cg_k,myid,nprocs)
    t2 = MPI_WTIME()
    
    ! Compute the time in CG up to now and the number of iterationss
    cg_total = cg_total + cg_k
    t_total = t_total + (t2-t1) 

    !--------------------------------------------------------------------------------
    ! U = U + s
    !--------------------------------------------------------------------------------
    call daxpy(u%iend-u%ibeg+1,1.0_8,u_last%xx(u_last%ibeg),1,u%xx(u%ibeg),1)    
    !--------------------------------------------------------------------------------
    ! Set r = F(U)
    !--------------------------------------------------------------------------------
    call Func(A,u,r,lambda,beta,ierr,myid,nprocs)
    if (ierr) goto 10
  
    last_norm = norm
    
    if (myid .eq. 0) then
        print '(AI5XXAX4XI5X14XAX5XES15.7XXA)',"|",k,"|",cg_k,"|",norm,"|"
    end if
  end do
 
  !--------------------------------------------------------------------------------
  ! End Newton's method
  !--------------------------------------------------------------------------------
  if(myid == 0) then  
    print '(AI6)', "The total CG iterations:", cg_total
    print '(AF15.7)', "Time spent in CG iterations:", t_total
    print '(AF15.7)', "Time for 1 CG iteration:", t_total/real(cg_total,8)
  end if
  
  10 deallocate(J%aa,J%ii,J%jj,r%xx,u_last%xx)

end subroutine Newton 
