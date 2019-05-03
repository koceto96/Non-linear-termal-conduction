subroutine sparsegather(u,bw,myid,nprocs)

  use header
  implicit none
  include "mpif.h"

  type(Vector), intent(inout) :: u 
  integer, intent(in) :: myid,nprocs,bw
  integer:: req1,req2,req3,req4,stat,ierror
  
  !--------------------------------------------------------------------------------
  ! Exchange the bandwith chunks from the left to the right neighbouring processes
  !--------------------------------------------------------------------------------
  if(myid .ne. 0) then
    call MPI_IRecv(u%xx(u%ibeg-bw:u%ibeg-1),bw,MPI_DOUBLE_PRECISION,myid-1,myid,MPI_COMM_WORLD,req1,ierror) 
  end if
  if(myid .lt. nprocs-1) then
    call MPI_ISend(u%xx(u%iend-bw+1:u%iend),bw,MPI_DOUBLE_PRECISION,myid+1,myid+1,MPI_COMM_WORLD,req2,ierror)
  end if
   
  !--------------------------------------------------------------------------------
  ! Wait for the operation to finish
  !--------------------------------------------------------------------------------
    
  if(myid .ne. 0) then
    call MPI_WAIT(req1,stat,ierror)
  end if
    
  !--------------------------------------------------------------------------------
  ! Exchange the bandwith chunks from the right to the left neighbouring processes
  !--------------------------------------------------------------------------------
    

  if(myid .lt. nprocs-1) then
     call MPI_IRecv(u%xx(u%iend+1:u%iend+bw),bw,MPI_DOUBLE_PRECISION,myid+1,myid,MPI_COMM_WORLD,req3,ierror)
  end if
  if(myid .ne. 0) then
    call MPI_ISend(u%xx(u%ibeg:u%ibeg+bw-1),bw,MPI_DOUBLE_PRECISION,myid-1,myid-1,MPI_COMM_WORLD,req4,ierror)
  end if

  !--------------------------------------------------------------------------------
  ! Wait for the operation to finish
  !--------------------------------------------------------------------------------
    
  if(myid .lt. nprocs-1) then
      call MPI_WAIT(req3,stat,ierror)
  end if
 
end subroutine sparsegather
