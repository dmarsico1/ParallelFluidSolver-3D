module mpi_interface

use grid
use mpi

integer :: ierr,total_sum,comm,myid,num_proc,comm2d,ndim,plane
real*8 :: rold,rnew,rnew_tot,rold_tot,alpha_tot,gam,gam_tot,rhs_sum,rhs_sum_tot
integer :: status(mpi_status_size)
integer, dimension(2) :: dims
logical, dimension(2) :: isperiodic,outperiods
integer, dimension(2) :: outdims,outcoords,coords
integer, dimension(3) :: starts,sizes,subsizes
logical :: reorder
integer :: left,right,front,back,local_disp,recvsubarray,resizedrecvsubarray
integer, dimension(:), allocatable :: counts,disps
integer, dimension(:), allocatable :: reqs
real*8, dimension(:,:,:), allocatable :: p
real*8 :: c
integer(kind=mpi_address_kind) :: begin, extent

contains 

subroutine mpi_begin


dims(1)=0
dims(2)=0

isperiodic(1)=.true.
isperiodic(2)=.true.
reorder=.true.
ndim=2

call mpi_init(ierr)

call mpi_comm_size(mpi_comm_world,num_proc,ierr)

call mpi_dims_create(num_proc,ndim,dims,ierr)

call mpi_cart_create(mpi_comm_world,ndim,dims,isperiodic,reorder,comm2d,ierr)

call mpi_cart_get(comm2d,ndim,outdims,outperiods,outcoords,ierr)
!write(*,*) outdims

if(outcoords(2)==0)then
	x_dim=N_x-(outdims(2)-1)*(N_x/outdims(2))
else
	x_dim=N_x/outdims(2)
endif

if(outcoords(1)==0)then
	y_dim=N_y-(outdims(1)-1)*(N_y/outdims(1))
else
	y_dim=N_y/outdims(1)
endif

call mpi_comm_rank(comm2d,myid,ierr)


allocate(var(2*Ngz+N_vert,2*NGx+x_dim,2*NGy+y_dim,6),qvs(N_vert),ql(N_vert,x_dim,y_dim),qt_vec(N_vert),dqt_vec(N_vert),&
		dtheta_vec(N_vert),dtheta_e_vec(N_Vert),xm(x_dim),xt(x_dim),ym(y_dim),yt(y_dim),zm(N_vert),zt(N_vert),&
		spge(N_vert,x_dim,y_dim),sf(N_vert,x_dim,y_dim))

allocate(p(2+N_vert,2+x_dim,2+y_dim))

allocate(reqs(4*num_proc))

end subroutine mpi_begin


subroutine neighbors

implicit none

call mpi_cart_shift(comm2d,0,1,back,front,ierr)

call mpi_cart_shift(comm2d,1,1,left,right,ierr)

end subroutine neighbors


subroutine var_exchange(var)

implicit none
real*8, dimension(N_vert+4,x_dim+4,y_dim+4,6) :: var

integer :: i

do i=1,6
	
	call mpi_type_vector(x_dim,N_vert,2*NGz+N_vert,mpi_double_precision,plane,ierr)
	call mpi_type_commit(plane,ierr)
	
	call mpi_sendrecv(var(3,3,y_dim+2,i),1,plane,front,0,&
						var(3,3,2,i),1,plane,back,0,comm2d,status,ierr)
	
	call mpi_type_vector(x_dim,N_vert,2*NGz+N_vert,mpi_double_precision,plane,ierr)
	call mpi_type_commit(plane,ierr)

	call mpi_sendrecv(var(3,3,3,i),1,plane,back,1,&
						var(3,3,NGy+y_dim+1,i),1,plane,front,1,comm2d,status,ierr)
	
	call mpi_type_vector(x_dim,N_vert,2*NGz+N_vert,mpi_double_precision,plane,ierr)
	call mpi_type_commit(plane,ierr)

	call mpi_sendrecv(var(3,3,y_dim+1,i),1,plane,front,2,&
						var(3,3,1,i),1,plane,back,2,comm2d,status,ierr)
	
	call mpi_type_vector(x_dim,N_vert,2*NGz+N_vert,mpi_double_precision,plane,ierr)
	call mpi_type_commit(plane,ierr)

	call mpi_sendrecv(var(3,3,4,i),1,plane,back,3,&
						var(3,3,2*NGy+y_dim,i),1,plane,front,3,comm2d,status,ierr)

enddo

do i=1,6
	
	call mpi_type_vector(y_dim,N_vert,(x_dim+2*Ngx)*(2*NGz+N_vert),mpi_double_precision,plane,ierr)
	call mpi_type_commit(plane,ierr)

	call mpi_sendrecv(var(3,x_dim+2,3,i),1,plane,right,0,&
						var(3,2,3,i),1,plane,left,0,comm2d,status,ierr)

	call mpi_type_vector(y_dim,N_vert,(x_dim+2*Ngx)*(2*NGz+N_vert),mpi_double_precision,plane,ierr)
	call mpi_type_commit(plane,ierr)

	call mpi_sendrecv(var(3,3,3,i),1,plane,left,1,&
						var(3,NGx+x_dim+1,3,i),1,plane,right,1,comm2d,status,ierr)

	call mpi_type_vector(y_dim,N_vert,(x_dim+2*Ngx)*(2*NGz+N_vert),mpi_double_precision,plane,ierr)
	call mpi_type_commit(plane,ierr)

	call mpi_sendrecv(var(3,x_dim+1,3,i),1,plane,right,2,&
						var(3,1,3,i),1,plane,left,2,comm2d,status,ierr)

	call mpi_type_vector(y_dim,N_vert,(x_dim+2*Ngx)*(2*NGz+N_vert),mpi_double_precision,plane,ierr)
	call mpi_type_commit(plane,ierr)

	call mpi_sendrecv(var(3,4,3,i),1,plane,left,3,&
						var(3,2*NGx+x_dim,3,i),1,plane,right,3,comm2d,status,ierr)

enddo

end subroutine var_exchange


subroutine conj_exchange(p,k)

implicit none

real*8, dimension(N_vert+2,x_dim+2,y_dim+2) :: p
integer :: k
integer :: local_reqs(4)

call mpi_type_vector(x_dim,N_vert,2+N_vert,mpi_double_precision,plane,ierr)
call mpi_type_commit(plane,ierr)

call mpi_sendrecv(p(2,2,1+y_dim),1,plane,front,0,&
					p(2,2,1),1,plane,back,0,comm2d,status,ierr)

call mpi_type_vector(x_dim,N_vert,2+N_vert,mpi_double_precision,plane,ierr)
call mpi_type_commit(plane,ierr)

call mpi_sendrecv(p(2,2,2),1,plane,back,1,&
					p(2,2,y_dim+2),1,plane,front,1,comm2d,status,ierr)

call mpi_type_vector(y_dim,N_vert,(2+N_vert)*(2+x_dim),mpi_double_precision,plane,ierr)
call mpi_type_commit(plane,ierr)

call mpi_sendrecv(p(2,2,2),1,plane,left,2,&
					p(2,x_dim+2,2),1,plane,right,2,comm2d,status,ierr)

call mpi_type_vector(y_dim,N_vert,(2+N_vert)*(2+x_dim),mpi_double_precision,plane,ierr)
call mpi_type_commit(plane,ierr)

call mpi_sendrecv(p(2,x_dim+1,2),1,plane,right,3,&
					p(2,1,2),1,plane,left,3,comm2d,status,ierr)





end subroutine conj_exchange


subroutine pressure_exchange(p,k)

implicit none

real*8, dimension(N_vert+2,x_dim+2,y_dim+2) :: p
integer :: k

call mpi_type_vector(x_dim,N_vert,2+N_vert,mpi_double_precision,plane,ierr)
call mpi_type_commit(plane,ierr)

call mpi_sendrecv(p(2,2,1+y_dim),1,plane,front,0,&
					p(2,2,1),1,plane,back,0,comm2d,status,ierr)

call mpi_type_vector(x_dim,N_vert,2+N_vert,mpi_double_precision,plane,ierr)
call mpi_type_commit(plane,ierr)

call mpi_sendrecv(p(2,2,2),1,plane,back,1,&
					p(2,2,y_dim+2),1,plane,front,1,comm2d,status,ierr)

call mpi_type_vector(y_dim,N_vert,(2+N_vert)*(2+x_dim),mpi_double_precision,plane,ierr)
call mpi_type_commit(plane,ierr)

call mpi_sendrecv(p(2,2,2),1,plane,left,2,&
					p(2,x_dim+2,2),1,plane,right,2,comm2d,status,ierr)

call mpi_type_vector(y_dim,N_vert,(2+N_vert)*(2+x_dim),mpi_double_precision,plane,ierr)
call mpi_type_commit(plane,ierr)

call mpi_sendrecv(p(2,x_dim+1,2),1,plane,right,3,&
					p(2,1,2),1,plane,left,3,comm2d,status,ierr)

end subroutine pressure_exchange

subroutine get_ranks

implicit none

integer :: i

if(outcoords(2)>0)then
	local_disp=(N_x-(outdims(2)-1)*(N_x/outdims(2)))*N_vert+(outcoords(2)-1)*N_vert*(N_x/outdims(2))
endif

if(outcoords(1) > 0)then
	local_disp=local_disp+(N_y-(outdims(1)-1)*(N_y/outdims(1)))*N_x*N_vert+(outcoords(1)-1)*(N_y/outdims(1))*N_vert*N_x
endif

if(myid==0)then
	allocate(counts(num_proc),disps(num_proc))
	counts=1
endif

call mpi_gather(local_disp,1,mpi_int,disps,1,mpi_int,0,comm2d,ierr)


end subroutine get_ranks


subroutine get_full_array(var,ql)

implicit none
real*8, dimension(N_vert,x_dim,y_dim,6) :: var
real*8, dimension(N_Vert,x_dim,y_dim) :: ql
integer :: i

do i=1,6
	sizes(1)=N_vert; sizes(2)=N_x; sizes(3)=N_y
	subsizes(1)=N_vert; subsizes(2)=x_dim; subsizes(3)=y_dim
	starts(1)=0; starts(2)=0; starts(3)=0
	
	call mpi_type_create_subarray(3, sizes, subsizes, starts, mpi_order_fortran, &
									mpi_double_precision, recvsubarray, ierr)
	
	call mpi_type_commit(recvsubarray,ierr)
	
	extent = sizeof(c)
	begin = 0
	call mpi_type_create_resized(recvsubarray, begin, extent, resizedrecvsubarray, ierr)
	call mpi_type_commit(resizedrecvsubarray,ierr)
	!counts=1
	call mpi_gatherv(var(1,1,1,i),x_dim*y_dim*N_vert,mpi_double_precision,write_array(:,:,:,i),counts,disps,&
						resizedrecvsubarray,0,comm2d,ierr)

enddo

	sizes(1)=N_vert; sizes(2)=N_x; sizes(3)=N_y
	subsizes(1)=N_vert; subsizes(2)=x_dim; subsizes(3)=y_dim
	starts(1)=0; starts(2)=0; starts(3)=0
	
	call mpi_type_create_subarray(3, sizes, subsizes, starts, mpi_order_fortran, &
									mpi_double_precision, recvsubarray, ierr)
	
	call mpi_type_commit(recvsubarray,ierr)
	
	extent = sizeof(c)
	begin = 0
	call mpi_type_create_resized(recvsubarray, begin, extent, resizedrecvsubarray, ierr)
	call mpi_type_commit(resizedrecvsubarray,ierr)

call mpi_gatherv(ql(1,1,1),x_dim*y_dim*N_vert,mpi_double_precision,write_array(:,:,:,7),counts,disps,&
						resizedrecvsubarray,0,comm2d,ierr)

end subroutine get_full_array

subroutine mpi_end

implicit none

call mpi_finalize(ierr)

end subroutine mpi_end


end module mpi_interface
