module step

use grid
use mpi_interface
use netcdf_write
use therm
use rk_time_step

contains

subroutine time_stepper

implicit none

real*8 :: time=0.0
integer :: m=1
integer :: i,j


do while(time<Tend)
	call mpi_barrier(comm2d,ierr)
	if((write_freq_vec(m)-write_freq_dt <= time) .and. (time <= write_freq_vec(m)+write_freq_dt))then
		call calc_ql(ql,var(3:N_vert+2,3:x_dim+2,3:y_dim+2,5),qvs,N_Vert,x_dim,y_dim)
		call get_full_array(var(3:N_vert+2,3:x_dim+2,3:y_dim+2,:),ql)
		call write_netcdf_vars(write_array(:,:,:,1),write_array(:,:,:,2),write_array(:,:,:,3),&
									write_array(:,:,:,4),write_array(:,:,:,5),write_array(:,:,:,6),write_array(:,:,:,7),m)
		m=m+1
		write(*,*) m
	endif
	write(*,*) "Time =", time,myid
	call var_exchange(var)

	var(1,:,:,1)=var(4,:,:,1)
	var(2,:,:,1)=var(3,:,:,1)
	var(1,:,:,2)=var(4,:,:,2)
	var(2,:,:,2)=var(3,:,:,2)
	var(1,:,:,4)=var(4,:,:,4)
	var(2,:,:,4)=var(3,:,:,4)
	var(1,:,:,5)=var(4,:,:,5)
	var(2,:,:,5)=var(3,:,:,5)
	var(1,:,:,3)=0.0!-var(4,:,:,3)
	var(2,:,:,3)=0.0

	var(N_vert+3,:,:,1)=var(N_vert+2,:,:,1)
	var(N_vert+4,:,:,1)=var(N_vert+1,:,:,1)
	var(N_vert+3,:,:,2)=var(N_vert+2,:,:,2)
	var(N_vert+4,:,:,2)=var(N_vert+1,:,:,2)
	var(N_vert+3,:,:,4)=var(N_vert+2,:,:,4)
	var(N_vert+4,:,:,4)=var(N_vert+1,:,:,4)
	var(N_vert+3,:,:,5)=var(N_vert+2,:,:,5)
	var(N_vert+4,:,:,5)=var(N_vert+1,:,:,5)
	var(N_vert+3,:,:,3)=0.0
	var(N_vert+4,:,:,3)=0.0

	call rk_stepper(var(:,:,:,1),var(:,:,:,2),var(1:N_vert+3,:,:,3),var(:,:,:,4),&
					var(:,:,:,5),var(2:N_vert+3,2:x_dim+3,2:y_dim+3,6),N_vert,x_dim,y_dim,dt,hx,hy,hz)
	time=time+dt
enddo


end subroutine time_stepper

end module step
