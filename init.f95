module init

use grid
use mpi_interface

contains

subroutine init_grids

implicit none

integer :: i,j

do i=1,N_vert
	zt(i)=hz*(dble(i)-0.5)
	zm(i)=hz*dble(i)
enddo

if(outcoords(2) > 0)then
	do i=1,x_dim
		xm(i)=(N_x-(outdims(2)-1)*(N_x/outdims(2)))*hx+x_dim*(outcoords(2)-1)*hx+hx*(dble(i)-1.0)
		xt(i)=(N_x-(outdims(2)-1)*(N_x/outdims(2)))*hx+x_dim*(outcoords(2)-1)*hx+hx*(dble(i)-0.5)
	enddo
else
	do i=1,x_dim
		xm(i)=hx*(dble(i)-1.0)
		xt(i)=hx*(dble(i)-0.5)
	enddo
endif

if(outcoords(1)>0)then
	do i=1,y_dim
		ym(i)=(N_y-(outdims(1)-1)*(N_y/outdims(1)))*hy+y_dim*(outcoords(1)-1)*hy+hy*(dble(i)-1.0)
		yt(i)=(N_y-(outdims(1)-1)*(N_y/outdims(1)))*hy+y_dim*(outcoords(1)-1)*hy+hy*(dble(i)-0.5)
	enddo
else
	do i=1,y_dim
		ym(i)=hy*(dble(i)-1.0)
		yt(i)=hy*(dble(i)-0.5)
	enddo
endif

end subroutine init_grids

subroutine init_conds

implicit none

integer :: i,j,k

do i=1,N_vert
	do j=1,x_dim
		do k=1,y_dim
			var(2+i,2+j,2+k,1)=0.0
			var(2+i,2+j,2+k,2)=0.0
			var(2+i,2+j,2+k,3)=0.0!3.0*exp((-1.0/(2.0*450.0**2))*((xt(j)- 3200.0)**2 + (zt(i)- 1800.0)**2))
			var(2+i,2+j,2+k,4)=3.0*exp((-1.0/(2.0*550.0**2))*((xt(j)- (N_x*hx/2.0))**2 + (yt(k)- (N_y*hy/2.0))**2 &
								+ (zt(i)- (N_vert*hz/2.0))**2))
!			if(zt(i)<500.0)then
!				var(2+i,2+j,2+k,4)=0.2*(2.0*rand()-1.0)
!			else
!				var(2+i,2+j,2+k,4)=0.0
!			endif
			var(2+i,2+j,2+k,4)=0.0
			var(2+i,2+j,2+k,5)=0.0!0.9*qvs(k)%p_var(z_dim_vec_ext(k)-z_dim_vec(k)+i)
			var(2+i,2+j,2+k,6)=0.0
		enddo
	enddo
enddo

!call random_number(var(3:7+2,3:x_dim+2,3:y_dim+2,4))

!var(3:10+2,3:x_dim+2,3:y_dim+2,4)=0.2*(2.0*var(3:10+2,3:x_dim+2,3:y_dim+2,4)-1.0)

end subroutine init_conds

subroutine calc_write_dim_freq_vec

implicit none

integer :: dim_freq,i

dim_freq=dint(Tend/write_freq)

allocate(write_freq_vec(dim_freq))

do i=1,dim_freq
	write_freq_vec(i)=(i-1.0)*write_freq
enddo

end subroutine calc_write_dim_freq_vec


end module init
