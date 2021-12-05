module therm

use grid
use mpi_interface

real*8, dimension(:), allocatable :: qvs_tot

contains

subroutine init_qvs

implicit none

integer :: i

allocate(qvs_tot(N_vert))

call init_qvs_grid_k(1)

call calc_qt_vec
call calc_dqt_vec

call init_qvs_grid_k(2)


end subroutine init_qvs

subroutine init_qvs_grid_k(m)

implicit none

integer :: m

if(m==1)then
	call calc_qvs(qvs_tot,zt,N_vert)
else
	call calc_qvs_pert(qvs,qvs_tot,zt,N_vert)
endif


end subroutine init_qvs_grid_k


subroutine calc_qvs(qvs,zt,M)

!calculate the total saturation mixing ratio

implicit none

real*8 :: p0=10**5
real*8 :: e_s0=3500.0
integer, intent(in) :: M
real*8, dimension(M) :: qvs,zt
real*8, dimension(M) :: theta_tilde,p_tilde,T_tilde
integer :: i

do i=1,M
	theta_tilde(i)=zt(i)*dtheta+th0
enddo

do i=1,M
	p_tilde(i)=p0*(1-(g/(cp*dtheta))*log(1+(dtheta/th0)*(zt(i))))**(cp/Rd)
enddo

do i=1,M
	T_tilde(i)=theta_tilde(i)*(p_tilde(i)/p0)**(Rd/cp)
enddo

do i=1,M
	qvs(i)=(Rd/Rv)*(e_s0/p_tilde(i))*exp(-(Lv/Rv)*((1.0/T_tilde(i))-(1.0/th0)))
enddo


end subroutine calc_qvs

subroutine calc_qvs_pert(qvs_pert,qvs_tot,zt,M)

!calculate the perturbation saturation mixing ratio by subtracting off the backround total water

implicit none

integer, intent(in) :: M
real*8, dimension(M) :: qvs_pert,qvs_tot,zt
integer :: i


do i=1,M
	qvs_pert(i)=qvs_tot(i)-qt_vec(i)
enddo


end subroutine calc_qvs_pert

function VT_func(z)

implicit none

real*8, intent(in) :: z
real*8 :: rho_0,rho
real*8 :: gam=-0.3654
real*8 :: f0=14.34
real*8 :: p0=10.0**5
real*8 :: p_tilde,T_tilde,theta_tilde
integer :: i
real*8 :: VT_func

rho_0=p0/(Rd*th0)

theta_tilde=z*dtheta+th0

p_tilde=p0*(1-(g/(cp*dtheta))*log(1+(dtheta/th0)*(z)))**(cp/Rd)

T_tilde=theta_tilde*(p_tilde/p0)**(Rd/cp)

rho=p_tilde/(Rd*T_tilde)

VT_func=f0*sqrt(rho_0)*rho**gam

end function VT_func


function qvs_func(z)

implicit none

real*8 :: z
real*8 :: p0=10**5
real*8 :: e_s0=3500.0
real*8 :: theta_tilde,p_tilde,T_tilde,qvs,qvs_pert
real*8 :: qvs_func

theta_tilde=z*dtheta+th0

p_tilde=p0*(1-(g/(cp*dtheta))*log(1+(dtheta/th0)*(z)))**(cp/Rd)

T_tilde=theta_tilde*(p_tilde/p0)**(Rd/cp)

qvs=(Rd/Rv)*(e_s0/p_tilde)*exp(-(Lv/Rv)*((1.0/T_tilde)-(1.0/th0)))

qvs_func=qvs-0.5*qvs


end function

!subroutine init_sponge

!implicit none

!call sponge_grid_k(N_vert,x_dim)

!end subroutine init_sponge

subroutine sponge_grid(M,N,L)

implicit none

integer :: i,j,M,N,L,k

do i=1,M
	if(zt(i)<N_vert*hz*spge_h)then
		spge(i,:,:)=0.0
	else
	do j=1,N
		do k=1,L
			spge(i,j,k)=(1.0/3600.0)*(sin((0.5*4.0*atan(1.0)/(N_vert*hz-spge_h*N_vert*hz))*(zt(i)&
										-N_vert*hz*spge_h))**2)
		enddo
	enddo
	endif
enddo


end subroutine sponge_grid

subroutine sponge_layer(u,t,M,N,L)

implicit none

real*8, dimension(M,N,L), intent(out) :: t
real*8, dimension(M,N,L), intent(in) :: u
integer :: i,j,M,N,k,L

do i=1,M
	do j=1,N
		do k=1,L
			t(i,j,k)=-spge(i,j,k)*u(i,j,k)
		enddo
	enddo
enddo

end subroutine sponge_layer

subroutine force(f,Bu,Bs,qvs,M,N,L)

implicit none

real*8, intent(out), dimension(M-1,N,L) :: f
real*8, dimension(M,N,L) :: Bu,Bs
real*8, dimension(M) :: qvs
integer, intent(in) :: M,N,L
integer :: i,j,k
real*8 :: row_sum,row_sum_tot



do i=1,M-1
	do j=1,N
		do k=1,L
			if(Bs(i,j,k)<qvs(i))then
				f(i,j,k)=0.5*g*((1.0/th0)*(Bu(i,j,k)+Bu(i+1,j,k))+(Rvd-(Lv/(cp*th0)))*(Bs(i,j,k)+Bs(i+1,j,k)))
			else
				f(i,j,k)=0.5*g*((1.0/th0)*(Bu(i,j,k)+Bu(i+1,j,k))+(Rvd-(Lv/(cp*th0))+1)*(qvs(i)+qvs(i+1))-(Bs(i,j,k)+Bs(i+1,j,k)))
			endif
		enddo
	enddo
enddo


!f=(g/th0)*Bu(1:M-1,:)

do i=1,M-1
	row_sum=sum(f(i,:,:))/(N_X*N_y)
	call mpi_allreduce(row_sum,row_sum_tot,1,mpi_real8,mpi_sum,comm2d,ierr)
	!if(row_sum >=0)then
		do j=1,N
			do k=1,L
				f(i,j,k)=f(i,j,k)-row_sum_tot
			enddo
	!	enddo
	!else
	!	do j=1,N
	!		f(i,j)=f(i,j)-row_sum
		enddo
	!endif
enddo


end subroutine force

subroutine calc_dthetae

!calculate derivative of background equivalent potential temperature

implicit none

dthetae=dtheta+(Lv/cp)*dqt

end subroutine calc_dthetae


subroutine calc_dtheta_e_vec

implicit none

integer :: i

do i=1,N_vert
	if(zt(i)<=BL_top)then
		dtheta_e_vec(i)=(Lv/cp)*dqt_vec(i)
	else
		dtheta_e_vec(i)=dtheta+(Lv/cp)*dqt_vec(i)
	endif
enddo


end subroutine calc_dtheta_e_vec


subroutine calc_qt_vec

implicit none

integer :: i

if(qt_tilde_type==1)then
	do i=1,N_vert
		!if(zt(i)<=BL_top)then
			!qt_vec(i)=qt0
		!else
			qt_vec(i)=0.0!0.87*qvs_tot(i)
		!endif
	enddo

elseif(qt_tilde_type==0)then
		do i=1,N_vert
			if(zt(i)<=BL_top)then
				qt_vec(i)=qt0
			else
				qt_vec(i)=qt0+dqt*(zt(i)-BL_top)
			endif
		enddo
endif


end subroutine calc_qt_vec


subroutine calc_dqt_vec

implicit none

integer :: i

if(qt_tilde_type==1)then
	do i=1,N_vert
		!if(zt(i)<=BL_top)then
		!	dqt_vec(i)=0.0
		!else
			if(i==1)then
				dqt_vec(i)=(qt_vec(i+1)-qt_vec(i))/hz
			elseif(i==N_vert)then
				dqt_vec(i)=(3.0*qt_vec(i)-4.0*qt_vec(i-1)+qt_vec(i-2))/(2.0*hz)
			else
				dqt_vec(i)=(qt_vec(i+1)-qt_vec(i-1))/(2.0*hz)
			!dqt_vec(i)%p_var(j)=(dqt/qt0)*qt0*exp((dqt/qt0)*(1.0/nls*hz)*(z_grid_ext(i)%p_var(j)-BL_top))
			endif
			!dqt_vec(i)%p_var(j)=(dqt/qt0)*(1.0/(nls*hz))*qt0*exp((dqt/qt0)*(1.0/nls*hz)*(z_grid_ext(i)%p_var(j)-BL_top))
		!endif
	enddo

elseif(qt_tilde_type==0)then
		do i=1,N_vert
			if(zt(i)<=BL_top)then
				dqt_vec(i)=0.0!qt0
			else
				dqt_vec(i)=dqt
			endif
		enddo
endif

end subroutine calc_dqt_vec

subroutine calc_ql(ql,qt,qvs,M,N,L)

implicit none

real*8, dimension(M,N,L) :: qt,ql
real*8, dimension(M) :: qvs
integer :: M,N,L
integer :: i,j,k

do i=1,M
	do j=1,N
		do k=1,L
			ql(i,j,k)=max(qt(i,j,k)-qvs(i),0.0)
		enddo
	enddo
enddo

end subroutine calc_ql

!subroutine init_surf_flux

!integer :: k

!call surf_flux(sf(k)%p_var(:,:,1),z_dim_vec_ext(k),x_dim_vec(k))


!end subroutine init_surf_flux

!subroutine surf_flux(sf,M,N)

!implicit none

!integer, intent(in) :: M,N
!real*8, intent(out), dimension(M,N) :: sf
!integer :: i,j

!do i=1,M
!	do j=1,N
!		if(zt(i)<BL_top)then
!			sf(i,j)=0.0*(1/10800)*(sin((0.5*4.0*atan(1.0)/(BL_top))*(zt(i)&
!										-BL_top))**2)!(1.0/3600.0)*exp(-z_grid_ext(grid_num)%p_var(i)/(N_vert*hz))
!		else
!			sf(i,j)=0.0
!		endif
!	enddo
!enddo

!end subroutine surf_flux



subroutine ls_force_vel(ls_fu,U,M,N,L)

implicit none

integer :: i,M,N,L
real*8, dimension(M,N,L) :: ls_fu,U

do i=1,M
	ls_fu(i,:,:)=-(sum(U(i,:,:))/(N*L))/relax_u
enddo

end subroutine ls_force_vel

subroutine ls_force_therm(ls_fBu,ls_fBs,Bu,Bs,M,N,L)

implicit none

integer :: i,M,N,L,j,k
real*8, dimension(M,N,L) :: ls_fBu,ls_fBs,Bu,Bs


do i=1,M
	do j=1,N
		do k=1,L
			ls_fbs(i,j,k)=max_moist*(1.0/therm_relax)*(1.0/(max_moist_z*N_vert*hz))*zt(i)*&
						exp(-(1.0/2.0)*((1.0/(max_moist_z*N_vert*hz))**2)*zt(i)**2)
		enddo
	enddo
enddo


do i=1,M
	do j=1,N
		do k=1,L
			ls_fbu(i,j,k)=max_t*(1.0/therm_relax)*(1.0/(max_t_z*N_vert*hz))*zt(i)*&
						exp(-(1.0/2.0)*((1.0/(max_t_z*N_vert*hz))**2)*zt(i)**2)&
						+(Lv/cp)*ls_fbs(i,j,k)
		enddo
	enddo
enddo


end subroutine ls_force_therm

end module therm

