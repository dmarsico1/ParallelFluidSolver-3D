module rk_time_step

use grid
use mpi_interface
use therm
use advect
use diffuse
use poisson
use matrix

contains

subroutine rk_stepper(U,V,W,Bu,Bs,P,M,N,L,dt,hx,hy,hz)


implicit none

integer :: M,N,L,grid_num
real*8, dimension(M+2*NGz,N+2*NGx,L+2*NGy), intent(inout) :: U,V,Bu,Bs
real*8, dimension(M+2,N+2,L+2), intent(inout) :: P
real*8, dimension(M+2*NGz-1,N+2*NGx,L+2*NGy), intent(inout) :: W
real*8, dimension(M,N,L) :: Px,Py
real*8, dimension(M-1,N,L) :: Pz
real*8 :: U1(M+2*NGz,N+2*NGx,L+2*NGy,3),V1(M+2*NGz,N+2*NGy,L+2*NGy,3),W1(M+2*NGz-1,N+2*NGx,L+2*NGy,3),Bu1(M+2*NGz,N+2*NGx,L+2*NGy,3)
real*8 :: Bs1(M+2*NGz,N+2*NGx,L+2*NGy,3)
real*8 :: dt,hx,hz,Ubar,hy
integer :: i,j,k

do i=1,3
	U1(:,:,:,i)=U
	V1(:,:,:,i)=V
	W1(:,:,:,i)=W
	Bu1(:,:,:,i)=Bu
	Bs1(:,:,:,i)=Bs
enddo

call vel_rk_step(U1(NGz+1:NGz+M,NGx+1:NGx+N,NGy+1:NGy+L,2),U1(:,:,:,1),U1(:,:,:,1),U1(:,:,:,1),&
				V1(NGz+1:NGz+M,NGx+1:NGx+N,NGy+1:NGy+L,2),V1(:,:,:,1),V1(:,:,:,1),V1(:,:,:,1),&
				W1(NGz+1:NGz+M-1,NGx+1:NGx+N,NGy+1:NGy+L,2),W1(:,:,:,1),W1(:,:,:,1),W1(:,:,:,1),&
				Bu1(:,:,:,1),Bs1(:,:,:,1),P,M,N,L,dt,hx,hy,hz)



!call PeriodicBoundary(U1(:,:,2),W1(:,:,2),Bu1(:,:,2),Bs1(:,:,2),P,M,N)

W1(NGz+M,NGx+1:NGx+N,NGy+1:NGy+L,2)=0.0

call PoissonSolve(U1(NGz+1:NGz+M,NGx+1:NGx+N+1,NGy+1:NGy+L,2),V1(NGz+1:NGz+M,NGx+1:NGx+N,NGy+1:NGy+L+1,2),&
					W1(NGz:NGz+M,NGx+1:NGx+N,NGy+1:NGy+NGy+L,2),P,N,M,L,dt,dble(10.0**(-3)),hx,hy,hz)

call pressure_exchange(P,1)

call GradP(Px,Py,Pz,P(2:M+1,1:N+1,1:L+1),hx,hy,hz,M,N,L)

U1(NGz+1:NGz+M,NGx+1:NGx+N,NGy+1:NGy+L,2)=U1(NGz+1:NGz+M,NGx+1:NGx+N,NGy+1:NGy+L,2)-dt*Px
V1(NGz+1:NGz+M,NGx+1:NGx+N,NGy+1:NGy+L,2)=V1(NGz+1:NGz+M,NGx+1:NGx+N,NGy+1:NGy+L,2)-dt*Py
W1(NGz+1:NGz+M-1,NGx+1:NGx+N,NGy+1:NGy+L,2)=W1(NGz+1:NGz+M-1,NGx+1:NGx+N,NGy+1:NGy+L,2)-dt*Pz



call B_rk_time_step(Bu1(NGz+1:NGz+M,NGx+1:NGx+N,NGy+1:NGy+L,2),Bu1(:,:,:,1),Bu1(:,:,:,1),Bu1(:,:,:,1),&
					Bs1(NGz+1:NGz+M,NGx+1:NGx+N,NGy+1:NGy+L,2),Bs1(:,:,:,1),Bs1(:,:,:,1),Bs1(:,:,:,1),&
					U1(:,:,:,2),V1(:,:,:,2),W1(:,:,:,2),M,N,L,dt,hx,hy,hz)


!call PeriodicBoundary(U1(:,:,2),W1(:,:,2),Bu1(:,:,2),Bs1(:,:,2),P,M,N)


U(NGz+1:NGz+M,NGx+1:NGx+N,NGy+1:NGy+L)=U1(NGz+1:NGz+M,NGx+1:NGx+N,NGy+1:NGy+L,2)
V(NGz+1:NGz+M,NGx+1:NGx+N,NGy+1:NGy+L)=V1(NGz+1:NGz+M,NGx+1:NGx+N,NGy+1:NGy+L,2)
W(NGz+1:NGz+M,NGx+1:NGx+N,NGy+1:NGy+L)=W1(NGz+1:NGz+M,NGx+1:NGx+N,NGy+1:NGy+L,2)

Bu(NGz+1:NGz+M,NGx+1:NGx+N,NGy+1:NGy+L)=Bu1(NGz+1:NGz+M,NGx+1:NGx+N,NGy+1:NGy+L,2)
Bs(NGz+1:NGz+M,NGx+1:NGx+N,NGy+1:NGy+L)=Bs1(NGz+1:NGz+M,NGx+1:NGx+N,NGy+1:NGy+L,2)

end subroutine rk_stepper


subroutine vel_rk_step(U,U1,U2,U3,V,V1,V2,V3,W,W1,W2,W3,Bu1,Bs1,P,M,N,L,dt,hx,hy,hz)

implicit none
real*8, dimension(M,N,L) :: U,V
real*8, dimension(M-1,N,L) :: W
real*8, dimension(M+2,N+2,L+2) :: P
real*8, dimension(M+4,N+4,L+4) :: U1,U2,U3,V1,V2,V3
real*8, dimension(M+3,N+4,L+4) :: W1,W2,W3
real*8, dimension(M+4,N+4,L+4) :: Bu1,Bs1
real*8, dimension(M,N,L) :: Ua,Ud,Va,Vd,Px,Py,ls_fu,Uhd,Vhd,U_spg,W_spg
real*8, dimension(M-1,N,L) :: f
real*8, dimension(M-1,N,L) :: Wa,Wd,Pz,Whd
integer :: M,N,L
integer :: i,j,k
real*8 :: hx,hz,hy,dt

call advectVel(Ua,Va,Wa,U2,V2,W2,M,N,L,hx,hy,hz)
call diffVel(Uhd,U1(NGz:M+NGz+1,NGx:N+NGx+1,NGy:NGy+L+1),Vhd,V1(NGz:M+NGz+1,NGx:N+NGx+1,NGy:NGy+L+1),&
			Whd,W1(NGz:M+Ngz,NGx:NGx+N+1,NGy:NGy+L+1),M,N,L,hx,hy,hz)
call Normal_viscVel(Ud,U3(NGz:M+NGz+1,NGx+1:N+NGx,NGy+1:NGy+L),Vd,V3(NGz:M+NGz+1,NGx+1:N+NGx,NGy+1:NGy+L),&
					Wd,W3(NGz:M+Ngz,NGx+1:NGx+N,NGy+1:NGy+L),M,N,L,hz)


call force(f,Bu1(3:M+NGz,3:NGx+N,3:NGy+L),Bs1(3:M+NGz,3:NGx+N,3:NGy+L),qvs,M,N,L)
call sponge_layer(U,U_spg,M,N,L)
call sponge_layer(W,W_spg,M,N,L)

U=U1(NGz+1:NGz+M,NGx+1:NGx+N,NGy+1:NGy+L)-dt*Ua+dt*Ud+dt*Uhd!+dt*U_spg!
V=V1(NGz+1:NGz+M,NGx+1:NGx+N,NGy+1:NGy+L)-dt*Va+dt*Vd+dt*Vhd
W=W1(NGz+1:NGz+M-1,NGx+1:NGx+N,NGy+1:NGy+L)-dt*Wa+dt*f+dt*Wd+dt*Whd!+dt*W_spg!


end subroutine vel_rk_step

subroutine B_rk_time_step(Bu,Bu1,Bu2,Bu3,Bs,Bs1,Bs2,Bs3,U1,V1,W1,M,N,L,dt,hx,hy,hz)

implicit none
real*8, dimension(M,N,L) :: Bu,Bs
real*8, dimension(M+4,N+4,L+4) :: Bu1,Bs1,Bu2,Bs2,Bu3,Bs3
real*8, dimension(M+4,N+4,L+4) :: U1,V1
real*8, dimension(M+3,N+4,L+4) :: W1
real*8, dimension(M,N,L) :: Bua,Bud,Bsa,Bsd,Wavg,Uavg,ls_fbu,ls_fbs,Buhd,Bshd,Bu_spg,Bs_spg
real*8, dimension(M+4,N+4,L+4) :: ql
real*8, dimension(M,N,L) :: qlz
integer :: M,N,L
real*8 :: hx,hz,hy,dt
integer :: i,j,k

call advectB(Bua,Bsa,Bu2,Bs2,U1,V1,W1,M,N,L,hx,hy,hz)
call diffB(Buhd,Bu1(Ngz:NGz+M+1,NGx:NGx+N+1,NGy:NGy+L+1),Bshd,Bs1(Ngz:NGz+M+1,NGx:NGx+N+1,NGy:NGy+L+1),M,N,L,hx,hy,hz)
call Normal_viscB(Bud,Bu3(Ngz:NGz+M+1,NGx:NGx+N+1,NGy:NGy+L+1),Bsd,Bs3(Ngz:NGz+M+1,NGx:NGx+N+1,NGy:NGy+L+1),M,N,L,hz)
call ls_force_therm(ls_fBu,ls_fBs,Bu1(NGz+1:NGz+M,NGx+1:NGx+N,NGy+1:NGy+L),Bs1(NGz+1:NGz+M,NGx+1:NGx+N,NGy+1:NGy+L),M,N,L)
!call sponge_layer(Bu,Bu_spg,M,N)
!call sponge_layer(Bs,Bs_spg,M,N)

do i=NGz+1,M+NGz
	do j=NGx+1,N+NGx
		do k=NGy+1,L+Ngy
			Wavg(i-2,j-2,k-2)=0.5*(W1(i-1,j,k)+W1(i,j,k))
		enddo
	enddo
enddo


!call calc_ql(ql(3:M+2,3:N+2),Bs,qvs,M,N)

Bu=Bu1(NGz+1:NGz+M,NGx+1:NGx+N,NGy+1:NGy+L)-dt*Bua+dt*ls_fbu+dt*Bud+dt*Buhd!+dt*Bu_spg!-Wavg*dthetae*dt+dt*ls_fbu!+dt*Bud
Bs=Bs1(NGz+1:NGz+M,NGx+1:NGx+N,NGy+1:NGy+L)-dt*Bsa+dt*ls_fbs+dt*Bsd+dt*Bshd!+dt*Bs_spg!+dt*Bshd+dt*Bsd!-Wavg*dqt*dt+dt*ls_fbs!+dt*Bsd

do i=1,M
	do j=1,N
		do k=1,L
			Bu(i,j,k)=Bu(i,j,k)-dt*Wavg(i,j,k)*dtheta_e_vec(i)
			Bs(i,j,k)=Bs(i,j,k)-dt*Wavg(i,j,k)*dqt_vec(i)
		enddo	
	enddo
enddo

!ql(1,3:N+2)=ql(4,3:N+2)
!ql(2,3:N+2)=ql(3,3:N+2)
!ql(M+3,3:N+2)=ql(M+2,3:N+2)
!ql(M+4,3:N+2)=ql(M+1,3:N+2)

!do i=1,N
!	ql(M+3,2+i)=max(Bs1(NGz+M+1,NGx+i)-qvs_func(zt(N)+hz),0.0)
!	ql(M+4,2+i)=max(Bs1(NGz+M+2,Ngx+i)-qvs_func(zt(N)+2.0*hz),0.0)
!enddo


!ql(:,1)=ql(:,N+1)
!ql(:,2)=ql(:,N+2)
!ql(:,N+3)=ql(:,3)
!ql(:,N+4)=ql(:,4)

!do i=3,M+2
!	do j=3,N+2
!		qlz(i-2,j-2)=(1.0/(12.0*hz))*(-ql(i+2,j)+8.0*(ql(i+1,j)-ql(i-1,j))+ql(i-2,j)) &
!					+(1.0/(12.0*hz))*sign(dble(1.0),Wavg(i-2,j-2))*(ql(i+2,j)-4.0*ql(i+1,j)+6.0*ql(i,j)-4.0*ql(i-1,j)+ql(i-2,j))
!	enddo
!enddo


Bu=Bu!-VT*dt*(Lv/cp)*qlz+dt*sf
Bs=Bs!+VT*dt*qlz


end subroutine B_rk_time_step




end module rk_time_step


