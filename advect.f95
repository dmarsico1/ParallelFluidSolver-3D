module advect

use grid

contains

subroutine advectB(Bua,Bsa,Bu,Bs,U,V,W,M,N,L,hx,hy,hz)

implicit none

real*8, dimension(:,:,:) :: Bu,Bs,U,V
real*8, dimension(:,:,:) :: W
real*8, intent(out) :: Bua(M,N,L),Bsa(M,N,L)
real*8 , dimension(M,N,L) :: U_Bux,W_Buz,U_Bsx,W_Bsz,V_Buy,V_Bsy
integer, intent(in) :: M,N,L
real*8, intent(in) :: hx,hz,hy
real*8 :: Uavg,Wavg,Vavg
integer :: i,j,k

do j=3,N+2
	do i=3,M+2
		do k=3,L+2
			Uavg=0.5*(U(i,j+1,k)+U(i,j,k))
			U_Bux(i-2,j-2,k-2)=Uavg*(Bu(i,j+1,k)-Bu(i,j-1,k))/(2.0*hx)
			U_Bux(i-2,j-2,k-2)=Uavg*(1.0/(12.0*hx))*(-Bu(i,j+2,k)+8.0*(Bu(i,j+1,k)-Bu(i,j-1,k))+Bu(i,j-2,k)) &
						+(1.0/(12.0*hx))*abs(Uavg)*(Bu(i,j+2,k)-4.0*Bu(i,j+1,k)+6.0*Bu(i,j,k)-4.0*Bu(i,j-1,k)+Bu(i,j-2,k))
		enddo
	enddo
enddo

do j=3,N+2
	do i=3,M+2
		do k=3,L+2
			Vavg=0.5*(V(i,j,k+1)+V(i,j,k))
			V_Buy(i-2,j-2,k-2)=Vavg*(Bu(i,j,k+1)-Bu(i,j,k-1))/(2.0*hx)
			V_Buy(i-2,j-2,k-2)=Vavg*(1.0/(12.0*hy))*(-Bu(i,j,k+2)+8.0*(Bu(i,j,k+1)-Bu(i,j,k-1))+Bu(i,j,k-2)) &
						+(1.0/(12.0*hy))*abs(Vavg)*(Bu(i,j,k+2)-4.0*Bu(i,j,k+1)+6.0*Bu(i,j,k)-4.0*Bu(i,j,k-1)+Bu(i,j,k-2))
		enddo
	enddo
enddo


do j=3,N+2
	do i=3,M+2
		do k=3,L+2
			Wavg=0.5*(W(i-1,j,k)+W(i,j,k))
			W_Buz(i-2,j-2,k-2)=Wavg*(Bu(i+1,j,k)-Bu(i-1,j,k))/(2.0*hz)
			W_Buz(i-2,j-2,k-2)=Wavg*(1.0/(12.0*hz))*(-Bu(i+2,j,k)+8.0*(Bu(i+1,j,k)-Bu(i-1,j,k))+Bu(i-2,j,k)) &
						+(1.0/(12.0*hz))*abs(Wavg)*(Bu(i+2,j,k)-4.0*Bu(i+1,j,k)+6.0*Bu(i,j,k)-4.0*Bu(i-1,j,k)+Bu(i-2,j,k))
		enddo
	enddo
enddo

do j=3,N+2
	do i=3,M+2
		do k=3,L+2
		Uavg=0.5*(U(i,j+1,k)+U(i,j,k))
		U_Bsx(i-2,j-2,k-2)=Uavg*(Bs(i,j+1,k)-Bs(i,j-1,k))/(2.0*hx)
		U_Bsx(i-2,j-2,k-2)=Uavg*(1.0/(12.0*hx))*(-Bs(i,j+2,k)+8.0*(Bs(i,j+1,k)-Bs(i,j-1,k))+Bs(i,j-2,k)) &
					+(1.0/(12.0*hx))*abs(Uavg)*(Bs(i,j+2,k)-4.0*Bs(i,j+1,k)+6.0*Bs(i,j,k)-4.0*Bs(i,j-1,k)+Bs(i,j-2,k))
		enddo
	enddo
enddo

do j=3,N+2
	do i=3,M+2
		do k=3,L+2
			Vavg=0.5*(V(i,j,k+1)+V(i,j,k))
			V_Bsy(i-2,j-2,k-2)=Vavg*(Bs(i,j,k+1)-Bs(i,j,k-1))/(2.0*hx)
			V_Bsy(i-2,j-2,k-2)=Vavg*(1.0/(12.0*hy))*(-Bs(i,j,k+2)+8.0*(Bs(i,j,k+1)-Bs(i,j,k-1))+Bs(i,j,k-2)) &
						+(1.0/(12.0*hy))*abs(Vavg)*(Bs(i,j,k+2)-4.0*Bs(i,j,k+1)+6.0*Bs(i,j,k)-4.0*Bs(i,j,k-1)+Bs(i,j,k-2))
		enddo
	enddo
enddo

do j=3,N+2
	do i=3,M+2
		do k=3,L+2
			Wavg=0.5*(W(i-1,j,k)+W(i,j,k))
			W_Bsz(i-2,j-2,k-2)=Wavg*(Bs(i+1,j,k)-Bs(i-1,j,k))/(2.0*hz)
			W_Bsz(i-2,j-2,k-2)=Wavg*(1.0/(12.0*hz))*(-Bs(i+2,j,k)+8.0*(Bs(i+1,j,k)-Bs(i-1,j,k))+Bs(i-2,j,k)) &
						+(1.0/(12.0*hz))*abs(Wavg)*(Bs(i+2,j,k)-4.0*Bs(i+1,j,k)+6.0*Bs(i,j,k)-4.0*Bs(i-1,j,k)+Bs(i-2,j,k))
		enddo
	enddo
enddo


Bua=U_Bux+V_Buy+W_Buz
Bsa=U_Bsx+V_Bsy+W_Bsz

end subroutine advectB

subroutine advectVel(Ua,Va,Wa,U,V,W,M,N,L,hx,hy,hz)

!U,W have the appropriate layers of ghost cells surruonding them

implicit none

real*8, dimension(M+4,N+4,L+4) :: U,V
real*8, dimension(M+3,N+4,L+4) :: W
real*8, intent(out) :: Ua(M,N,L),Wa(M-1,N,L),Va(M,N,L)
real*8 :: U_Ux(M,N,L),V_Uy(M,N,L),W_Uz(M,N,L),U_Wx(M-1,N,L),V_Wy(M-1,N,L),W_Wz(M-1,N,L)
real*8 :: U_Vx(M,N,L),V_Vy(M,N,L),W_Vz(M,N,L)
real*8 :: hx,hz,hy
integer :: M,N,L
integer :: i,j,k
real*8 :: Wbar,Ubar,Vbar


do j=3,N+2
	do i=3,M+2
		do k=3,L+2
			U_Ux(i-2,j-2,k-2)=U(i,j,k)*(U(i,j+1,k)-U(i,j-1,k))/(2.0*hx)
			U_Ux(i-2,j-2,k-2)=U(i,j,k)*(1.0/(12.0*hx))*(-U(i,j+2,k)+8.0*(U(i,j+1,k)-U(i,j-1,k))+U(i,j-2,k)) &
						+(1.0/(12.0*hx))*abs(U(i,j,k))*(U(i,j+2,k)-4.0*U(i,j+1,k)+6.0*U(i,j,k)-4.0*U(i,j-1,k)+U(i,j-2,k))
		enddo
	enddo
enddo

do j=3,N+2
	do i=3,M+2
		do k=3,L+2
			Vbar=(1.0/4.0)*(V(i,j,k)+V(i,j-1,k)+V(i,j,k+1)+V(i,j-1,k+1))
			V_Uy(i-2,j-2,k-2)=Vbar*(U(i,j,k+1)-U(i,j,k-1))/(2.0*hy)
			V_Uy(i-2,j-2,k-2)=Vbar*(1.0/(12.0*hy))*(-U(i,j,k+2)+8.0*(U(i,j,k+1)-U(i,j,k-1))+U(i,j,k-2)) &
						+(1.0/(12.0*hy))*abs(Vbar)*(U(i,j,k+2)-4.0*U(i,j,k+1)+6.0*U(i,j,k)-4.0*U(i,j,k-1)+U(i,j,k-2))
		enddo
	enddo
enddo

do j=3,N+2
	do i=3,M+2
		do k=3,L+2
			Wbar=(1.0/4.0)*(W(i-1,j,k)+W(i-1,j-1,k)+W(i,j,k)+W(i,j-1,k))
			W_Uz(i-2,j-2,k-2)=Wbar*(U(i+1,j,k)-U(i-1,j,k))/(2.0*hz)
			W_Uz(i-2,j-2,k-2)=Wbar*(1.0/(12.0*hz))*(-U(i+2,j,k)+8.0*(U(i+1,j,k)-U(i-1,j,k))+U(i-2,j,k)) &
						+(1.0/(12.0*hz))*abs(Wbar)*(U(i+2,j,k)-4.0*U(i+1,j,k)+6.0*U(i,j,k)-4.0*U(i-1,j,k)+U(i-2,j,k))
		enddo
	enddo
enddo


!V HERE

do j=3,N+2
	do i=3,M+2
		do k=3,L+2
			Ubar=(1.0/4.0)*(U(i,j,k)+U(i,j+1,k)+U(i,j,k-1)+U(i,j+1,k-1))
			U_Vx(i-2,j-2,k-2)=Ubar*(V(i,j,k+1)-V(i,j,k-1))/(2.0*hy)
			U_Vx(i-2,j-2,k-2)=Ubar*(1.0/(12.0*hx))*(-V(i,j+2,k)+8.0*(V(i,j+1,k)-V(i,j-1,k))+V(i,j-2,k)) &
						+(1.0/(12.0*hx))*abs(Ubar)*(V(i,j+2,k)-4.0*V(i,j+1,k)+6.0*V(i,j,k)-4.0*V(i,j-1,k)+V(i,j-2,k))
		enddo
	enddo
enddo

do j=3,N+2
	do i=3,M+2
		do k=3,L+2
			Vbar=V(i,j,k)
			V_Vy(i-2,j-2,k-2)=Vbar*(V(i,j,k+1)-V(i,j,k-1))/(2.0*hy)
			V_Vy(i-2,j-2,k-2)=V(i,j,k)*(1.0/(12.0*hy))*(-V(i,j,k+2)+8.0*(V(i,j,k+1)-V(i,j,k-1))+V(i,j,k-2)) &
						+(1.0/(12.0*hy))*abs(V(i,j,k))*(V(i,j,k+2)-4.0*V(i,j,k+1)+6.0*V(i,j,k)-4.0*V(i,j,k-1)+V(i,j,k-2))
		enddo
	enddo
enddo

do j=3,N+2
	do i=3,M+2
		do k=3,L+2
			Wbar=(1.0/4.0)*(W(i,j,k)+W(i-1,j,k)+W(i,j,k-1)+W(i-1,j,k-1))
			W_Vz(i-2,j-2,k-2)=Wbar*(V(i+1,j,k)-V(i-1,j,k))/(2.0*hz)
			W_Vz(i-2,j-2,k-2)=Wbar*(1.0/(12.0*hz))*(-V(i+2,j,k)+8.0*(V(i+1,j,k)-V(i-1,j,k))+V(i-2,j,k)) &
						+(1.0/(12.0*hz))*abs(Wbar)*(V(i+2,j,k)-4.0*V(i+1,j,k)+6.0*V(i,j,k)-4.0*V(i-1,j,k)+V(i-2,j,k))
		enddo
	enddo
enddo

!W HERE

do j=3,N+2
	do i=3,M+1
		do k=3,L+2
			Ubar=(1.0/4.0)*(U(i,j,k)+U(i,j+1,k)+U(i+1,j,k)+U(i+1,j+1,k))
			U_Wx(i-2,j-2,k-2)=Ubar*(W(i,j+1,k)-W(i,j-1,k))/(2.0*hx)
			U_Wx(i-2,j-2,k-2)=Ubar*(1.0/(12.0*hx))*(-W(i,j+2,k)+8.0*(W(i,j+1,k)-W(i,j-1,k))+W(i,j-2,k)) &
						+(1.0/(12.0*hx))*abs(Ubar)*(W(i,j+2,k)-4.0*W(i,j+1,k)+6.0*W(i,j,k)-4.0*W(i,j-1,k)+W(i,j-2,k))
		enddo
	enddo
enddo

do j=3,N+2
	do i=3,M+1
		do k=3,L+2
			Vbar=(1.0/4.0)*(V(i,j,k)+V(i,j,k+1)+V(i+1,j,k)+V(i+1,j,k+1))
			V_Wy(i-2,j-2,k-2)=Vbar*(1.0/(12.0*hy))*(-W(i,j,k+2)+8.0*(W(i,j,k+1)-W(i,j,k-1))+W(i,j,k-2)) &
							+(1.0/(12.0*hy))*abs(Vbar)*(W(i,j,k+2)-4.0*W(i,j,k+1)+6.0*W(i,j,k)-4.0*W(i,j,k-1)+W(i,j,k-2))
		enddo
	enddo
enddo

do j=3,N+2
	do i=3,M+1
		do k=3,L+2
			W_Wz(i-2,j-2,k-2)=0.5*(W(i,j+1,k)+W(i,j-1,k))*(W(i+1,j,k)-W(i-1,j,k))/(2.0*hz)
			W_Wz(i-2,j-2,k-2)=W(i,j,k)*(1.0/(12.0*hz))*(-W(i+2,j,k)+8.0*(W(i+1,j,k)-W(i-1,j,k))+W(i-2,j,k)) &
						+(1.0/(12.0*hz))*abs(W(i,j,k))*(W(i+2,j,k)-4.0*W(i+1,j,k)+6.0*W(i,j,k)-4.0*W(i-1,j,k)+W(i-2,j,k))
		enddo	
	enddo
enddo


Ua=U_Ux+V_Uy+W_Uz
Va=U_Vx+V_Vy+W_Vz
Wa=U_Wx+V_Wy+W_Wz

end subroutine advectVel


end module advect
