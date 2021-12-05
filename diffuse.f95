module diffuse

use grid

contains

subroutine diffVel(Ud,U,Vd,V,Wd,W,M,N,L,hx,hy,hz)

!This subroutine assumes that U,W already have the appropriate layers of ghost cells

implicit none

real*8, intent(out) :: Ud(M,N,L),Wd(M-1,N,L),Vd(M,N,L)
real*8, intent(inout) :: W(:,:,:),U(:,:,:),V(:,:,:)
real*8, intent(in) :: hx,hz,hy
integer, intent(in) :: M,N,L
integer :: i,j,k

do i=1,M
	do j=1,N
		do k=1,L
			Ud(i,j,k)=mu2*(1.0/hx**2)*(U(i+1,j+2,k+1)-2.0*U(i+1,j+1,k+1)+U(i+1,j,k+1)) &
					+mu2*(1.0/hy**2)*(U(i+1,j+1,k+2)-2.0*U(i+1,j+1,k+1)+U(i+1,j+1,k)) !&
 					!+mu1*(1.0/hz**2)*(U(i+2,j+1,k+1)-2.0*U(i+1,j+1,k+1)+U(i,j+1,k+1))
		enddo
	enddo
enddo

do i=1,M
	do j=1,N
		do k=1,L
			Vd(i,j,k)=mu2*(1.0/hx**2)*(V(i+1,j+2,k+1)-2.0*V(i+1,j+1,k+1)+V(i+1,j,k+1)) &
					+mu2*(1.0/hy**2)*(V(i+1,j+1,k+2)-2.0*V(i+1,j+1,k+1)+V(i+1,j+1,k))! &
 					!+mu1*(1.0/hz**2)*(V(i+2,j+1,k+1)-2.0*V(i+1,j+1,k+1)+V(i,j+1,k+1))
		enddo
	enddo
enddo

do i=1,M-1
	do j=1,N
		do k=1,L
			Wd(i,j,k)=mu2*(1.0/hx**2)*(W(i+1,j+2,k+1)-2.0*W(i+1,j+1,k+1)+W(i+1,j,k+1)) &
					+mu2*(1.0/hy**2)*(W(i+1,j+1,k+2)-2.0*W(i+1,j+1,k)+W(i+1,j+1,k)) !&
 					!+mu1*(1.0/hz**2)*(W(i+2,j+1,k+1)-2.0*W(i+1,j+1,k+1)+W(i,j+1,k+1))
		enddo
	enddo
enddo

end subroutine diffVel


subroutine diffB(Bud,Bu,Bsd,Bs,M,N,L,hx,hy,hz)

!This subroutine already assumes that Bu,Bs already have the approprate layer of ghost cells

implicit none

real*8, intent(out) :: Bud(M,N,L),Bsd(M,N,L)
real*8, intent(inout) :: Bu(M+2,N+2,L+2),Bs(M+2,N+2,L+2)
real*8 :: hx,hz,hy
integer, intent(in) :: M,N,L
integer :: i,j,k

do i=1,M
	do j=1,N
		do k=1,L
			Bud(i,j,k)=mu2*(1.0/hx**2)*(Bu(i+1,j+2,k+1)-2.0*Bu(i+1,j+1,k+1)+Bu(i+1,j,k+1)) &
						+mu2*(1.0/hy**2)*(Bu(i+1,j+1,k+2)-2.0*Bu(i+1,j+1,k+1)+Bu(i+1,j+1,k))! &
 						!+eb*(1.0/hz**2)*(Bu(i+2,j+1,k+1)-2.0*Bu(i+1,j+1,k+1)+Bu(i,j+1,k+1))
			Bsd(i,j,k)=mu2*(1.0/hx**2)*(Bs(i+1,j+2,k+1)-2.0*Bs(i+1,j+1,k+1)+Bs(i+1,j,k+1)) &
						+mu2*(1.0/hy**2)*(Bs(i+1,j+1,k+2)-2.0*Bs(i+1,j+1,k+1)+Bs(i+1,j+1,k))! &
 						!+eb*(1.0/hz**2)*(Bs(i+2,j+1,k+1)-2.0*Bs(i+1,j+1,k+1)+Bs(i,j+1,k+1))
		enddo
	enddo
enddo

end subroutine diffB

subroutine Horiz_HyperDiff(U,V,W,Uhd,Vhd,Whd,M,N,L,hx)

implicit none

real*8, dimension(M,N,L) :: Uhd,Vhd
real*8, dimension(M,N+4,L+4) :: U,V
real*8, dimension(M-1,N+4,L+4) :: W
real*8, dimension(M-1,N,L) :: Whd
integer, intent(in) :: M,N,L
real*8 :: hx
integer :: i,j,k

do i=1,M
	do j=NGx+1,NGx+N
		do k=NGy+1,NGy+L
			Uhd(i,j-NGx,k-NGy)=-mu1*(1.0/hx**4)*(U(i,j-2,k)-4.0*U(i,j-1,k)+6.0*U(i,j,k)-4.0*U(i,j+1,k)+U(i,j+2,k)) &
									-mu1*(1.0/hy**4)*(U(i,j,k-2)-4.0*U(i,j,k-1)+6.0*U(i,j,k)-4.0*U(i,j,k+1)+U(i,j,k+2))
		enddo
	enddo
enddo

do i=1,M
	do j=NGx+1,NGx+N
		do k=NGy+1,NGy+L
			Vhd(i,j-NGx,k-NGy)=-mu1*(1.0/hx**4)*(V(i,j-2,k)-4.0*V(i,j-1,k)+6.0*V(i,j,k)-4.0*V(i,j+1,k)+V(i,j+2,k)) &
									-mu1*(1.0/hy**4)*(V(i,j,k-2)-4.0*V(i,j,k-1)+6.0*V(i,j,k)-4.0*V(i,j,k+1)+V(i,j,k+2))
		enddo
	enddo
enddo

do i=1,M-1
	do j=NGx+1,NGx+N
		do k=NGy+1,NGy+L
			Whd(i,j-NGx,k-NGy)=-mu1*(1.0/hx**4)*(W(i,j-2,k)-4.0*W(i,j-1,k)+6.0*W(i,j,k)-4.0*W(i,j+1,k)+W(i,j+2,k)) &
									-mu1*(1.0/hy**4)*(W(i,j,k-2)-4.0*W(i,j,k-1)+6.0*W(i,j,k)-4.0*W(i,j,k+1)+W(i,j,k+2))
		enddo
	enddo
enddo

end subroutine Horiz_HyperDiff


subroutine Horiz_HyperDiffB(Bu,Bs,Buhd,Bshd,M,N,L,hx)

real*8, dimension(M,N,L) :: Buhd,Bshd
real*8, dimension(M,N+4,L+4) :: Bu,Bs
integer, intent(in) :: M,N,L
real*8 :: hx
integer :: i,j,k

do i=1,M
	do j=NGx+1,NGx+N
		do k=NGy+1,NGy+L
			Buhd(i,j-NGx,k-NGy)=-mu1*(1.0/hx**4)*(Bu(i,j-2,k)-4.0*Bu(i,j-1,k)+6.0*Bu(i,j,k)-4.0*Bu(i,j+1,k)+Bu(i,j+2,k)) &
								-mu1*(1.0/hy**4)*(Bu(i,j,k-2)-4.0*Bu(i,j,k-1)+6.0*Bu(i,j,k)-4.0*Bu(i,j,k+1)+Bu(i,j,k+2))
			Bshd(i,j-NGx,k-NGy)=-mu1*(1.0/hx**4)*(Bs(i,j-2,k)-4.0*Bs(i,j-1,k)+6.0*Bs(i,j,k)-4.0*Bs(i,j+1,k)+Bs(i,j+2,k)) &
								-mu1*(1.0/hy**4)*(Bs(i,j,k-2)-4.0*Bs(i,j,k-1)+6.0*Bs(i,j,k)-4.0*Bs(i,j,k+1)+Bs(i,j,k+2))
		enddo
	enddo
enddo


end subroutine Horiz_HyperDiffB


subroutine Normal_viscVel(Ud,U,Vd,V,Wd,W,M,N,L,hz)

!This subroutine assumes that U,W already have the appropriate layers of ghost cells

implicit none

real*8, intent(out) :: Ud(M,N,L),Wd(M-1,N,L),Vd(M,N,L)
real*8, intent(inout) :: W(M+1,N,L),U(M+2,N,L),V(M+2,N,L)
real*8, intent(in) :: hz
integer, intent(in) :: M,N,L
integer :: i,j,k

do i=1,M
	do j=1,N
		do k=1,L
			Ud(i,j,k)=mu2*(1.0/hz**2)*(U(i+2,j,k)-2.0*U(i+1,j,k)+U(i,j,k))
		enddo
	enddo
enddo

do i=1,M
	do j=1,N
		do k=1,L
			Vd(i,j,k)=mu2*(1.0/hz**2)*(V(i+2,j,k)-2.0*V(i+1,j,k)+V(i,j,k))
		enddo
	enddo
enddo

do i=1,M-1
	do j=1,N
		do k=1,L
			Wd(i,j,k)=mu2*(1.0/hz**2)*(W(i+2,j,k)-2.0*W(i+1,j,k)+W(i,j,k))
		enddo
	enddo
enddo


end subroutine Normal_viscVel

subroutine Normal_viscB(Bud,Bu,Bsd,Bs,M,N,L,hz)

!This subroutine already assumes that Bu,Bs already have the approprate layer of ghost cells

implicit none

real*8, intent(out) :: Bud(M,N,L),Bsd(M,N,L)
real*8, intent(inout) :: Bu(M+2,N+2,L+2),Bs(M+2,N+2,L+2)
real*8 :: hz
integer, intent(in) :: M,N,L
integer :: i,j,k

do i=1,M
	do j=1,N
		do k=1,L
			Bud(i,j,k)=mu2*(1.0/hz**2)*(Bu(i+2,j+1,k+1)-2.0*Bu(i+1,j+1,k+1)+Bu(i,j+1,k+1))
			Bsd(i,j,k)=mu2*(1.0/hz**2)*(Bs(i+2,j+1,k+1)-2.0*Bs(i+1,j+1,k+1)+Bs(i,j+1,k+1))
		enddo
	enddo
enddo

end subroutine Normal_viscB

end module diffuse
