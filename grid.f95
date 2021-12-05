module grid

implicit none

real*8, dimension(:,:,:,:), allocatable :: var
real*8, dimension(:,:,:), allocatable :: ql
real*8, dimension(:,:,:), allocatable :: spge,sf
real*8, dimension(:), allocatable :: xm,xt,zm,zt,ym,yt
real*8, dimension(:), allocatable :: qvs
real*8, dimension(:), allocatable :: qt_vec,dqt_vec,dtheta_vec,dtheta_e_vec
real*8, dimension(:), allocatable :: write_freq_vec
real*8, dimension(:,:,:,:), allocatable :: write_array
integer :: x_dim,y_dim !dimensions of the grids
integer :: N_vert=40 !number of coarsest vertical momentum points.
!integer, parameter :: nls=150 !max number of vertical levels
!integer, dimension(nls) :: ref_levels !vertical levels at which the grid is refined
integer :: rat=2 !grid refinement ratio (always 2)
integer :: N_x=60 !number of coarsest grid horizontal cells
integer :: N_y=60
integer :: N_scalars=6 !number of variables, including pressure
real*8 :: hz=100 !vertical grid spacing
real*8 :: hx=100 !horizontal grid spacing
real*8 :: hy=100
integer :: NGx=2 !number of ghost cells in horizontal direction
integer :: NGy=2
integer :: NGz=2 !number of ghost cells in vertical direction
character(len=305) :: file_name='grid1'
real*8 :: Tend=15.0 !stopping time
real*8 :: dt=2.0 !coarsest grid time step
real*8 :: eb=0.0!1 !buoyancy diffusion coefficient
real*8 :: ev=0.0!1 !velocity diffusion coefficient
real*8 :: mu1=0.0!.25
real*8 :: mu2=3.0!2.0!.25 !horizontal viscosity coefficient
real*8 :: dtheta=3.0/1000.0 !Note that for grid14, this was 2/1000
real*8 :: dthetae
real*8 :: dqt=(-20.0/15000.0)*10.0**(-3)  !background total water gradient
real*8 :: qt0=22.0*10.0**(-3)
real*8 :: VT=2.0
integer :: qt_tilde_type=1
real*8 :: th0=300.0 !constant background potential temperature
real*8 :: max_moist=10.0*10.0**(-3) !for test_cold_bubble this was 0.15
real*8 :: max_moist_z=0.10 !for test_cold_bubble this was 0.15
real*8 :: max_t=-10.0!-12.0 !originally -2 !I think for grid13 this was -20
real*8 :: max_t_z=0.1!was originally 0.1
real*8 :: spge_h=0.6
real*8 :: BL_top=0.0 !Height of top of lower boundary layer
real*8 :: cp=1000
real*8 :: Rvd=0.61
real*8 :: Rv=461.5
real*8 :: Rd=287.1
real*8 :: Lv=2500.0*10.0**3
real*8 :: g=9.8
real*8 :: relax_u=86400.0 !relaxation parameter for horizontal velocity
real*8 :: therm_relax=86400.0 
real*8 :: write_freq=100.0
real*8 :: write_freq_dt=0.000001

contains


subroutine define_name_vars

implicit none

namelist /model/ &
	N_x,N_y,N_vert,hx,hy,hz,dt,Tend,write_freq,file_name
open(unit=1,file='NAMELIST')
read(1,nml=model)
 close(1)

end subroutine define_name_vars

end module grid

