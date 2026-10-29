program project1

implicit none

integer, parameter :: dp = kind(1.0d0) ! double precision
integer, parameter :: jmax=5000 ! need jmax to be a known value at compile time cause we use it to declare array sizes later

! define constants from slides
real(dp), parameter :: msol = 1.989d33 ! g
real(dp), parameter :: kpc = 3.086d21 ! cm
real(dp), parameter :: G = 6.672d-8 ! cgs
real(dp), parameter :: k_B = 1.380649d-16 ! erg/K
real(dp), parameter :: m_p = 1.6726219d-24 ! g
real(dp), parameter :: mu = 0.61_dp
real(dp), parameter :: pi = 3.14159265358979323846_dp

real(dp) :: r(jmax), rr(jmax)
real(dp) :: rmin, rmax, dr
integer  :: j, iu

real(dp) :: vol(jmax)
real(dp) :: rhodm(jmax)
real(dp) :: mnfw(jmax), mdark(jmax)
real(dp) :: mhern(jmax)
real(dp) :: x

real(dp), parameter :: rho0dm = 7.35d-26  !Navarro-Frank-White Dark Matter density
real(dp), parameter :: rs = 435.7_dp*kpc
real(dp), parameter :: rho0 = 2.882d-26
real(dp), parameter :: ticm = 8.9d7       !Temp in ICM
real(dp), parameter :: rvir = 2797.*kpc   !Virial radius
real(dp), parameter :: fc = 1.138799_dp     !Normalization factor
real(dp), parameter :: mvir = 1.3d15*msol !Virial mass
real(dp), parameter :: mbgc = 1.d12*msol
real(dp), parameter :: ahern=10.*kpc/(1.+2.**(0.5))

! set the grid

rmin = 0.0_dp*kpc
rmax = 3000.0_dp*kpc
dr = (rmax-rmin)/real(jmax-1, dp) ! define grid spacing between nodes: the total length divided by the number of intervals
do j=1,jmax
    r(j)  = rmin + real(j-1, dp)*dr ! nodes r_j
    rr(j) = r(j) + 0.5_dp*dr ! half-point grid, shifted by half a step, r_{j+1/2}, defined for every j
enddo

! check that the grid ends at rmax
if (abs(r(jmax) - rmax) > 1.0d-10*rmax) stop 'grid does not end at rmax' ! stop code if grid does not end at rmax

! write to file (in kpc)
open(newunit=iu, file='grid.dat', status='replace')
do j = 1, jmax
   write(iu,'(2(1pe14.6))') r(j)/kpc, rr(j)/kpc
end do
close(iu)

!! Density Profile Considering Dark Matter Dominance !!


vol(1) = 4.1888*r(1)**3 
do j=2, jmax
    vol(j) = 4.1888*(r(j)**3 - r(j-1)**3) !Volume between two shells
enddo

            !NFW Dark Matter Profile
do j=1,jmax                       
    x=rr(j)/rs
    rhodm(j)=rho0dm/(x*(1.+x)**2)
enddo

open(20,file='mass_profile.dat', status='replace')  !Store mass profile data into file/this overwrites file btw!

mnfw(1) = 0   !Navarro Frank White Mass, no mass in 1st
do j=2, jmax
    x = r(j)/rs
    mnfw(j) = mnfw(j-1) + rhodm(j-1)*vol(j)    !DM mass
    mdark(j) = mvir*(log(1.+x) - x/(1.+x))/fc  !Analytical DM mass
    mhern(j) = mbgc*r(j)**2/(r(j)+ahern)**2    !Analytical stellar mass
    write(20,100)r(j)/kpc,mnfw(j)/msol,mdark(j)/msol,mhern(j)/msol
enddo

100 format(4(1pe12.4))
close(20)

end program project1
