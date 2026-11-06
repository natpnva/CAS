program project1_part2
    implicit none

    ! Part II: turbulent diffusion of iron in a galaxy cluster
    ! First test: diffusion only, without iron sources
    ! Uses the gas density from Model C of Part I

    integer, parameter :: dp = kind(1.0d0)
    integer, parameter :: n = 5000
    integer, parameter :: nc = n - 1

    ! Physical constants in cgs units
    real(dp), parameter :: pi = 3.141592653589793_dp
    real(dp), parameter :: kpc = 3.086d21
    real(dp), parameter :: gyr = 3.15576d16

    ! Constant turbulent diffusion coefficient (cm^2/s)
    real(dp), parameter :: D = 1.32d29

    ! Radial grid and gas density
    real(dp) :: radius(n), rho_gas(n)
    real(dp) :: rnode(n), volume(nc), area(n)
    real(dp) :: dr, temp

    ! Iron abundance and diffusive mass flux
    real(dp) :: x, zfe(nc), znew(nc), zfe_observed
    real(dp) :: flux(n), rho_wall
    real(dp) :: mass_initial, mass_final
    real(dp), parameter :: abundance_factor = 1.4_dp
    real(dp), parameter :: zfe_sol = 0.0018_dp  
    real(dp), parameter :: zfe_outer = 0.4_dp * zfe_sol

    ! Time integration
    real(dp) :: time, dt, dt_stable, target_time
    real(dp) :: snapshots(4)
    integer :: j, k, iu, ios
    character(len=256) :: header
    character(len=40) :: filename

    ! Diffusion timescale
    real(dp) :: mass_100_initial, mass_100
    real(dp) :: diffusion_time
    logical :: half_mass_reached

    ! ---------------------------------------------------------
    ! 1. Read gas density from Part I
    ! ---------------------------------------------------------

    open(newunit=iu, file="part1.dat", status="old", action="read", iostat=ios)

    if (ios /= 0) stop "Cannot open part1.dat"

    read(iu,'(A)') header   ! Skip column headings

    do j = 1, n
        ! Read radius, then five other columns, then T
        ! The sixth column is Model C gas density
        read(iu,*,iostat=ios) radius(j), temp, temp, temp, temp, rho_gas(j), temp

        if (ios /= 0) stop "Error reading part1.dat"
    end do

    close(iu)

    ! ---------------------------------------------------------
    ! 2. Reconstruct radial grid and shell geometry
    ! ---------------------------------------------------------

    dr = (radius(2) - radius(1)) * kpc

    do j = 1, n
        rnode(j) = real(j-1,dp) * dr
        area(j) = 4.0_dp * pi * rnode(j)**2
    end do

    do j = 1, nc
        ! Volume of the spherical shell containing rho_gas(j)
        volume(j) = 4.0_dp * pi / 3.0_dp * (rnode(j+1)**3 - rnode(j)**3)
    end do

    ! ---------------------------------------------------------
    ! 3. Set initial iron abundance
    ! ---------------------------------------------------------

    do j = 1, nc
        ! Gaussian test profile, concentrated near cluster center
        ! zfe is a dimensionless iron mass fraction
        ! zfe(j) = 1.0d-3 * exp(-0.5_dp * (radius(j) / 30.0_dp)**2)

        x = radius(j)/80.0_dp
        zfe_observed = zfe_sol * (0.3_dp)*(1.4_dp)*(1.15_dp)*(2.2_dp +x**3)/(1.0_dp + x**3)/1.15_dp 

        ! Subtract background iron abundance
        zfe_observed = max(zfe_observed - zfe_outer, 0.0_dp)

    !NEW: Calculate observed iron abundance at radius(j) using the observational profile from the slides
        zfe(j) = zfe_observed / abundance_factor
    end do

    ! Calculate initial total iron mass
    mass_initial = sum(rho_gas(1:nc) * zfe * volume) 

    mass_100_initial = sum(rho_gas(1:nc) * zfe * volume, mask=radius(1:nc)<100.0_dp)
    half_mass_reached = .false.
    diffusion_time = -1.0_dp

    print *, "Central iron mass fraction:", zfe(1)
    print *, "Outer iron mass fraction:", zfe(nc)
    print *, "Initial iron mass [g]:", mass_initial
    print *, "Initial M_Fe(<100 kpc) [g]:", mass_100_initial

    ! ---------------------------------------------------------
    ! 4. Set time step and output times
    ! ---------------------------------------------------------

    ! Conservative explicit diffusion time step
    dt_stable = 0.2_dp*dr**2/D

    snapshots = [0.0_dp, 1.0_dp, 2.0_dp, 5.0_dp]*gyr

    time = 0.0_dp

    print *, "Grid spacing [kpc]:", dr/kpc
    print *, "Time step [Myr]:", dt_stable/gyr*1000.0_dp

    ! Save initial abundance profile
    call save_profile("iron_0gyr.dat")

    ! ---------------------------------------------------------
    ! 5. Evolve iron diffusion
    ! ---------------------------------------------------------

    do k = 2, 4

        target_time = snapshots(k)

        do while (time < target_time)

            ! Shorten final step to hit snapshot time exactly
            dt = min(dt_stable, target_time-time)

            ! No-flux boundary conditions at both ends
            flux(1) = 0.0_dp
            flux(n) = 0.0_dp

            ! Calculate diffusive iron flux through each
            ! internal spherical shell boundary
            do j = 2, nc

                ! Average gas density on the boundary
                rho_wall = 0.5_dp*(rho_gas(j-1)+rho_gas(j))

                ! Fick's law: iron flows down abundance gradient
                flux(j) = -D*rho_wall * (zfe(j)-zfe(j-1))/dr
            end do

            ! Update iron abundance in each shell
            do j = 1, nc

                ! Net iron mass flowing into/out of the shell
                znew(j) = zfe(j) - dt * (area(j+1)*flux(j+1)-area(j)*flux(j)) / (rho_gas(j)*volume(j))
            end do

            ! Move to the next time step
            zfe = znew
            time = time + dt

            if (.not. half_mass_reached) then

                mass_100 = sum(rho_gas(1:nc)*zfe*volume, mask=radius(1:nc)<100.0_dp)

            if (mass_100 <= 0.5_dp * mass_100_initial) then
                diffusion_time = time
                half_mass_reached = .true.

                print *, "Diffusion timescale [Gyr]:", diffusion_time/gyr
            end if

end if

        end do

        ! Save iron abundance at requested time
        select case(k)
        case(2)
            filename = "iron_1gyr.dat"
        case(3)
            filename = "iron_2gyr.dat"
        case(4)
            filename = "iron_5gyr.dat"
        end select

        call save_profile(trim(filename))

        ! Check iron mass conservation
        mass_final = sum(rho_gas(1:nc)*zfe*volume)

        print *, "M_Fe(<100 kpc) [g]:", sum(rho_gas(1:nc)*zfe*volume, mask=radius(1:nc)<100.0_dp)
        print *, "Time [Gyr]:", time/gyr
        print *, "Relative iron mass change:", (mass_final-mass_initial)/mass_initial

    end do

    if (.not. half_mass_reached) then
        print *, "Iron mass within 100 kpc did not decrease by 50% within 5 Gyr."
    end if

    print *, "Part II diffusion test completed."

contains

    ! ---------------------------------------------------------
    ! Save radius and iron abundance for Python plotting
    ! ---------------------------------------------------------
    subroutine save_profile(fname)
        character(len=*), intent(in) :: fname
        integer :: unit, i

        open(newunit=unit, file=fname, status="replace", action="write")

        write(unit,'(A)') "# Radius [kpc]   Iron mass fraction"

        do i = 1, nc
            write(unit,'(2(1pe18.9))') radius(i), zfe(i)
        end do

        close(unit)
    end subroutine save_profile

end program project1_part2
