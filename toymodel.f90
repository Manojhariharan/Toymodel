!======================================================================!
! Program: toymodel Model
! Purpose: Simulates redistribution of SOM and mineral soil mass 
!          across fixed-depth layers while maintaining constant 
!          interface depths (except for the bottom).
!======================================================================!

program toymodel
    implicit none

    !------------------------------------------------------------------!
    ! Parameters and constants
    !------------------------------------------------------------------!
    integer, parameter  :: nlayers = 3                   ! Number of soil layers
    real, parameter     :: rho_SOM = 1.0                 ! SOM density (g/cm3)
    real, parameter     :: rho_min = 1.5                 ! Mineral soil density (g/cm3)
    integer, parameter  :: nsteps = 100                  ! Number of time steps
    real, parameter     :: eps = 1.0e-5                  ! Tolerance for conservation check

    !------------------------------------------------------------------!
    ! Variables
    !------------------------------------------------------------------!
    real, dimension(nlayers)     :: MSOM                 ! SOM mass per layer (g/cm2)
    real, dimension(nlayers)     :: Mmin                 ! Mineral soil mass per layer (g/cm2)
    real, dimension(nlayers)     :: dMSOM                ! SOM mass increment per layer (g/cm2)
    real, dimension(nlayers)     :: dMmin                ! Mineral mass increment per layer (g/cm2)
    real, dimension(nlayers+1)   :: z                    ! Interface depths (cm)
    real, dimension(nlayers+1)   :: z_initial            ! To store initial depths

    real :: dz, dz_p                                     ! Actual and potential layer thicknesses (cm)
    real :: f, mMSOM, mMmin                              ! Redistribution fractions and masses
    real :: total_SOM_before, total_Mmin_before          ! Mass totals before update
    real :: total_SOM_after, total_Mmin_after            ! Mass totals after update
    real :: som_diff, min_diff                           ! Difference in mass before and after each timestep (g/cm2)
    real :: total_SOM_initial, total_Mmin_initial        ! Mass at the beginning of the simulation (g/cm2)
    real :: total_SOM_final, total_Mmin_final            ! Mass at the end of the simulation (g/cm2)
    logical :: som_ok, min_ok                            ! Flag indicating whether mass was conserved
    integer :: i, j                                      ! Loop indices

    !------------------------------------------------------------------!
    ! Initialization
    !------------------------------------------------------------------!
    MSOM = 1.0                                           ! Initial SOM mass in each layer (g/cm2)
    Mmin = 1.0                                           ! Initial mineral soil mass in each layer (g/cm2)
    dMSOM = 0.0                                          ! Initialize SOM mass additions per layer to zero (g/cm2/timestep)
    dMmin = 0.0                                          ! Initialize mineral mass additions per layer to zero (g/cm2/timestep)
    dMSOM(1) = 1.0                                       ! Add 1 g SOM to the top layer per time step (g/cm2/timestep)

    total_SOM_initial = sum(MSOM)
    total_Mmin_initial = sum(Mmin)

    ! Compute initial layer interface depths
    z(1) = 0.0
    do i = 2, nlayers
        z(i) = z(i-1) + MSOM(i)/rho_SOM + Mmin(i)/rho_min
    end do
    z(nlayers+1) = sum(MSOM)/rho_SOM + sum(Mmin)/rho_min
    z_initial = z

    som_ok = .true.
    min_ok = .true.

    !------------------------------------------------------------------!
    ! Time step simulation loop
    !------------------------------------------------------------------!
    do j = 1, nsteps

        ! Track initial total mass for conservation check
        total_SOM_before = sum(MSOM) + sum(dMSOM)
        total_Mmin_before = sum(Mmin) + sum(dMmin)

        ! Apply additions and redistribute SOM and mineral mass
        do i = 1, nlayers - 1

            ! Apply input additions
            MSOM(i) = MSOM(i) + dMSOM(i)
            Mmin(i) = Mmin(i) + dMmin(i)

            ! Calculate potential and actual thickness
            dz_p = MSOM(i)/rho_SOM + Mmin(i)/rho_min
            dz = z(i+1) - z(i)

            ! Compute overflow fraction if potential exceeds actual
            f = (dz_p - dz) / dz_p

            ! Determine mass to redistribute
            mMSOM = f * MSOM(i)
            mMmin = f * Mmin(i)

            ! Update current layer
            MSOM(i) = MSOM(i) - mMSOM
            Mmin(i) = Mmin(i) - mMmin

            ! Transfer overflow to next layer
            MSOM(i+1) = MSOM(i+1) + mMSOM
            Mmin(i+1) = Mmin(i+1) + mMmin
        end do

        ! Update total profile depth (only bottom interface changes)
        z(nlayers+1) = sum(MSOM)/rho_SOM + sum(Mmin)/rho_min

        ! Track final total mass for conservation check
        total_SOM_after = sum(MSOM)
        total_Mmin_after = sum(Mmin)

        som_diff = abs(total_SOM_after - total_SOM_before)
        min_diff = abs(total_Mmin_after - total_Mmin_before)

        if (som_diff > eps) then
            write (*,'(a,i3,a,f10.6)') 'WARNING: SOM mass not conserved at step ', j, &
                                        ', Δ = ', total_SOM_after - total_SOM_before
            som_ok = .false.
        end if

        if (min_diff > eps) then
            write (*,'(a,i3,a,f10.6)') 'WARNING: Mineral mass not conserved at step ', j, &
                                        ', Δ = ', total_Mmin_after - total_Mmin_before
            min_ok = .false.
        end if

    end do

    ! Final totals
    total_SOM_final = sum(MSOM)
    total_Mmin_final = sum(Mmin)

    !------------------------------------------------------------------!
    ! Final output summary
    !------------------------------------------------------------------!
    write(*,*) '==================================================='
    write(*,*) '                  FINAL SUMMARY                    '
    write(*,*) '==================================================='
    write(*,*) 'Initial interface depths (z):'
    write(*,'(a,*(f8.2))') ' z_initial = ', z_initial

    write(*,*) 'Final interface depths (z):'
    write(*,'(a,*(f8.2))') ' z_final   = ', z

    write(*,'(a,f10.4)') ' Initial total SOM mass   = ', total_SOM_initial
    write(*,'(a,f10.4)') ' Final total SOM mass     = ', total_SOM_final
    write(*,'(a,f10.4)') ' Initial total Mineral    = ', total_Mmin_initial
    write(*,'(a,f10.4)') ' Final total Mineral      = ', total_Mmin_final

    if (som_ok) then
        write(*,*) 'SUCCESS: SOM mass conserved across all timesteps.'
    else
        write(*,*) 'FAILURE: SOM mass not conserved in some timesteps.'
    end if

    if (min_ok) then
        write(*,*) 'SUCCESS: Mineral mass conserved across all timesteps.'
    else
        write(*,*) 'FAILURE: Mineral mass not conserved in some timesteps.'
    end if
    write(*,*) '==================================================='

end program toymodel
!======================================================================!
