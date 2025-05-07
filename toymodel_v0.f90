!======================================================================!
! Program: toymodel_v0
! Purpose: Minimal 1-layer, 1-pool SOM model with variable timestep
!          and correct mass conservation check before state update
!======================================================================!
program toymodel_v0
    implicit none

    !------------------------------------------------------------------!
    ! Precision definition
    !------------------------------------------------------------------!
    integer, parameter :: dp = selected_real_kind(15, 307)             ! Double precision

    !------------------------------------------------------------------!
    ! Simulation control parameters
    !------------------------------------------------------------------!
    integer, parameter 	:: nyr = 6000                                  ! Total simulation years
    real(dp), parameter :: dt = 1.0_dp                                 ! Timestep length (fraction of year; e.g., 0.25 = quarterly)
    integer, parameter 	:: nsteps = nint(1.0_dp / dt)                  ! Number of steps per year
    real(dp), parameter :: eps = 1.0d-8                                ! Mass conservation tolerance (kg C/m2)

    !------------------------------------------------------------------!
    ! Model parameters
    !------------------------------------------------------------------!
    real(dp), parameter :: input_rate = 1.05_dp                        ! Annual litter input (kg C/m2/year)
    real(dp), parameter :: k_decay = 0.007_dp                          ! Annual decay rate of SOM (/year)

    !------------------------------------------------------------------!
    ! State variables
    !------------------------------------------------------------------!
    real(dp) :: SOM                                                    ! Soil organic matter pool (kg C/m2)
    real(dp) :: SOM_old                                                ! SOM at start of timestep (kg C/m2)
    real(dp) :: litter, resp, dSOM                                     ! Fluxes and net change (kg C/m2)
    real(dp) :: mass_start, mass_end, mass_error                       ! Mass conservation check (kg C/m2)

    !------------------------------------------------------------------!
    ! Loop counters
    !------------------------------------------------------------------!
    integer :: kyr, it                                                 ! kyr = year, it = timestep within year
    real(dp) :: time                                                   ! Continuous simulation time (years)

    !------------------------------------------------------------------!
    ! Initialization
    !------------------------------------------------------------------!
    SOM = 0.0_dp                                                       ! Initial SOM (kg C/m2)
    write(*,'(a)') 'Year    SOM_C     Input    Respired'               ! Output header

    !------------------------------------------------------------------!
    ! Simulation loop over years and sub-timesteps
    !------------------------------------------------------------------!
    do kyr = 1, nyr
        do it = 1, nsteps

            time = (kyr - 1) + it * dt                                 ! Current simulation time (years)
            SOM_old = SOM                                              ! Save SOM before update

            !----------------------------------------------------------!
            ! Step 1: Compute gross fluxes from SOM_old
            !----------------------------------------------------------!
            litter = input_rate * dt                                   ! Litter input (kg C/m2)
            resp   = k_decay * SOM_old * dt                            ! Respiration loss (kg C/m2)
            dSOM   = litter - resp                                     ! Net change in SOM (kg C/m2)

            !----------------------------------------------------------!
            ! Step 2: Mass conservation check BEFORE state update
            ! mass_start = SOM before + input
            ! mass_end   = new SOM + respired carbon
            !----------------------------------------------------------!
            mass_start = SOM_old + litter
            mass_end   = SOM_old + dSOM + resp                         ! SOM is NOT yet updated
            mass_error = abs(mass_end - mass_start)

            if (mass_error > eps) then
                write(*,'(a,f6.2)') 'Mass conservation error at: ', time
                write(*,'(a,f12.5)') 'Start mass = ', mass_start
                write(*,'(a,f12.5)') 'End mass   = ', mass_end
                write(*,'(a,f12.5)') 'Difference = ', mass_error
                stop 'Mass not conserved'
            end if

            !----------------------------------------------------------!
            ! Step 3: Update SOM only after mass check passes
            !----------------------------------------------------------!
            SOM = SOM_old + dSOM                                       ! Update SOM (kg C/m2)

        end do                                                         ! End of sub-timesteps

        !--------------------------------------------------------------!
        ! Output annual status of SOM and fluxes
        !--------------------------------------------------------------!
        write(*,'(i5,3f12.5)') kyr, SOM, litter, resp
    end do                                                             ! End of yearly loop

end program toymodel_v0
!======================================================================!
