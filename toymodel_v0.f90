!======================================================================!
! Program: toymodel_v0
! Purpose: Minimal 1-layer, 1-pool SOM model with structured timestep,
!          correct dimensional units, and mass conservation check
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
    integer, parameter 	:: nyr = 6000                                  ! Total simulation duration (years)
    real(dp), parameter :: dt = 1.0_dp                                 ! Timestep (fraction of year; e.g., 0.25 = quarterly)
    integer, parameter 	:: nsteps = nint(1.0_dp / dt)                  ! Number of sub-steps per year
    real(dp), parameter :: eps = 1.0d-8                                ! Mass conservation tolerance (kg C/m2)

    !------------------------------------------------------------------!
    ! Model parameters
    !------------------------------------------------------------------!
    real(dp), parameter :: input_rate = 1.05_dp                        ! Annual litter input rate (kg C/m2/year)
    real(dp), parameter :: k_decay = 0.007_dp                          ! Annual decay rate of SOM (/year)

    !------------------------------------------------------------------!
    ! State and flux variables
    !------------------------------------------------------------------!
    real(dp) :: SOM                                                    ! Soil organic matter pool (kg C/m2)
    real(dp) :: SOM_before, SOM_after                                  ! SOM before and after update (kg C/m2)
    real(dp) :: dSOM                                                   ! Net SOM rate of change (kg C/m2/year)
    real(dp) :: litter, resp                                           ! Fluxes per timestep (kg C/m2/timestep)
    real(dp) :: mass_start, mass_end, mass_error                       ! Mass conservation diagnostics (kg C/m2)

    !------------------------------------------------------------------!
    ! Loop counters and time
    !------------------------------------------------------------------!
    integer :: kyr, it                                                 ! Loop counters: kyr = year, it = timestep
    real(dp) :: time                                                   ! Simulation time (years)

    !------------------------------------------------------------------!
    ! Initialization
    !------------------------------------------------------------------!
    SOM = 0.0_dp                                                       ! Initial SOM (kg C/m2)
    write(*,'(a)') 'Year    SOM_C     Input    Respired'               ! Output header

    !------------------------------------------------------------------!
    ! Simulation loop: years and sub-steps
    !------------------------------------------------------------------!
    do kyr = 1, nyr
        do it = 1, nsteps
            time = (kyr - 1) + it * dt                                 ! Current simulation time (years)

            !----------------------------------------------------------!
            ! Step 1: Store SOM before update
            !----------------------------------------------------------!
            SOM_before = SOM                                           ! SOM at start of timestep (kg C/m2)

            !----------------------------------------------------------!
            ! Step 2: Compute net rate of SOM change (kg C/m2/year)
            !----------------------------------------------------------!
            dSOM = input_rate - k_decay * SOM_before                   ! Net annual rate of change

            !----------------------------------------------------------!
            ! Step 3: Update SOM with timestep-adjusted change
            !----------------------------------------------------------!
            SOM = SOM_before + dt * dSOM                               ! SOM after update (kg C/m2)

            !----------------------------------------------------------!
            ! Step 4: Compute actual timestep fluxes for output/checks
            !----------------------------------------------------------!
            litter = input_rate * dt                                   ! Litter input this step (kg C/m2/timestep)
            resp   = k_decay * SOM_before * dt                         ! Respiration loss this step (kg C/m2/timestep)

            !----------------------------------------------------------!
            ! Step 5: Mass conservation check after SOM update
            !----------------------------------------------------------!
            SOM_after = SOM                                            ! SOM now updated
            mass_start = SOM_before + litter                           ! Total C before (kg C/m2)
            mass_end   = SOM_after + resp                              ! Total C after (kg C/m2)
            mass_error = abs(mass_end - mass_start)                    ! Error magnitude

            if (mass_error > eps) then
                write(*,'(a,f6.2)') 'Mass conservation error at: ', time
                write(*,'(a,f12.5)') 'Start mass = ', mass_start
                write(*,'(a,f12.5)') 'End mass   = ', mass_end
                write(*,'(a,f12.5)') 'Difference = ', mass_error
                stop 'Mass not conserved'
            end if

        end do                                                         ! End sub-timestep loop

        !--------------------------------------------------------------!
        ! Output at end of each year
        !--------------------------------------------------------------!
        write(*,'(i5,3f12.5)') kyr, SOM, litter, resp
    end do                                                             ! End year loop

end program toymodel_v0
!======================================================================!
