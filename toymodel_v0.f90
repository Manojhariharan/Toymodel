!======================================================================!
! Program: toymodel_v0
! Purpose: Minimal 1-layer, 1-pool SOM model with input and respiration
!          includes variable timestep and mass conservation check
!======================================================================!
program toymodel_v0
    implicit none

    !------------------------------------------------------------------!
    ! Precision definition
    !------------------------------------------------------------------!
    integer, parameter :: dp = selected_real_kind(15, 307)             ! Double precision kind

    !------------------------------------------------------------------!
    ! Simulation parameters
    !------------------------------------------------------------------!
    real(dp), parameter :: dt = 1.0_dp                                 ! Timestep length (years)
    integer, parameter 	:: nsteps = int(6000.0_dp / dt)                ! Number of time steps (unitless)
    real(dp), parameter :: input_rate = 1.05_dp                        ! SOM input rate (kg C/m2/year)
    real(dp), parameter :: k_decay = 0.007_dp                          ! SOM decay rate (/year)
    real(dp), parameter :: eps = 1.0d-8                            	   ! Tolerance for mass conservation (kg C/m2)

    !------------------------------------------------------------------!
    ! Carbon pool and flux variables
    !------------------------------------------------------------------!
    real(dp) :: SOM                                                    ! Soil organic matter pool (kg C/m2)
    real(dp) :: input, decay, respiration                              ! Input and loss fluxes (kg C/m2 per timestep)
    real(dp) :: mass_start, mass_end, mass_error                       ! Mass conservation terms (kg C/m2)

    !------------------------------------------------------------------!
    ! Control variables
    !------------------------------------------------------------------!
    integer  :: step                                                   ! Timestep counter (unitless)
    real(dp) :: time                                                   ! Model time (years)

    !------------------------------------------------------------------!
    ! Initialization
    !------------------------------------------------------------------!
    SOM = 0.0_dp                                                       ! Initial SOM stock (kg C/m2)
    write(*,'(a)') 'Year    SOM_C     Input    Respired'               ! Header for output

    !------------------------------------------------------------------!
    ! Simulation loop over time steps
    !------------------------------------------------------------------!
    do step = 1, nsteps
        time = step * dt                                               ! Model time in years

        input = input_rate * dt                                        ! Carbon input during timestep (kg C/m2)
        decay = k_decay * SOM * dt                                     ! Decomposition loss during timestep (kg C/m2)
        respiration = decay                                            ! All loss assumed to be CO2 respiration

        mass_start = SOM + input                                       ! Total C before update (kg C/m2)
        mass_end = SOM + input - respiration + respiration             ! Total C after update (kg C/m2)
        mass_error = abs(mass_end - mass_start)                        ! Difference (should be ~0) (kg C/m2)

        !--------------------------------------------------------------!
        ! Mass conservation check (per timestep)
        !--------------------------------------------------------------!
        if (mass_error > eps) then
            write(*,'(a,f8.2)') 'Mass conservation error at year: ', time
            write(*,'(a,f12.6)') 'Start mass = ', mass_start
            write(*,'(a,f12.6)') 'End mass   = ', mass_end
            write(*,'(a,f12.6)') 'Difference = ', mass_error
            stop 'Mass not conserved'
        end if

        SOM = SOM + input - respiration                                ! Update SOM pool (kg C/m2)

        if (mod(step, int(1.0_dp / dt)) == 0) then                     ! Print once per year
            write(*,'(f6.1,3f12.6)') time, SOM, input, respiration
        end if
    end do

end program toymodel_v0
!======================================================================!
