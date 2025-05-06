!======================================================================!
! Program: toymodel_v0
! Purpose: Minimal 1-layer, 1-pool SOM model with input and respiration
!          includes per-timestep mass conservation check using double precision
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
    integer, parameter :: nyears = 6000                                ! Total simulation years (years)
    real(dp), parameter :: input_rate = 1.05_dp                        ! SOM input (kg C/m2/year)
    real(dp), parameter :: k_decay = 0.007_dp                          ! Decay rate (/year)
    real(dp), parameter :: eps = 1.0e-8_dp                             ! Mass conservation tolerance

    !------------------------------------------------------------------!
    ! Carbon pool and flux variables
    !------------------------------------------------------------------!
    real(dp) :: SOM                                                    ! SOM pool (kg C/m2)
    real(dp) :: input, decay, respiration                              ! Annual input, decay, and respiration fluxes
    real(dp) :: mass_start, mass_end, mass_error                       ! Mass conservation diagnostics
    real(dp) :: total_input, total_respired                            ! Cumulative totals (kg C/m2)

    !------------------------------------------------------------------!
    ! Control variables
    !------------------------------------------------------------------!
    integer :: t                                                       ! Current year in the simulation loop

    !------------------------------------------------------------------!
    ! Initialization
    !------------------------------------------------------------------!
    SOM = 0.0_dp                                                       ! Initial SOM stock (kg C/m2)
    total_input = 0.0_dp                                               ! Initialize input sum
    total_respired = 0.0_dp                                            ! Initialize respiration sum
    write(*,'(a)') 'Year    SOM_C     Input    Respired'

    !------------------------------------------------------------------!
    ! Simulation loop over years
    !------------------------------------------------------------------!
    do t = 1, nyears
        input = input_rate                                             ! Annual litter input (kg C/m2)
        decay = k_decay * SOM                                          ! First-order decay loss (kg C/m2)
        respiration = decay                                            ! All decay results in respiration

        mass_start = SOM + input                                       ! Mass before update (kg C/m2)

        SOM = SOM + input - respiration                                ! Update SOM stock (kg C/m2)

        mass_end = SOM + respiration                                   ! Mass after update (kg C/m2)
        mass_error = abs(mass_end - mass_start)                        ! Difference from mass conservation (kg C/m2)

        !--------------------------------------------------------------!
        ! Mass conservation check (per timestep)
        !--------------------------------------------------------------!
        if (mass_error > eps) then
            write(*,'(a,i5)') 'Mass conservation error at year: ', t
            write(*,'(a,f12.4)') 'Start mass = ', mass_start
            write(*,'(a,f12.4)') 'End mass   = ', mass_end
            write(*,'(a,f12.4)') 'Difference = ', mass_error
            stop 'Mass not conserved'
        end if

        total_input = total_input + input                              ! Total accumulated input (kg C/m2)
        total_respired = total_respired + respiration                  ! Total accumulated respiration (kg C/m2)

        write(*,'(i4,3f12.4)') t, SOM, input, respiration
    end do

    !------------------------------------------------------------------!
    ! Final mass conservation summary
    !------------------------------------------------------------------!
    write(*,*)
    write(*,'(a,f12.4)') 'Total C input:      ', total_input
    write(*,'(a,f12.4)') 'Total C respired:   ', total_respired
    write(*,'(a,f12.4)') 'Final SOM stock:    ', SOM
    write(*,'(a,f12.4)') 'Residual: ', total_input - total_respired - SOM

end program toymodel_v0
!======================================================================!
