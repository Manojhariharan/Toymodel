!======================================================================!
! Program: toymodel_v0
! Purpose: Minimal 1-layer, 1-pool SOM model with clear structure and
!          mass conservation check (1 input, 1 decay, 1 pool)
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
    integer, parameter  :: nyr = 6000                                  ! Total simulation duration (years)
    real(dp), parameter :: dt = 1.0_dp                                 ! Timestep (fraction of year; e.g., 0.25 = quarterly)
    integer, parameter  :: nsteps = nint(1.0_dp / dt)                  ! Number of sub-steps per year
    real(dp), parameter :: eps = 1.0d-8                                ! Mass conservation tolerance (kg C/m2)

    !------------------------------------------------------------------!
    ! Model parameters
    !------------------------------------------------------------------!
    real(dp), parameter :: input_rate = 1.05_dp                        ! Annual litter input (kg C/m2/year)
    real(dp), parameter :: k_decay = 0.007_dp                          ! Annual decay rate of SOM (/year)

    !------------------------------------------------------------------!
    ! State and flux variables
    !------------------------------------------------------------------!
    real(dp) :: SOM                                                    ! Soil organic matter pool (kg C/m2)
    real(dp) :: litter, resp, dSOM                                     ! Fluxes and change (kg C/m2/year)
    real(dp) :: mass_start, mass_end, mass_error                       ! Mass conservation diagnostics (kg C/m2)

    !------------------------------------------------------------------!
    ! Loop counters
    !------------------------------------------------------------------!
    integer :: kyr, it                                                 ! kyr = year, it = sub-timestep within year

    !------------------------------------------------------------------!
    ! Output file parameters
    !------------------------------------------------------------------!
    integer, parameter :: unit_out = 20                                ! File unit number
    character(len=*), parameter :: outfile = 'Diagnostics.csv'         ! Output filename

    open(unit=unit_out, file=outfile, status='replace', action='write')! Open file for writing
    write(*,'(a)') 'Year    SOM_C     Input    Respired'               ! Output header - Screen
    write(unit_out,'(a)') 'Year,SOM_C,input_C,respired_C'              ! Output header - file

    !------------------------------------------------------------------!
    ! Initialization
    !------------------------------------------------------------------!
    SOM = 0.0_dp                                                       ! Initial SOM (kg C/m2)

    !------------------------------------------------------------------!
    ! Simulation loop
    !------------------------------------------------------------------!
    do kyr = 1, nyr
        do it = 1, nsteps

            !----------------------------------------------------------!
            ! Mass at start: SOM plus expected input
            !----------------------------------------------------------!
            mass_start = SOM + dt * input_rate                         ! Estimated total C before update (kg C/m2)

            !----------------------------------------------------------!
            ! Step 1: Compute gross fluxes
            !----------------------------------------------------------!
            litter = input_rate                                        ! Annual litter input rate (kg C/m2/year)
            resp   = k_decay * SOM                                     ! Annual respiration rate (kg C/m2/year)

            !----------------------------------------------------------!
            ! Step 2: Compute net rate of SOM change
            !----------------------------------------------------------!
            dSOM = litter - resp                                       ! Net annual rate of change  (kg C/m2/year)

            !----------------------------------------------------------!
            ! Step 3: Update SOM with timestep-adjusted change				
            !----------------------------------------------------------!
            SOM = SOM + dt * dSOM                                      ! SOM after update (kg C/m2)

            !----------------------------------------------------------!
            ! Mass at end: SOM after update plus CO2 loss
            !----------------------------------------------------------!
            mass_end = SOM + dt * resp                                 ! Total C after (kg C/m2)
            mass_error = abs(mass_end - mass_start)                    ! Error magnitude

            if (mass_error > eps) then
                write(*,'(a,i5)') 'Mass conservation error at year: ', kyr, ', step: ', it
                write(*,'(a,f12.5)') 'Start mass = ', mass_start
                write(*,'(a,f12.5)') 'End mass   = ', mass_end
                write(*,'(a,f12.5)') 'Difference = ', mass_error
                stop 'Mass not conserved'
            end if

        end do                                                         ! End sub-timestep loop

        !--------------------------------------------------------------!
        ! Output at end of each year
        !--------------------------------------------------------------!
        write(*,'(i5,3f12.5)') kyr, SOM, input_rate, k_decay * SOM
        write(unit_out,'(i0,",",f12.5,",",f12.5,",",f12.5)') kyr, SOM, input_rate, k_decay * SOM
    end do                                                             ! End year loop
    
        close(unit_out)                                                ! Close file after writing
        write(*,*) 'Simulation completed. Results saved to ', outfile

end program toymodel_v0
!======================================================================!
