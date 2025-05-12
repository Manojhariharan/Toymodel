!===============================================================================!
! Program: toymodel_v0
! Purpose: 3-layer, 1-pool SOM model with diagnostics for NEE and mass 
!          conservation, with layer redistribution to maintain fixed SOM 
!          mass-density targets
!===============================================================================!
program toymodel_v0
    implicit none

    !---------------------------------------------------------------------------!
    ! Precision definition
    !---------------------------------------------------------------------------!
    integer, parameter :: dp = selected_real_kind(15, 307)                      ! Double precision

    !---------------------------------------------------------------------------!
    ! Simulation control parameters
    !---------------------------------------------------------------------------!
    integer, parameter  :: nyr = 6000                                           ! Total simulation duration (years)
    real(dp), parameter :: dt = 1.00_dp                                         ! Timestep (fraction of year; e.g., 0.25 = quarterly)
    integer, parameter  :: nsteps = nint(1.0_dp / dt)                           ! Number of sub-steps per year
    real(dp), parameter :: eps = 1.0d-8                                         ! Mass conservation tolerance (kg C/m2)
    integer, parameter  :: nlayers = 3                                          ! Number of soil layers																						  
    !---------------------------------------------------------------------------!
    ! Model parameters
    !---------------------------------------------------------------------------!
    real(dp), parameter :: input_rate = 1.05_dp                                 ! Annual litter input (kg C/m2/year)
    real(dp), parameter :: k_decay = 0.007_dp                                   ! Annual decay rate of SOM (/year)
    real(dp), parameter :: rho_SOM = 0.03_dp                                    ! Target SOM density (kg C/m3)

    !---------------------------------------------------------------------------!
    ! Soil layer depth and thickness
    !---------------------------------------------------------------------------!
    real(dp), dimension(nlayers+1) :: z_interface                               ! Layer interface depths (mm)
    real(dp), dimension(nlayers)   :: dz                                        ! Layer thickness (mm)

    !---------------------------------------------------------------------------!
    ! State and flux variables (arrays for layers)
    !---------------------------------------------------------------------------!
    real(dp) :: SOM(nlayers)                                                    ! Soil organic matter per layer (kg C/m2)
    real(dp) :: litter(nlayers), dSOM(nlayers)                                  ! Fluxes and change (kg C/m2)
    real(dp) :: resp_total, resp                                                ! Total and layer respiration (kg C/m2)																												  
    real(dp) :: mass_start, mass_end, mass_error                                ! Mass conservation diagnostics (kg C/m2)
    real(dp) :: nee, total_nee                                                  ! Net ecosystem exchange (kg C/m2/year)																										  
    real(dp) :: SOM_want, SOM_delta                                             ! Redistribution diagnostics (kg C/m2)																												  
    !---------------------------------------------------------------------------!
    ! Loop counters
    !---------------------------------------------------------------------------!
    integer :: kyr, it, ilayer                                                  ! kyr = year, it = sub-timestep within year, ilayer = soil layer

    !---------------------------------------------------------------------------!
    ! Output file parameters
    !---------------------------------------------------------------------------!
    integer, parameter :: unit_out = 20                                         ! File unit number
    character(len=*), parameter :: outfile = 'Diagnostics.csv'                  ! Output filename

    open(unit=unit_out, file=outfile, status='replace', action='write')                    ! Open file for writing
    write(*,'(a)') 'Year    SOM_C_L1  SOM_C_L2   SOM_C_L3   Input    Respired   NEE'       ! Output header - screen
    write(unit_out,'(a)') 'Year,SOM_C_L1,SOM_C_L2,SOM_C_L3,input_C,respired_C,NEE'         ! Output header - file

    !---------------------------------------------------------------------------!
    ! Initialization
    !---------------------------------------------------------------------------!
    SOM(:) = 0.0_dp                                                             ! Initial SOM per layer
    total_nee = 0.0_dp                                                          ! Initialize NEE
    
    z_interface = [0.0_dp, 45.0_dp, 91.0_dp, -999.0_dp]                         ! Define layer interfaces (mm)
    dz(1:nlayers-1) = z_interface(2:nlayers) - z_interface(1:nlayers-1)         ! Compute thickness of first two layers
    dz(nlayers) = -1.0_dp                                                       ! Placeholder for bottomless layer

    !---------------------------------------------------------------------------!
    ! Simulation loop
    !---------------------------------------------------------------------------!
    do kyr = 1, nyr                                                             ! Start of year loop
        do it = 1, nsteps                                                       ! Start of sub-timestep loop

            !-------------------------------------------------------------------!
            ! Mass at start: SOM total plus expected input
            !-------------------------------------------------------------------!
            mass_start = sum(SOM(:)) + dt * input_rate                          ! Estimated total C before update (kg C/m2)

            !-------------------------------------------------------------------!
            ! Reset flux accumulation
            !-------------------------------------------------------------------!				
            resp_total = 0.0_dp                                                 ! Reset total respiration
            
            !-------------------------------------------------------------------!
            ! Layer-wise calculation
            !-------------------------------------------------------------------!
            do ilayer = 1, nlayers                                              ! Start of layer loop

                !---------------------------------------------------------------!
                ! Step 1: Compute gross fluxes
                !---------------------------------------------------------------!
                litter(ilayer) = input_rate / real(nlayers, dp)                 ! Split annual litter input rate per layer (kg C/m2/year)
                resp = k_decay * SOM(ilayer)                                    ! Annual respiration rate per layer (kg C/m2/year)

                !---------------------------------------------------------------!
                ! Step 2: Compute net rate of SOM change
                !---------------------------------------------------------------!
                dSOM(ilayer) = litter(ilayer) - resp                            ! Net annual rate of change per layer (kg C/m2/year)

                !---------------------------------------------------------------!
                ! Step 3: Update SOM with timestep-adjusted change				
                !---------------------------------------------------------------!
                SOM(ilayer) = SOM(ilayer) + dt * dSOM(ilayer)                   ! SOM after update (kg C/m2)
            
                !---------------------------------------------------------------!
                ! Step 4: Accumulate total respiration				
                !---------------------------------------------------------------!
                resp_total = resp_total + resp                                  ! Accumulate respiration (kg C/m2)

                !---------------------------------------------------------------!
                ! Step 5: SOM redistribution (skip bottom layer)
                !---------------------------------------------------------------!
                if (ilayer /= nlayers .and. dz(ilayer) > 0.0_dp) then
                    SOM_want = rho_SOM * dz(ilayer) / 1000.0_dp                 ! Target SOM (convert mm to m)
                    SOM_delta = SOM(ilayer) - SOM_want                          ! Move excess SOM downward
                    SOM(ilayer) = SOM(ilayer) - SOM_delta                       ! Subtract excess SOM from current layer
                    SOM(ilayer+1) = SOM(ilayer+1) + SOM_delta                   ! Add excess to next (deeper) layer
                end if

            end do                                                              ! End of layer loop

            !-------------------------------------------------------------------!
            ! Enforce non-negative SOM after redistribution
            !-------------------------------------------------------------------!
            do ilayer = 1, nlayers
                if (SOM(ilayer) < 0.0_dp) SOM(ilayer) = 0.0_dp
            end do

            !-------------------------------------------------------------------!
            ! Compute NEE: net exchange with atmosphere 
            !-------------------------------------------------------------------!
            nee = resp_total - input_rate                                       ! Net flux with atmosphere (kg C/m2)
            total_nee = total_nee + dt * nee                                    ! Accumulate NEE across time (kg C/m2)

            !-------------------------------------------------------------------!
            ! Mass at end: SOM after update plus CO2 loss
            !-------------------------------------------------------------------!
            mass_end = sum(SOM(:)) + dt * resp_total                            ! Total C after (kg C/m2)
            mass_error = abs(mass_end - mass_start)                             ! Error magnitude

            if (mass_error > eps) then
                write(*,'(a,i5)') 'Mass conservation error at year: ', kyr, ', step: ', it
                write(*,'(a,f12.5)') 'Start mass = ', mass_start
                write(*,'(a,f12.5)') 'End mass   = ', mass_end
                write(*,'(a,f12.5)') 'Difference = ', mass_error
                stop 'Mass not conserved'
            end if

        end do                                                                  ! End of sub-timestep loop

        !-----------------------------------------------------------------------!																  
        ! Output at end of each year
        !-----------------------------------------------------------------------!
        write(*,'(i5,6f12.5)') kyr, SOM(1), SOM(2), SOM(3), input_rate, resp_total, nee
        write(unit_out,'(i0,",",f12.5,",",f12.5,",",f12.5,",",f12.5,",",f12.5,",",f12.5)') &
                kyr, SOM(1), SOM(2), SOM(3), input_rate, resp_total, nee
    end do                                                                      ! End of year loop
    
    !---------------------------------------------------------------------------!
    ! Final diagnostics
    !---------------------------------------------------------------------------!
    write(*,*)
    write(*,'(a,f12.5)') ' Total SOM (kg C/m2)     : ', sum(SOM(:))             ! Accumulated SOM (kg C/m2)
    write(*,'(a,f12.5)') ' Cumulative NEE (kg C/m2): ', total_nee               ! Total atmosphere exchange (kg C/m2)
    write(*,'(a,f12.5)') ' Total (SOM + NEE)       : ', sum(SOM(:)) + total_nee ! Residual (should be ~0)
    write(*,*)

    close(unit_out)                                                             ! Close file after writing
    write(*,*) 'Simulation completed. Results saved to ', outfile

end program toymodel_v0
!===============================================================================!
