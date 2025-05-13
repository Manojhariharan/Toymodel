!===============================================================================!
! Program: toymodel_v0
! Purpose: Multi-layer, 1-pool SOM model with diagnostics for mass conservation,
!          with bidirectional redistribution and historical accumulation
!          followed by degradation (time-varying decay rate; k1, k2)
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
    integer, parameter  :: nlayers = 9                                          ! Number of soil layers																						  
    !---------------------------------------------------------------------------!
    ! Model parameters
    !---------------------------------------------------------------------------!
    real(dp), parameter :: input_rate = 1.05_dp                                 ! Annual litter input (kg C/m2/year)
    real(dp), parameter :: k1_decay = 0.007_dp                                  ! Annual decay rate of SOM; accumulation phase (/year)
    real(dp), parameter :: k2_decay = 0.021_dp                                  ! Annual decay rate of SOM; degradation phase  (/year)																										   
    real(dp), parameter :: rho_SOM = 50.0_dp                                    ! SOM density (kg C/m3)

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
    real(dp) :: final_depth                                                     ! Estimated SOM depth based on rho_SOM (mm)

    !---------------------------------------------------------------------------!
    ! Loop counters
    !---------------------------------------------------------------------------!
    integer :: kyr, it, ilayer                                                  ! kyr = year, it = sub-timestep within year, ilayer = soil layer

    !---------------------------------------------------------------------------!
    ! Output file parameters
    !---------------------------------------------------------------------------!
    integer, parameter :: unit_out = 20                                         ! File unit number
    character(len=*), parameter :: outfile = 'Diagnostics.csv'                  ! Output filename

    open(unit=unit_out, file=outfile, status='replace', action='write')                                                 ! Open file for writing
    write(*,'(a)') 'Year  SOM_L1  SOM_L2  SOM_L3 SOM_L4  SOM_L6  SOM_L6 SOM_L7  SOM_L8  SOM_L9 Input  Respired  NEE'    ! Output header - screen
    write(unit_out,'(a)') 'Year,SOM_L1,SOM_L2,SOM_L3,SOM_L4,SOM_L5,SOM_L6,SOM_L7,SOM_L8,SOM_L9,input_C,respired_C,NEE'  ! Output header - file

    !---------------------------------------------------------------------------!          
    ! Initialization
    !---------------------------------------------------------------------------!
    SOM(:) = 0.0_dp                                                             ! Initial SOM per layer
    total_nee = 0.0_dp                                                          ! Initialize NEE
    
    z_interface = [0.0_dp, 45.0_dp, 91.0_dp, 166.0_dp, 289.0_dp, 493.0_dp, &
                   829.0_dp, 1383.0_dp, 2296.0_dp, 5000.0_dp]                   ! Define layer interfaces (mm)
    dz(:) = z_interface(2:) - z_interface(1:nlayers)                            ! Compute thickness of layers

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
            ! Layer-wise SOM calculation
            !-------------------------------------------------------------------!
            do ilayer = 1, nlayers                                              ! Start of SOM update layer loop

                !---------------------------------------------------------------!
                ! Step 1: Compute gross fluxes
                !---------------------------------------------------------------!
                litter(ilayer) = input_rate / real(nlayers, dp)                 ! Split annual litter input rate per layer (kg C/m2/year)

                if (kyr <= 1700) then                                           ! Select decay rate based on historical degradation (kg C/m2/year)
                    resp = k1_decay * SOM(ilayer)                               ! Decay rate for accumulation phase
                else
                    resp = k2_decay * SOM(ilayer)                               ! Decay rate for degradation phase 
                end if

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

            end do                                                              ! End of SOM update layer loop
                

            !-------------------------------------------------------------------!
            ! Step 5: SOM redistribution (bidirectional, conservative)
            !-------------------------------------------------------------------!
            do ilayer = 1, nlayers - 1                                          ! Start of layer redistribution loop

                    SOM_want = rho_SOM * dz(ilayer) / 1000.0_dp                 ! Target SOM mass (kg C/m2) based on density (kg/m3) and layer thickness (mm to m)
                    SOM_delta = min(SOM_want - SOM(ilayer), SOM(ilayer + 1))    ! Amount to transfer: positive if deficit in current layer, limited by donor pool

                    SOM(ilayer) = SOM(ilayer) + SOM_delta                       ! Receive SOM if deficit, give SOM if excess (mass-conserving adjustment)
                    SOM(ilayer + 1) = SOM(ilayer + 1) - SOM_delta               ! Mirror adjustment to adjacent (deeper) layer

            end do                                                              ! End of layer redistribution loop

            !-------------------------------------------------------------------!
            ! Compute NEE: net exchange with atmosphere 
            !-------------------------------------------------------------------!
            nee = resp_total - input_rate                                       ! Net flux with atmosphere (kg C/m2)
            
            ! Accumulation of NEE will be handled post-simulation
            ! total_nee = total_nee + dt * nee                                  ! Accumulate NEE across time (kg C/m2); 

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
        write(*,'(i5,12f12.5)') kyr, SOM(1), SOM(2), SOM(3), SOM(4), SOM(5), SOM(6), &
                SOM(7), SOM(8), SOM(9), input_rate, resp_total, nee

        write(unit_out,'(i0,",",*(f12.5,:,","))') kyr, SOM(1), SOM(2), SOM(3), SOM(4), &
                SOM(5), SOM(6), SOM(7), SOM(8), SOM(9), input_rate, resp_total, nee

    end do                                                                      ! End of year loop
    
    !---------------------------------------------------------------------------!
    ! Final diagnostics
    !---------------------------------------------------------------------------!          
    final_depth = sum(SOM(:)) / rho_SOM * 1000.0_dp                             ! Effective SOM column depth (mm)

    write(*,*)
    write(*,'(a,f12.5,a)') ' Layer 1 depth   : ', SOM(1)/rho_SOM * 1000.0_dp, ' mm'    ! Layer 1 depth (mm)
    write(*,'(a,f12.5,a)') ' Layer 2 depth   : ', SOM(2)/rho_SOM * 1000.0_dp, ' mm'    ! Layer 2 depth (mm)
    write(*,'(a,f12.5,a)') ' Layer 3 depth   : ', SOM(3)/rho_SOM * 1000.0_dp, ' mm'    ! Layer 3 depth (mm)
    write(*,'(a,f12.5,a)') ' Layer 4 depth   : ', SOM(4)/rho_SOM * 1000.0_dp, ' mm'    ! Layer 4 depth (mm)
    write(*,'(a,f12.5,a)') ' Layer 5 depth   : ', SOM(5)/rho_SOM * 1000.0_dp, ' mm'    ! Layer 5 depth (mm)
    write(*,'(a,f12.5,a)') ' Layer 6 depth   : ', SOM(6)/rho_SOM * 1000.0_dp, ' mm'    ! Layer 6 depth (mm)
    write(*,'(a,f12.5,a)') ' Layer 7 depth   : ', SOM(7)/rho_SOM * 1000.0_dp, ' mm'    ! Layer 7 depth (mm)
    write(*,'(a,f12.5,a)') ' Layer 8 depth   : ', SOM(8)/rho_SOM * 1000.0_dp, ' mm'    ! Layer 8 depth (mm)
    write(*,'(a,f12.5,a)') ' Layer 9 depth   : ', SOM(9)/rho_SOM * 1000.0_dp, ' mm'    ! Layer 9 depth (mm)
    write(*,*)
    write(*,'(a,f12.5,a)') ' Total SOM depth : ', final_depth, ' mm'                   ! SOM column depth (mm)
    write(*,'(a,f12.5,a)') ' Total SOM       : ', sum(SOM(:)),  ' kg C/m2'             ! Accumulated SOM (kg C/m2)
    
    ! NEE diagnostics will be handled post-simulation
    ! write(*,'(a,f12.5)') ' Cumulative NEE (kg C/m2): ', total_nee                    ! Total atmosphere exchange (kg C/m2)
    ! write(*,'(a,f12.5)') ' Total (SOM + NEE)       : ', sum(SOM(:)) + total_nee      ! Residual (should be ~0)
    write(*,*)

    close(unit_out)                                                             ! Close file after writing
    write(*,*) 'Simulation completed. Results saved to ', outfile

end program toymodel_v0
!===============================================================================!
