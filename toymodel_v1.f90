!===============================================================================!
! Program: toymodel_v1
! Purpose: Multi-layer, 2-pool SOM model with mass conservation diagnostics,
!          CENTURY-style decomposition (fast and slow pools), and
!          redistribution to maintain target SOM density.
!===============================================================================!
program toymodel_v1
    implicit none

    !---------------------------------------------------------------------------!
    ! Precision definition
    !---------------------------------------------------------------------------!
    integer, parameter :: dp = selected_real_kind(15, 307)                      ! Double precision for high-precision floats

    !---------------------------------------------------------------------------!
    ! Simulation control parameters
    !---------------------------------------------------------------------------!
    integer, parameter  :: nyr = 6000                                           ! Total simulation years
    real(dp), parameter :: dt = 1.0_dp                                          ! Timestep (yearly)
    integer, parameter  :: nsteps = nint(1.0_dp / dt)                           ! Number of sub-timesteps (1 for yearly timestep)
    real(dp), parameter :: eps = 1.0d-8                                         ! Tolerance for mass conservation (kg C/m2)
    integer, parameter  :: nlayers = 9                                          ! Number of soil layers

    !---------------------------------------------------------------------------!
    ! Model parameters
    !---------------------------------------------------------------------------!
    real(dp), parameter :: input_rate = 1.05_dp                                 ! Annual carbon input from litterfall (kg C/m2/year)
    real(dp), parameter :: rho_SOM = 50.0_dp                                    ! Target SOM density (kg C/m3)
    real(dp), parameter :: k_fast = 0.007_dp                                    ! Decay rate for fast SOM pool (/year)
    real(dp), parameter :: k_slow = 0.007_dp                                    ! Decay rate for slow SOM pool (/year)
    real(dp), parameter :: fCO2_fast = 1.0_dp                                   ! Fraction of fast decomposition lost as respiration
    real(dp), parameter :: fCO2_slow = 1.0_dp                                   ! Fraction of slow decomposition lost as respiration
    real(dp), parameter :: EM = 1.0_dp                                          ! Environmental modifier

    !---------------------------------------------------------------------------!
    ! Soil layer depth and thickness (mm)
    !---------------------------------------------------------------------------!
    real(dp), dimension(nlayers+1) :: z_interface = [0.0_dp, 45.0_dp, 91.0_dp, &
        166.0_dp, 289.0_dp, 493.0_dp, 829.0_dp, 1383.0_dp, 2296.0_dp, 9999.0_dp]! Layer interfaces
    real(dp), dimension(nlayers)   :: dz                                        ! Layer thicknesses (computed from interfaces)

    !---------------------------------------------------------------------------!
    ! State variables and fluxes
    !---------------------------------------------------------------------------!
    real(dp), dimension(nlayers) :: SOM_fast, SOM_slow                          ! Fast and slow SOM pools (kg C/m2)
    real(dp), dimension(nlayers) :: litter, decomp_fast, decomp_slow            ! Inputs and decomposition fluxes (kg C/m2/year)
    real(dp), dimension(nlayers) :: resp_fast, resp_slow                        ! Respiration losses from fast/slow (kg C/m2/year)
    real(dp), dimension(nlayers) :: dSOM_fast, dSOM_slow                        ! Net SOM changes per timestep
    real(dp) :: mass_start, mass_end, mass_error                                ! Mass conservation diagnostics
    real(dp) :: resp_total, nee                                                 ! Total respiration and net ecosystem exchange
    real(dp) :: final_depth                                                     ! Effective depth of SOM column (mm)
    real(dp) :: SOM_want, SOM_delta                                             ! Redistribution diagnostics

    !---------------------------------------------------------------------------!
    ! Loop counters
    !---------------------------------------------------------------------------!
    integer :: kyr, it, ilayer                                                  ! Year, timestep, and layer indices

    !---------------------------------------------------------------------------!
    ! Output file setup
    !---------------------------------------------------------------------------!
    integer, parameter :: unit_out = 20                                         ! File unit for diagnostics
    character(len=*), parameter :: outfile = 'Diagnostics.csv'                  ! Output filename

    open(unit=unit_out, file=outfile, status='replace', action='write')         ! Open output file for writing
    
    write(unit_out,'(a)') 'Year,' // 'SOMf_L1,SOMf_L2,SOMf_L3,SOMf_L4,' // &    
        'SOMf_L5,SOMf_L6,SOMf_L7,SOMf_L8,SOMf_L9,SOMs_L1,SOMs_L2,SOMs_L3,' // &
        'SOMs_L4,SOMs_L5,SOMs_L6,SOMs_L7,SOMs_L8,SOMs_L9,Input,Respired,' // &
        'NEE,TotalDepth'                                                        ! Output header - file
    
    write(*,'(a)') 'Year    SOM_fast_L1 SOM_fast_L2 ... SOM_slow_L9 Input ' // &      
        'Respired NEE TotalDepth'                                               ! Output header - screen

    !---------------------------------------------------------------------------!
    ! Initialization
    !---------------------------------------------------------------------------!
    dz = z_interface(2:) - z_interface(1:nlayers)                               ! Compute layer thicknesses (mm)
    SOM_fast(:) = 0.0_dp                                                        ! Initialize fast pool (kg C/m2)
    SOM_slow(:) = 0.0_dp                                                        ! Initialize slow pool (kg C/m2)

    !---------------------------------------------------------------------------!
    ! Simulation loop
    !---------------------------------------------------------------------------!
    do kyr = 1, nyr
        do it = 1, nsteps

            !-------------------------------------------------------------------!
            ! Mass at start: SOM total plus expected input
            !-------------------------------------------------------------------!
            mass_start = sum(SOM_fast(:)) + sum(SOM_slow(:)) + dt * input_rate  ! Total expected C before update  (kg C/m2)
            
            !-------------------------------------------------------------------!
            ! Reset flux accumulation
            !-------------------------------------------------------------------!				
            resp_total = 0.0_dp                                                 ! Reset total respiration

            !-------------------------------------------------------------------!
            ! Layer wise SOM calculation
            !-------------------------------------------------------------------!
            do ilayer = 1, nlayers

                !---------------------------------------------------------------!
                ! Step 1: Compute decomposition fluxes
                !---------------------------------------------------------------!
                litter(ilayer)      = input_rate / real(nlayers, dp)            ! Evenly distribute input across layers  (kg C/m2/year)

                decomp_fast(ilayer) = EM * k_fast * SOM_fast(ilayer)            ! Fast pool decomposition flux (kg C/m2/year)
                decomp_slow(ilayer) = EM * k_slow * SOM_slow(ilayer)            ! Slow pool decomposition flux (kg C/m2/year)
                
                resp_fast(ilayer)   = fCO2_fast * decomp_fast(ilayer)           ! Fast pool respiration (kg C/m2)
                resp_slow(ilayer)   = fCO2_slow * decomp_slow(ilayer)           ! Slow pool respiration (kg C/m2)

                !---------------------------------------------------------------!
                ! Step 2: Compute net rate of SOM change (Net annual rate of change per layer; kg C/m2/year)
                !---------------------------------------------------------------!
                dSOM_fast(ilayer) = litter(ilayer) + (1.0_dp - fCO2_slow) * &
                    decomp_slow(ilayer) - decomp_fast(ilayer)                   ! Net change in fast SOM pool (kg C/m2/year)
                
                dSOM_slow(ilayer) = (1.0_dp - fCO2_fast) * decomp_fast(ilayer) - &
                    decomp_slow(ilayer)                                         ! Net change in slow SOM pool (kg C/m2/year)

                !---------------------------------------------------------------!
                ! Step 3: Update SOM state (Updated SOM pools; kg C/m2)
                !---------------------------------------------------------------!
                SOM_fast(ilayer) = SOM_fast(ilayer) + dt * dSOM_fast(ilayer)    ! Apply timestep-adjusted net change to fast SOM pool
                SOM_slow(ilayer) = SOM_slow(ilayer) + dt * dSOM_slow(ilayer)    ! Apply timestep-adjusted net change to slow SOM pool

                !---------------------------------------------------------------!
                ! Step 4: Accumulate total respiration				
                !---------------------------------------------------------------!
                resp_total = resp_total + resp_fast(ilayer) + resp_slow(ilayer) ! Accumulate respiration (kg C/m2)

                end do

            !-------------------------------------------------------------------!
            ! SOM redistribution ((bidirectional, mass-conserving)
            !-------------------------------------------------------------------!
            do ilayer = 1, nlayers - 1
                
                SOM_want  = rho_SOM * dz(ilayer) / 1000.0_dp                    ! Target SOM mass (kg C/m2) based on density (kg/m3) and layer thickness (mm to m)
                SOM_delta = min(SOM_want - (SOM_fast(ilayer) + SOM_slow(ilayer)), &
                    SOM_fast(ilayer+1) + SOM_slow(ilayer+1))                    ! SOM to redistribute: positive if deficit in current layer, limited by donor pool
                
                SOM_slow(ilayer)   = SOM_slow(ilayer) + SOM_delta               ! Receive SOM in slow pool if deficit, give SOM if excess (mass-conserving adjustment)
                SOM_slow(ilayer+1) = SOM_slow(ilayer+1) - SOM_delta             ! Mirror adjustment to adjacent (deeper) layer
            
            end do

            !-------------------------------------------------------------------!
            ! Compute NEE: net exchange with atmosphere 
            !-------------------------------------------------------------------!
            nee = resp_total - input_rate                                       ! Net flux with atmosphere (kg C/m2)

            !-------------------------------------------------------------------!
            ! Mass conservation diagnostics
            !-------------------------------------------------------------------!
            mass_end = sum(SOM_fast(:)) + sum(SOM_slow(:)) + dt * resp_total    ! Total C after update (kg C/m2)
            mass_error = abs(mass_end - mass_start)                             ! Error magnitude		
            
            if (mass_error > eps) then
                write(*,'(a,i5)') 'Mass conservation error at year: ', kyr
                write(*,'(a,f12.5)') 'Start mass = ', mass_start
                write(*,'(a,f12.5)') 'End mass   = ', mass_end
                write(*,'(a,f12.5)') 'Difference = ', mass_error
                stop 'Mass not conserved'
            
            end if

        end do                                                                  ! End of sub-timestep loop

        !-----------------------------------------------------------------------!
        ! Annual diagnostics and output
        !-----------------------------------------------------------------------!
        final_depth = sum(SOM_fast(:)+SOM_slow(:)) / rho_SOM * 1000.0_dp        ! Total SOM depth based on sum of pools (mm)

        write(*,'(i5,27f12.5)') kyr, SOM_fast(:), SOM_slow(:), input_rate, &
            resp_total, nee, final_depth
        
         write(unit_out,'(i0,",",27(f12.5,","),f12.5)') kyr, SOM_fast(:), &
            SOM_slow(:), input_rate, resp_total, nee, final_depth

    end do                                                                      ! End of year loop

    !---------------------------------------------------------------------------!
    ! Final diagnostics
    !---------------------------------------------------------------------------!
    write(*,*)
    do ilayer = 1, nlayers                                                      ! Depth of SOM in layer i (mm)
        write(*,'(a,i1,a,f12.5,a)') ' Layer ', ilayer, ' depth   : ', &
            (SOM_fast(ilayer)+SOM_slow(ilayer))/rho_SOM * 1000.0_dp, ' mm' 
    
    end do

    write(*,*)
    write(*,'(a,f12.5,a)') ' Total SOM depth : ', final_depth, ' mm'            ! SOM column depth (mm)
    write(*,'(a,f12.5,a)') ' Total SOM       : ', sum(SOM_fast(:) + SOM_slow(:)), &
        ' kg C/m2'                                                              ! Total SOM stock across all layers (kg C/m2)
    
    write(*,*)

    close(unit_out)                                                             ! Close output file
    write(*,*) 'Simulation complete. Results saved to ', outfile

end program toymodel_v1
!===============================================================================!
