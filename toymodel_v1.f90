!===============================================================================!
! Program: toymodel_v1
! Purpose: Multi-layer, 3-pool SOM model with mass conservation diagnostics,
!          CENTURY-style decomposition (fast/ slow/ passive pools),
!          and redistribution to maintain target SOM density.
!          Surface vs root litter logic based on Parton et al. (1993)
!===============================================================================!
program toymodel_v1
    implicit none

    !---------------------------------------------------------------------------!
    ! Precision definition
    !---------------------------------------------------------------------------!
    integer, parameter :: dp = selected_real_kind(15, 307)                                      ! Double precision for high-precision floats

    !---------------------------------------------------------------------------!
    ! Simulation control parameters
    !---------------------------------------------------------------------------!
    integer, parameter  :: nyr = 6000                                                           ! Total simulation years
    real(dp), parameter :: dt = 1.0_dp                                                          ! Timestep (yearly)
    integer, parameter  :: nsteps = nint(1.0_dp / dt)                                           ! Number of sub-timesteps (1 for yearly timestep)
    real(dp), parameter :: eps = 1.0d-8                                                         ! Tolerance for mass conservation (kg C/m2)
    integer, parameter  :: nlayers = 9                                                          ! Number of soil layers

    !---------------------------------------------------------------------------!
    ! Model parameters
    !---------------------------------------------------------------------------!
    real(dp), parameter :: input_rate = 1.05_dp                                                 ! Annual carbon input from litter (kg C/m2/year)
    real(dp), parameter :: rho_SOM = 50.0_dp                                                    ! SOM density (kg C/m3)
    real(dp), parameter :: k_fast = 1.75_dp                                                     ! Adjusted decay rate for fast SOM pool from 0.07 (/year)
    real(dp), parameter :: k_slow = 0.04175_dp                                                  ! Adjusted decay rate for slow SOM pool from 0.00167 (/year)
    real(dp), parameter :: k_passive = 0.0025_dp                                                ! Adjusted decay rate for passive SOM pool from 0.0001 (/year)
    real(dp), parameter :: fCO2_fast = 0.55_dp                                                  ! Fraction of fast decomposition lost as respiration
    real(dp), parameter :: fCO2_slow = 0.55_dp                                                  ! Fraction of slow decomposition lost as respiration
    real(dp), parameter :: fCO2_passive = 0.55_dp                                               ! Fraction of passive decomposition lost as respiration
    real(dp), parameter :: EM = 0.1_dp                                                          ! Environmental modifier revised from 2.5
    real(dp), parameter :: clay_frac = 0.2_dp                                                   ! Clay fraction (for CAP, CSP)
    real(dp), parameter :: sand_frac = 0.6_dp                                                   ! Sand fraction (if needed later)

    !---------------------------------------------------------------------------!
    ! Litter input decomposition parameters (for structural/metabolic)
    !---------------------------------------------------------------------------!
    real(dp), parameter :: L_N_ratio = 8.0_dp                                                   ! Lignin to Nitrogen ratio
    real(dp), parameter :: Fm = 0.99_dp - 0.018_dp * L_N_ratio                                  ! Metabolic fraction
    real(dp), parameter :: Ls = 0.25_dp                                                         ! Fraction of structural C that is lignin
    real(dp), parameter :: LC = exp(-3.0_dp * Ls)                                               ! Lignin control factor
    real(dp), parameter :: k_struct = 4.8_dp                                                    ! Structural decay rate (/year)
    real(dp), parameter :: k_metabolic = 18.5_dp                                                ! Metabolic decay rate (/year)
    real(dp), parameter :: fCO2_struct_lig = 0.3_dp                                             ! Respiration from lignin decay
    real(dp), parameter :: fCO2_struct_cel = 0.6_dp                                             ! Respiration from cellulose decay
    real(dp), parameter :: fCO2_metabolic = 0.55_dp                                             ! Respiration from metabolic decay

    !---------------------------------------------------------------------------!
    ! Century-style decomposition transfer fractions (computed each year)
    !---------------------------------------------------------------------------!
    real(dp) :: CAP                                                                             ! Fast -> passive pool (root only)
    real(dp) :: CSP                                                                             ! Slow -> passive pool
    real(dp) :: CSA                                                                             ! Slow -> fast pool (root only)

    !---------------------------------------------------------------------------!
    ! Soil layer depth and thickness (mm)
    !---------------------------------------------------------------------------!
    real(dp), dimension(nlayers+1) :: z_interface = [0.0_dp, 45.0_dp, 91.0_dp, &              
        166.0_dp, 289.0_dp, 493.0_dp, 829.0_dp, 1383.0_dp, 2296.0_dp, 9999.0_dp]                ! Layer interfaces
    real(dp), dimension(nlayers)   :: dz                                                        ! Layer thicknesses (computed from interfaces)

    !---------------------------------------------------------------------------!
    ! State variables and fluxes
    !---------------------------------------------------------------------------!
    real(dp), dimension(nlayers) :: litter                                                      ! Input fluxes (kg C/m2/year)
    real(dp), dimension(nlayers) :: litter_struct_lig, litter_struct_cel, litter_metabolic      ! Litter pools per layer (kg C/m2)
    real(dp), dimension(nlayers) :: SOM_fast, SOM_slow, SOM_passive                             ! Fast/slow/passive SOM pools (kg C/m2)
    real(dp), dimension(nlayers) :: decomp_fast, decomp_slow, decomp_passive                    ! Decomposition fluxes (kg C/m2/year)
    real(dp), dimension(nlayers) :: decomp_struct_lig, decomp_struct_cel, decomp_metabolic      ! Litter pool decomposition fluxes
    real(dp), dimension(nlayers) :: resp_fast, resp_slow, resp_passive                          ! Respiration losses from fast/slow/passive pools
    real(dp), dimension(nlayers) :: resp_struct_lig, resp_struct_cel, resp_metabolic            ! Respiration from structural/metabolic pools
    real(dp), dimension(nlayers) :: dSOM_fast, dSOM_slow, dSOM_passive                          ! Net SOM changes per timestep

    real(dp) :: mass_start, mass_end, mass_error                                                ! Mass conservation diagnostics
    real(dp) :: resp_total, nee                                                                 ! Total respiration and net ecosystem exchange
    real(dp) :: final_depth                                                                     ! Effective depth of SOM column
    real(dp) :: SOM_want, SOM_delta                                                             ! Redistribution diagnostics
    real(dp) :: SOM_current, SOM_below, move_amt                                                ! Local SOM stocks and actual redistributed amount
    real(dp) :: frac_fast_donor, frac_slow_donor, frac_passive_donor                            ! Fraction of fast/ slow/ passive SOM in donor layer
    real(dp) :: fast_to_slow, slow_to_fast, slow_to_passive, passive_to_slow                    ! Non-respired fraction of fast/slow/passive pool decomposition
    real(dp) :: fast_to_passive                                                                 ! NEW: Fast-to-passive flux in root layers
    real(dp) :: fast_gain, fast_lose                                                            ! Transfer between fast SOM pool decomposition
    real(dp) :: slow_gain, slow_lose                                                            ! Transfer between slow SOM pool decomposition
    real(dp) :: passive_gain, passive_lose                                                      ! Transfer between passive SOM pool decomposition
    real(dp) :: total_thickness                                                                 ! Total soil profile thickness (mm)

    !---------------------------------------------------------------------------!
    ! Loop counters
    !---------------------------------------------------------------------------!
    integer :: kyr, it, ilayer                                                                  ! Year, timestep, and layer indices
    !---------------------------------------------------------------------------!
    ! Inter-layer pool ratio warnings check
    !---------------------------------------------------------------------------!
    logical, parameter :: enable_ratio_check = .false.                                          ! Set to .true. to enable warnings

    !---------------------------------------------------------------------------!
    ! Output file setup
    !---------------------------------------------------------------------------!
    integer, parameter :: unit_out = 20                                                         ! File unit for diagnostics
    character(len=*), parameter :: outfile = 'Diagnostics.csv'                                  ! Diagnostics for total output

    integer, parameter :: unit_layer = 21                                                       ! Diagnostics for per-layer fast/slow pool evolution
    character(len=*), parameter :: layerfile = 'LayerwisePools.csv'

    open(unit=unit_out, file=outfile, status='replace', action='write')                         ! Open diagnostics output file for writing
    open(unit=unit_layer, file=layerfile, status='replace', action='write')                     ! Open pool output file for writing

    write(unit_out,'(a)') 'Year,' // 'SOMf_L1,SOMf_L2,SOMf_L3,SOMf_L4,' // &    
        'SOMf_L5,SOMf_L6,SOMf_L7,SOMf_L8,SOMf_L9,SOMs_L1,SOMs_L2,' // &
        'SOMs_L3,SOMs_L4,SOMs_L5,SOMs_L6,SOMs_L7,SOMs_L8,SOMs_L9,' // &
        'SOMp_L1,SOMp_L2,SOMp_L3,SOMp_L4,SOMp_L5,SOMp_L6,SOMp_L7,' // &
        'SOMp_L8,SOMp_L9,Input,Respired,NEE,TotalDepth'                                         ! Output header - Diagnostics file

    write(*,'(a)') 'Year  ' // 'SOMf_L1 SOMf_L2 SOMf_L3 SOMf_L4 SOMf_L5 ' // &
        'SOMf_L6 SOMf_L7 SOMf_L8 SOMf_L9 SOMs_L1 SOMs_L2 SOMs_L3 SOMs_L4 ' // &
        'SOMs_L5 SOMs_L6 SOMs_L7 SOMs_L8 SOMs_L9 SOMp_L1 SOMp_L2 SOMp_L3 ' // &
        'SOMp_L4 SOMp_L5 SOMp_L6 SOMp_L7 SOMp_L8 SOMp_L9 Input Respired ' // &
        'NEE TotalDepth'                                                                        ! Output header - screen

    write(unit_layer,*) 'Year,Layer,SOM_fast,SOM_slow,SOM_passive,' // &
        'SOM_total,FastSlow_ratio,SlowPassive_ratio,FastPassive_ratio'                          ! Output header - Pool file

    !---------------------------------------------------------------------------!
    ! Initialization
    !---------------------------------------------------------------------------!
    dz = z_interface(2:) - z_interface(1:nlayers)                                               ! Compute layer thicknesses (mm)
    total_thickness = sum(dz)                                                                   ! Compute total thickness of all soil layers (mm)

    SOM_fast(:) = 0.0_dp                                                                        ! Initialize fast pool (kg C/m2)
    SOM_slow(:) = 0.0_dp                                                                        ! Initialize slow pool (kg C/m2)
    SOM_passive(:) = 0.0_dp                                                                     ! Initialize passive pool (kg C/m2)
    litter_struct_lig(:) = 0.0_dp                                                               ! Initialize structural lignin litter pool (kg C/m2)
    litter_struct_cel(:) = 0.0_dp                                                               ! Initialize structural cellulose litter pool (kg C/m2)
    litter_metabolic(:)  = 0.0_dp                                                               ! Initialize metabolic litter pool (kg C/m2)

    !---------------------------------------------------------------------------!
    ! Check for depth boundary
    !---------------------------------------------------------------------------!
    if (z_interface(nlayers+1) > 9999.0_dp) then
        write(*,*) 'Warning: Bottom layer interface exceeds 9 m; check boundary condition'
    
    end if

    !---------------------------------------------------------------------------!
    ! Simulation loop
    !---------------------------------------------------------------------------!
    do kyr = 1, nyr
        do it = 1, nsteps

            !-------------------------------------------------------------------!
            ! Compute Century-style transfer coefficients (root litter logic)
            !-------------------------------------------------------------------!
            CAP = 0.003_dp + 0.032_dp * clay_frac                                               ! Fast pool -> passive pool (root litter logic)
            CSP = 0.003_dp - 0.009_dp * clay_frac                                               ! Slow pool -> passive pool (for all layers)
            CSA = 1.0_dp - CAP - 0.55_dp                                                        ! Slow pool -> fast pool (root litter logic)

            !-------------------------------------------------------------------!
            ! Reset flux accumulation
            !-------------------------------------------------------------------!
            resp_total = 0.0_dp                                                                 ! Reset total respiration (kg C/m2)

            !-------------------------------------------------------------------!
            ! Layer wise SOM calculation
            !-------------------------------------------------------------------!
            do ilayer = 1, nlayers

                litter(ilayer) = input_rate * dz(ilayer) / total_thickness                      ! Distribute litter input proportional to each layer's thickness

                !---------------------------------------------------------------!
                ! Step 0: Partition litter input into structural and metabolic pools
                !---------------------------------------------------------------!
                if (ilayer == 1) then
                    
                    ! Surface litter decomposition logic (surface microbes)
                    litter_struct_lig(ilayer) = litter_struct_lig(ilayer) + Ls * &
                        (1.0_dp - Fm) * litter(ilayer)                                          ! Structural lignin input

                    litter_struct_cel(ilayer) = litter_struct_cel(ilayer) + &
                        (1.0_dp - Ls) * (1.0_dp - Fm) * litter(ilayer)                          ! Structural cellulose input

                    litter_metabolic(ilayer)  = litter_metabolic(ilayer) + &
                        Fm * litter(ilayer)                                                     ! Metabolic input

                else
                    
                    ! Root litter decomposition logic (soil microbes)
                    litter_struct_lig(ilayer) = litter_struct_lig(ilayer) + Ls * &
                        (1.0_dp - Fm) * litter(ilayer)                                          ! Structural lignin input

                    litter_struct_cel(ilayer) = litter_struct_cel(ilayer) + &
                        (1.0_dp - Ls) * (1.0_dp - Fm) * litter(ilayer)                          ! Structural cellulose input

                    litter_metabolic(ilayer)  = litter_metabolic(ilayer) + &
                        Fm * litter(ilayer)                                                     ! Metabolic input
                end if
 
                !---------------------------------------------------------------!
                ! Step 1: Structural lignin decomposition -> slow pool
                !---------------------------------------------------------------!
                decomp_struct_lig(ilayer) = EM * k_struct * LC * litter_struct_lig(ilayer)      ! Structural lignin litter decay rate (kg C/m2/year)
                resp_struct_lig(ilayer)   = fCO2_struct_lig * decomp_struct_lig(ilayer)         ! Respiration from structural lignin decay (kg C/m2/year)
                
                SOM_slow(ilayer) = SOM_slow(ilayer) + (1.0_dp - fCO2_struct_lig) * &
                    decomp_struct_lig(ilayer)                                                   ! Transfer to slow pool (non-respired fraction)

                litter_struct_lig(ilayer) = litter_struct_lig(ilayer) - &
                    decomp_struct_lig(ilayer)                                                   ! Update structural lignin pool after decay

                !---------------------------------------------------------------!
                ! Step 2: Structural cellulose decomposition -> fast pool
                !---------------------------------------------------------------!
                decomp_struct_cel(ilayer) = EM * k_struct * LC * litter_struct_cel(ilayer)      ! Structural cellulose litter decay rate (kg C/m2/year)
                resp_struct_cel(ilayer)   = fCO2_struct_cel * decomp_struct_cel(ilayer)         ! Respiration from structural cellulose decay (kg C/m2/year)
                
                SOM_fast(ilayer) = SOM_fast(ilayer) + (1.0_dp - fCO2_struct_cel) * &
                    decomp_struct_cel(ilayer)                                                   ! Transfer to fast pool (non-respired fraction)
                
                litter_struct_cel(ilayer) = litter_struct_cel(ilayer) - &
                    decomp_struct_cel(ilayer)                                                   ! Update structural cellulose pool after decay

                !---------------------------------------------------------------!
                ! Step 3: Metabolic litter decomposition -> fast pool
                !---------------------------------------------------------------!
                decomp_metabolic(ilayer) = EM * k_metabolic * litter_metabolic(ilayer)          ! Metabolic litter decay rate (kg C/m2/year)
                resp_metabolic(ilayer)   = fCO2_metabolic * decomp_metabolic(ilayer)            ! Respiration from metabolic decay (kg C/m2/year)
                
                SOM_fast(ilayer) = SOM_fast(ilayer) + (1.0_dp - fCO2_metabolic) * &
                    decomp_metabolic(ilayer)                                                    ! Transfer to fast pool (non-respired fraction)

                litter_metabolic(ilayer) = litter_metabolic(ilayer) - decomp_metabolic(ilayer)  ! Update metabolic pool after decay

                !---------------------------------------------------------------!
                ! Step 4: Compute decomposition fluxes
                !---------------------------------------------------------------!
                decomp_fast(ilayer) = EM * k_fast * SOM_fast(ilayer)                            ! Fast pool decomposition flux (kg C/m2/year)
                decomp_slow(ilayer) = EM * k_slow * SOM_slow(ilayer)                            ! Slow pool decomposition flux (kg C/m2/year)
                decomp_passive(ilayer)= EM * k_passive * SOM_passive(ilayer)                    ! Passive pool decomposition flux (kg C/m2/year)

                resp_fast(ilayer)   = fCO2_fast * decomp_fast(ilayer)                           ! Fast pool respiration (kg C/m2)
                resp_slow(ilayer)   = fCO2_slow * decomp_slow(ilayer)                           ! Slow pool respiration (kg C/m2)
                resp_passive(ilayer)  = fCO2_passive * decomp_passive(ilayer)                   ! Passive pool respiration (kg C/m2)

                !---------------------------------------------------------------!
                ! Step 5: Define decomposition-based carbon transfers between pools
                !---------------------------------------------------------------!
                if (ilayer == 1) then
                    
                    ! Surface layer logic (original)
                    fast_to_slow = (1.0_dp - fCO2_fast) * decomp_fast(ilayer)
                    slow_to_passive = (1.0_dp - fCO2_slow) * CSP * decomp_slow(ilayer)
                    slow_to_fast = (1.0_dp - fCO2_slow) * (1.0_dp - CSP) * decomp_slow(ilayer)
                    passive_to_slow = (1.0_dp - fCO2_passive) * decomp_passive(ilayer)
                    fast_to_passive = 0.0_dp                                                    ! No CAP transfer in surface layer

                else
                    
                    ! Root layer logic (CENTURY-style)
                    fast_to_passive = (1.0_dp - fCO2_fast) * CAP * decomp_fast(ilayer)
                    fast_to_slow    = (1.0_dp - fCO2_fast) * (1.0_dp - CAP) * decomp_fast(ilayer)
                    slow_to_passive = (1.0_dp - fCO2_slow) * CSP * decomp_slow(ilayer)
                    slow_to_fast    = (1.0_dp - fCO2_slow) * CSA * decomp_slow(ilayer)
                    passive_to_slow = (1.0_dp - fCO2_passive) * decomp_passive(ilayer)

                end if

                !---------------------------------------------------------------!
                ! Step 6: Advance SOM pools based on gain and loss fluxes
                !---------------------------------------------------------------!
                fast_gain = slow_to_fast                                                       ! Fast gains from slow decay
                
                if (ilayer > 1) then
                    fast_gain = fast_gain + (1.0_dp - fCO2_struct_cel) * decomp_struct_cel(ilayer) + &
                        (1.0_dp - fCO2_metabolic)   * decomp_metabolic(ilayer)
                
                end if
                
                fast_lose = decomp_fast(ilayer)                                                ! Fast pool total outflow

                slow_gain = fast_to_slow + passive_to_slow                                     ! Slow pool gains
                slow_lose = decomp_slow(ilayer)                                                ! Slow pool total decay

                passive_gain  = slow_to_passive + fast_to_passive                              ! Passive pool gains
                passive_lose  = decomp_passive(ilayer)                                         ! Passive pool decay

                dSOM_fast(ilayer) = fast_gain - fast_lose                                      ! Net change in fast SOM
                dSOM_slow(ilayer) = slow_gain - slow_lose                                      ! Net change in slow SOM
                dSOM_passive(ilayer) = passive_gain - passive_lose                             ! Net change in passive SOM

                !---------------------------------------------------------------!
                ! Step 7: Update SOM state (Updated SOM pools; kg C/m2)
                !---------------------------------------------------------------!
                SOM_fast(ilayer) = SOM_fast(ilayer) + dt * dSOM_fast(ilayer)
                SOM_slow(ilayer) = SOM_slow(ilayer) + dt * dSOM_slow(ilayer)
                SOM_passive(ilayer) = SOM_passive(ilayer) + dt * dSOM_passive(ilayer)

                !---------------------------------------------------------------!
                ! Step 8: Accumulate total respiration
                !---------------------------------------------------------------!
                resp_total = resp_total + resp_fast(ilayer) + resp_slow(ilayer) + &
                    resp_passive(ilayer) + resp_struct_lig(ilayer) + &
                    resp_struct_cel(ilayer) + resp_metabolic(ilayer)

            end do  ! End of ilayer loop

            !-------------------------------------------------------------------!
            ! Mass at start: SOM total plus expected input
            !-------------------------------------------------------------------!
            mass_start = sum(SOM_fast(:)) + sum(SOM_slow(:)) + sum(SOM_passive(:)) + &
                sum(litter_struct_lig(:)) + sum(litter_struct_cel(:)) + &
                sum(litter_metabolic(:))                                                        ! Total expected C before update (kg C/m2)

            !-------------------------------------------------------------------!
            ! SOM redistribution (bidirectional, 3-pool, mass-conserving)
            !-------------------------------------------------------------------!
            do ilayer = 1, nlayers - 1

                SOM_want    = rho_SOM * dz(ilayer) / 1000.0_dp                                  ! Target SOM mass in layer (kg C/m2)
                SOM_current = SOM_fast(ilayer) + SOM_slow(ilayer) + SOM_passive(ilayer)         ! Current SOM in layer (kg C/m2)
                SOM_below   = SOM_fast(ilayer+1) + SOM_slow(ilayer+1) + SOM_passive(ilayer+1)   ! SOM in layer below (kg C/m2)

                SOM_delta = SOM_want - SOM_current                                              ! Net deficit or excess in current layer

                if (SOM_delta > 0.0_dp .and. SOM_below > 0.0_dp) then
                    !-----------------------------------------------------------!
                    ! Case 1: Upward SOM movement (deficit in current layer)
                    !-----------------------------------------------------------!
                    frac_fast_donor    = SOM_fast(ilayer+1) / SOM_below                         ! Fast fraction in donor layer
                    frac_slow_donor    = SOM_slow(ilayer+1) / SOM_below                         ! Slow fraction in donor layer
                    frac_passive_donor = SOM_passive(ilayer+1) / SOM_below                      ! Passive fraction in donor layer
                    move_amt = min(SOM_delta, SOM_below)                                        ! Move up only what donor can afford

                    SOM_fast(ilayer)     = SOM_fast(ilayer)     + move_amt * frac_fast_donor    ! Update fast pool (recipient layer)
                    SOM_slow(ilayer)     = SOM_slow(ilayer)     + move_amt * frac_slow_donor    ! Update slow pool (recipient layer)
                    SOM_passive(ilayer)  = SOM_passive(ilayer)  + move_amt * frac_passive_donor ! Update passive pool (recipient layer)

                    SOM_fast(ilayer+1)   = SOM_fast(ilayer+1)   - move_amt * frac_fast_donor    ! Update fast pool (donor layer)
                    SOM_slow(ilayer+1)   = SOM_slow(ilayer+1)   - move_amt * frac_slow_donor    ! Update slow pool (donor layer)
                    SOM_passive(ilayer+1)= SOM_passive(ilayer+1)- move_amt * frac_passive_donor ! Update passive pool (donor layer)

                else if (SOM_delta < 0.0_dp .and. SOM_current > 0.0_dp) then

                    !-----------------------------------------------------------!
                    ! Case 2: Downward SOM movement (excess in current layer)
                    !-----------------------------------------------------------!
                    frac_fast_donor    = SOM_fast(ilayer) / SOM_current                         ! Fast fraction in donor layer
                    frac_slow_donor    = SOM_slow(ilayer) / SOM_current                         ! Slow fraction in donor layer
                    frac_passive_donor = SOM_passive(ilayer) / SOM_current                      ! Passive fraction in donor layer
                    move_amt = min(-SOM_delta, SOM_current)                                     ! Move down only what donor can give

                    SOM_fast(ilayer)     = SOM_fast(ilayer)     - move_amt * frac_fast_donor    ! Update fast pool (donor layer)
                    SOM_slow(ilayer)     = SOM_slow(ilayer)     - move_amt * frac_slow_donor    ! Update slow pool (donor layer)
                    SOM_passive(ilayer)  = SOM_passive(ilayer)  - move_amt * frac_passive_donor ! Update passive pool (donor layer)

                    SOM_fast(ilayer+1)   = SOM_fast(ilayer+1)   + move_amt * frac_fast_donor    ! Update fast pool (recipient layer)
                    SOM_slow(ilayer+1)   = SOM_slow(ilayer+1)   + move_amt * frac_slow_donor    ! Update slow pool (recipient layer)
                    SOM_passive(ilayer+1)= SOM_passive(ilayer+1)+ move_amt * frac_passive_donor ! Update passive pool (recipient layer)

                end if

            end do

            !-------------------------------------------------------------------!
            ! Compute NEE: net exchange with atmosphere
            !-------------------------------------------------------------------!
            nee = resp_total - input_rate                                                       ! Net flux with atmosphere (kg C/m2)

            !-------------------------------------------------------------------!
            ! Mass conservation diagnostics
            !-------------------------------------------------------------------!
            mass_end = sum(SOM_fast(:)) + sum(SOM_slow(:)) + sum(SOM_passive(:)) + &
                sum(litter_struct_lig(:)) + sum(litter_struct_cel(:)) + &
                sum(litter_metabolic(:))                                                        ! Total C after update (kg C/m2)

            mass_error = abs(mass_end - mass_start)                                             ! Error magnitude

            if (mass_error > eps) then
                write(*,'(a,i5)') 'Mass conservation error at year: ', kyr
                write(*,'(a,f12.5)') 'Start mass = ', mass_start
                write(*,'(a,f12.5)') 'End mass   = ', mass_end
                write(*,'(a,f12.5)') 'Difference = ', mass_error
                stop 'Mass not conserved'
            end if

            !-------------------------------------------------------------------!
            ! Unrealistic fast/slow ratio divergence diagnostic check
            !-------------------------------------------------------------------!
            do ilayer = 1, nlayers - 1

                !---------------------------------------------------------------!
                ! Check fast/slow ratio mismatch
                !---------------------------------------------------------------!
                if (enable_ratio_check) then
                    if (SOM_slow(ilayer) > 0.0_dp .and. SOM_slow(ilayer+1) > 0.0_dp) then
                        if (abs((SOM_fast(ilayer)/SOM_slow(ilayer)) - (SOM_fast(ilayer+1)/SOM_slow(ilayer+1))) > 1.0_dp) then
                            write(*,'(a,i4)') 'Warning: Large fast/slow ratio mismatch between layers ', ilayer
                            write(*,'(a,f12.5)') '  L', real(ilayer, dp), '  Fast = ', SOM_fast(ilayer)
                            write(*,'(a,f12.5)') '  L', real(ilayer, dp), '  Slow = ', SOM_slow(ilayer)
                            write(*,'(a,f12.5)') '  L', real(ilayer, dp), '  Ratio = ', SOM_fast(ilayer)/SOM_slow(ilayer)
                            write(*,'(a,f12.5)') '  L', real(ilayer+1, dp), '  Fast = ', SOM_fast(ilayer+1)
                            write(*,'(a,f12.5)') '  L', real(ilayer+1, dp), '  Slow = ', SOM_slow(ilayer+1)
                            write(*,'(a,f12.5)') '  L', real(ilayer+1, dp), '  Ratio = ', SOM_fast(ilayer+1)/SOM_slow(ilayer+1)
                        end if
                    end if
                end if

                !---------------------------------------------------------------!
                ! Check slow/passive ratio mismatch
                !---------------------------------------------------------------!
                if (enable_ratio_check) then
                    if (SOM_passive(ilayer) > 0.0_dp .and. SOM_passive(ilayer+1) > 0.0_dp) then
                        if (abs((SOM_slow(ilayer)/SOM_passive(ilayer)) - (SOM_slow(ilayer+1)/SOM_passive(ilayer+1))) > 1.0_dp) then
                            write(*,'(a,i4)') 'Warning: Large slow/passive ratio mismatch between layers ', ilayer
                            write(*,'(a,f12.5)') '  L', real(ilayer, dp), '  Slow    = ', SOM_slow(ilayer)
                            write(*,'(a,f12.5)') '  L', real(ilayer, dp), '  Passive = ', SOM_passive(ilayer)
                            write(*,'(a,f12.5)') '  L', real(ilayer, dp), '  Ratio   = ', SOM_slow(ilayer)/SOM_passive(ilayer)
                            write(*,'(a,f12.5)') '  L', real(ilayer+1, dp), '  Slow    = ', SOM_slow(ilayer+1)
                            write(*,'(a,f12.5)') '  L', real(ilayer+1, dp), '  Passive = ', SOM_passive(ilayer+1)
                            write(*,'(a,f12.5)') '  L', real(ilayer+1, dp), '  Ratio   = ', SOM_slow(ilayer+1)/SOM_passive(ilayer+1)
                        end if
                    end if
                end if

            end do                                                                            ! End of ratio diagnostics per layer

        end do                                                                                ! End of sub-timestep loop

        !-----------------------------------------------------------------------!
        ! Annual diagnostics and output
        !-----------------------------------------------------------------------!
        final_depth = sum(SOM_fast(:) + SOM_slow(:) + SOM_passive(:)) / &
                      rho_SOM * 1000.0_dp                                                     ! Total SOM depth based on pool sum (mm)

        write(*,'(i5,36f12.5)') kyr, SOM_fast(:), SOM_slow(:), SOM_passive(:), &
            input_rate, resp_total, nee, final_depth

        write(unit_out,'(i0,",",36(f12.5,","),f12.5)') kyr, SOM_fast(:), &
            SOM_slow(:), SOM_passive(:), input_rate, resp_total, nee, final_depth

        do ilayer = 1, nlayers
            write(unit_layer,'(i0,",",i0,",",f12.5,",",f12.5,",",f12.5,",",f12.5,",",es12.5,",",es12.5,",",es12.5)') &
                kyr, ilayer, SOM_fast(ilayer), SOM_slow(ilayer), SOM_passive(ilayer), &
                SOM_fast(ilayer) + SOM_slow(ilayer) + SOM_passive(ilayer), &
                merge(SOM_fast(ilayer)/SOM_slow(ilayer), 0.0_dp, SOM_slow(ilayer) > 0.0_dp), &
                merge(SOM_slow(ilayer)/SOM_passive(ilayer), 0.0_dp, SOM_passive(ilayer) > 0.0_dp), &
                merge(SOM_fast(ilayer)/SOM_passive(ilayer), 0.0_dp, SOM_passive(ilayer) > 0.0_dp)
        end do

    end do  ! End of year loop

    !---------------------------------------------------------------------------!
    ! Final diagnostics
    !---------------------------------------------------------------------------!
    write(*,*)
    do ilayer = 1, nlayers                                                                      ! Depth of SOM in layer i (mm)
        write(*,'(a,i1,a,f12.5,a)') ' Layer ', ilayer, ' depth   : ', &
            (SOM_fast(ilayer)+SOM_slow(ilayer)+SOM_passive(ilayer))/rho_SOM * 1000.0_dp, ' mm'        
    
    end do

    write(*,*)
    write(*,'(a,f12.5,a)') ' Total SOM depth : ', final_depth, ' mm'                            ! SOM column depth (mm)
    write(*,'(a,f12.5,a)') ' Total SOM       : ', sum(SOM_fast(:) + SOM_slow(:) &
        + SOM_passive(:)), ' kg C/m2'                                                           ! Total SOM stock across all layers (kg C/m2)
    
    write(*,*)

    close(unit_out)                                                                             ! Close output file
    close(unit_layer)                                                                           ! Close output file
    write(*,*) 'Simulation complete. Results saved to ', outfile

end program toymodel_v1
!===============================================================================!
