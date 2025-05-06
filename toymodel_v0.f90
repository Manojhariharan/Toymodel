!======================================================================!
! Program: toymodel_v0
! Purpose: Minimal 1-layer, 1-pool SOM model with input and respiration
!          includes check for conservation of carbon mass
!======================================================================!
program toymodel_v0
    implicit none

    !------------------------------------------------------------------!
    ! Parameters
    !------------------------------------------------------------------!
    integer, parameter :: nyears = 6000                                ! Simulation length (years)
    integer, parameter :: nlayers = 1                                  ! Number of layers
    real, parameter :: dt = 1.0                                        ! Timestep (year)
    real, parameter :: input_rate = 1.05                               ! SOM input (kg C/m2/year)
    real, parameter :: k_decay = 0.007                                 ! Decay rate (/year)

    !------------------------------------------------------------------!
    ! Variables
    !------------------------------------------------------------------!
    real :: input, decay, respiration                                  ! Fluxes per year
    real :: total_input, total_respired, final_SOM                     ! Totals for mass balance
    real, dimension(nlayers) :: SOM                                    ! SOM pool (kg C/m2)
    integer :: t, i

    !------------------------------------------------------------------!
    ! Initialization
    !------------------------------------------------------------------!
    SOM = 0.0                                                          ! Initial SOM for each layer
    total_input = 0.0                                                  ! Initialize input sum
    total_respired = 0.0                                               ! Initialize respiration sum
    write(*,'(a)') 'Year    SOM_C     Input    Respired'

    !------------------------------------------------------------------!
    ! Time loop
    !------------------------------------------------------------------!
    do t = 1, nyears
        do i = 1, nlayers
            input = input_rate                                         ! Constant litter input
            decay = k_decay * SOM(i)                                   ! First-order decay
            respiration = decay                                        ! All decay becomes CO2

            SOM(i) = SOM(i) + input - respiration                      ! Update SOM pool

            total_input = total_input + input                          ! Accumulate total input
            total_respired = total_respired + respiration              ! Accumulate total respired
        end do

        write(*,'(i4,3f10.4)') t, SOM(1), input, respiration
    end do

    !------------------------------------------------------------------!
    ! Mass conservation check
    !------------------------------------------------------------------!
    final_SOM = SOM(1)                                                 ! Final SOM stock
    write(*,*)
    write(*,'(a,f12.4)') 'Total C input:      ', total_input
    write(*,'(a,f12.4)') 'Total C respired:   ', total_respired
    write(*,'(a,f12.4)') 'Final SOM stock:    ', final_SOM
    write(*,'(a,f12.4)') 'Input - (resp + SOM): ', total_input - (total_respired + final_SOM)

end program toymodel_v0
!======================================================================!
