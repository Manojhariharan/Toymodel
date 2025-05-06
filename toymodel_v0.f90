!======================================================================!
! Program: toymodel_v0
! Purpose: Minimal 1-layer, 1-pool SOM model with input and respiration
!======================================================================!
program toymodel_v0
    implicit none

    !------------------------------------------------------------------!
    ! Parameters
    !------------------------------------------------------------------!
    integer, parameter :: nyears = 6000                ! Simulation length (years)
    integer, parameter :: nlayers = 1                  ! Number of layers
    real, parameter :: dt = 1.0                        ! Timestep (year)
    real, parameter :: input_rate = 1.05               ! SOM input (kg C/m2/year)
    real, parameter :: k_decay = 0.007                 ! Updated decay rate (/year)

    !------------------------------------------------------------------!
    ! Variables
    !------------------------------------------------------------------!
    real :: input, decay, respiration
    real, dimension(nlayers) :: SOM                    ! SOM pool (kg C/m2)
    integer :: t, i

    !------------------------------------------------------------------!
    ! Initialization
    !------------------------------------------------------------------!
    SOM = 0.0                                          ! Initial SOM for each layer
    write(*,'(a)') 'Year    SOM_C     Input    Respired'

    !------------------------------------------------------------------!
    ! Time loop
    !------------------------------------------------------------------!
    do t = 1, nyears
        do i = 1, nlayers
            input = input_rate
            decay = k_decay * SOM(i)
            respiration = decay                        ! All decay becomes CO2

            SOM(i) = SOM(i) + input - respiration
        end do

        write(*,'(i4,3f10.4)') t, SOM(1), input, respiration
    end do

end program toymodel_v0
!======================================================================!
