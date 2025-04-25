!======================================================================!
! Program: toymodel_v0
! Purpose: Minimal 1-layer, 1-pool SOM model with input and respiration
!======================================================================!

program toymodel_v0
    implicit none

    !------------------------------------------------------------------!
    ! Parameters
    !------------------------------------------------------------------!
    integer, parameter :: nyears = 100               ! Simulation length (years)
    real, parameter :: dt = 1.0                      ! Timestep (year)
    real, parameter :: input_rate = 0.2              ! SOM input (kg C/m2/year)
    real, parameter :: k_decay = 0.05                ! Decay rate (/year)

    !------------------------------------------------------------------!
    ! Variables
    !------------------------------------------------------------------!
    real :: SOM                                      ! Soil organic matter (kg C/m2)
    real :: input, decay, respiration
    integer :: t

    !------------------------------------------------------------------!
    ! Initialization
    !------------------------------------------------------------------!
    SOM = 1.0                                         ! Initial SOM (kg C/m2)
    write(*,'(a)') 'Year    SOM_C     Input    Respired'

    !------------------------------------------------------------------!
    ! Time loop
    !------------------------------------------------------------------!
    do t = 1, nyears
        input = input_rate
        decay = k_decay * SOM
        respiration = decay                          ! All decayed C is respired

        SOM = SOM + input - respiration

        write(*,'(i4,3f10.4)') t, SOM, input, respiration
    end do

end program toymodel_v0
!======================================================================!