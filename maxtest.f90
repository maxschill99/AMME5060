program test

    ! Using modules
    USE jacobi
    USE variablemodule

    IMPLICIT NONE

    ! Calling subroutines from modules
    CALL intialise()

    !-----------------------------------------------------------------------------------------------------!
    !-----------------------------------------------------------------------------------------------------!

    ! Initialising variables
    Real(kind = 8) :: Tin(nx,ny), Tout(nx,ny)
    Integer :: it,jt

    !-----------------------------------------------------------------------------------------------------!
    !-----------------------------------------------------------------------------------------------------!

    ! Inittial temperature distribution

    !                 f(x)
    !        ------------------------
    !        |                       |
    !        |                       |
    !        |                       |
    !    0   |                       |   0
    !        |                       |
    !        |                       |
    !        |                       |
    !        ------------------------
    !                   0


    !-----------------------------------------------------------------------------------------------------!
    !-----------------------------------------------------------------------------------------------------!

    do it = 1,nx
        ! Tin(1,i) = sin((pi*(i-1))/Lx)
        
    end do
    ! Tin(:,1) = 0
    ! Tin(:,nx) = 0
    ! Tin(ny,:) = 0

    call solver(Tin,Tout)

end program test