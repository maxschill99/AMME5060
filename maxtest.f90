program test

    ! Using modules
    USE jacobi
    USE variablemodule

    IMPLICIT NONE

    !-----------------------------------------------------------------------------------------------------!
    !-----------------------------------------------------------------------------------------------------!

    ! Initialising variables

    ! Variables newly defined in this program
    Integer :: it,jt
    Real(kind = 4), ALLOCATABLE :: Tin(:,:), Tout(:,:)

    ! Calling subroutines from modules
    CALL intialise()

    !-----------------------------------------------------------------------------------------------------!
    !-----------------------------------------------------------------------------------------------------!
    ALLOCATE(Tin(nx,ny))
    ALLOCATE(Tout(nx,ny))
    ALLOCATE(x(nx))

    Tin(:,:) = 0
    
    do i = 1,nx
        x = (i-1)*dx
        Tin(1,i) = sin((pi*x(i))/Lx)
    end do
    Tin(:,1) = 0
    Tin(:,nx) = 0
    Tin(ny,:) = 0

    call solver(Tin,Tout)

    ! write(*,1600) Tout

    
    
    
    ! NOTE: Need to update this to match the number of spatial divisions (nx)
    1100 FORMAT(8(F20.10,1x))
    1600 FORMAT(33(F14.8,1x))
    1400 FORMAT(8(F14.8,1x))
    1200 FORMAT(I2.1,I2.1,6(F8.4,1x))
    1300 FORMAT(I2.1,6(F10.8,1x))

end program test