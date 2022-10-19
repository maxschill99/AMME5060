program test

    ! Using modules
    USE jacobi
    USE variablemodule

    IMPLICIT NONE

    !-----------------------------------------------------------------------------------------------------!
    !-----------------------------------------------------------------------------------------------------!

    ! Initialising variables

    ! ! Variables newly defined in this program
    Real(kind = 8), allocatable :: Tin(:,:), Tout(:,:)
    Real(kind = 8), allocatable :: an(:,:), as(:,:), ae(:,:), aw(:,:), ap(:,:), b(:,:)


    ! ! Calling subroutines from modules
    CALL intialise()

    ! call allocatevar()
    ! Allocation of variable sizes
    allocate(an(nx,ny))
    allocate(as(nx,ny))
    allocate(ae(nx,ny))
    allocate(aw(nx,ny))
    allocate(ap(nx,ny))
    allocate(b(nx,ny))
    allocate(x(nx))
    allocate(y(ny))
    allocate(Tin(nx,ny))
    allocate(Tout(nx,ny))

    ! !-----------------------------------------------------------------------------------------------------!
    ! !-----------------------------------------------------------------------------------------------------!
    ! Computation

    call solninit(an,as,ae,aw,ap,b,Tin,Tout)
    call solver(an,as,ae,aw,ap,b,x,y,Tin,Tout)

    write(*,1600) Tout

    
    
    
    ! NOTE: Need to update this to match the number of spatial divisions (nx)
    1100 FORMAT(8(F20.10,1x))
    1600 FORMAT(5(F14.8,1x))
    1400 FORMAT(8(F14.8,1x))
    1200 FORMAT(I2.1,I2.1,6(F8.4,1x))
    1300 FORMAT(I2.1,6(F10.8,1x))

end program test