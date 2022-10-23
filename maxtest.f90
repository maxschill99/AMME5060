program test

    ! Using modules
    USE jacobi
    USE variablemodule
    USE cgradient

    IMPLICIT NONE
    
    ! Initialising MPI Module
    include 'mpif.h'

  	! INITIALISE MPI
	CALL MPI_INIT(ierr)
	CALL MPI_COMM_RANK(MPI_COMM_WORLD, pid, ierr) ! Getting processor ID number
	CALL MPI_COMM_SIZE(MPI_COMM_WORLD, Nprocs, ierr) ! Getting number of total processors in global communicator

    !-----------------------------------------------------------------------------------------------------!
    !-----------------------------------------------------------------------------------------------------!

    ! Initialising variables

    ! ! Variables newly defined in this program
    Real(kind = 8), allocatable :: T(:,:)
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
    allocate(T(nx,ny))

    ! !-----------------------------------------------------------------------------------------------------!
    ! !-----------------------------------------------------------------------------------------------------!
    ! Computation

    ! Calculation of solution using only jacobi solver
    ! call solninit(an,as,ae,aw,ap,b,T)
    ! call jacobisolv(an,as,ae,aw,ap,b,x,y,T)

    ! write(*,1600) T


    ! Calculation of solution using Conjugate Gradient Method
    call solninit(an,as,ae,aw,ap,b,T)
    call CGSolve(an,as,ae,aw,ap,b,T)
    
    write(*,1600) T
    
    
    ! NOTE: Need to update this to match the number of spatial divisions (nx)
    1100 FORMAT(8(F20.10,1x))
    1600 FORMAT(5(F14.8,1x))
    1400 FORMAT(8(F14.8,1x))
    1200 FORMAT(I2.1,I2.1,6(F8.4,1x))
    1300 FORMAT(I2.1,6(F10.8,1x))


    CALL MPI_FINALIZE(ierr)

end program test