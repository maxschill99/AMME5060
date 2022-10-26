program test

    ! Using modules
    ! USE jacobi
    USE variablemodule
    ! USE cgradient

    IMPLICIT NONE
    
    ! Initialising MPI Module
    include 'mpif.h'

    !-----------------------------------------------------------------------------------------------------!
    !-----------------------------------------------------------------------------------------------------!

    ! Initialising variables

    ! ! Variables newly defined in this program
    Real(kind = 8), allocatable :: T(:,:)
    Real(kind = 8), allocatable :: an(:,:), as(:,:), ae(:,:), aw(:,:), ap(:,:), b(:,:)
    Integer(kind = 4) :: il, ih, jl, jh, npp


    il = 0; ih = 0; jl = 0; jh = 0


    ! ! Calling subroutines from modules
    CALL intialise()

    ! call allocatevar()
    ! Allocation of variable sizes

    ! do j = jl,jh
    !     do i = il,ih
    !         allocate(an(il:ih,jl:jh))
    !         allocate(as(il:ih,jl:jh))
    !         allocate(ae(il:ih,jl:jh))
    !         allocate(aw(il:ih,jl:jh))
    !         allocate(ap(il:ih,jl:jh))
    !         allocate(b(il:ih,jl:jh))
    !         allocate(T(il:ih,jl:jh))
    !     end do
    ! end do

    ! INITIALISE MPI
	CALL MPI_INIT(ierr)
	CALL MPI_COMM_RANK(MPI_COMM_WORLD, pid, ierr) ! Getting processor ID number
	CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr) ! Getting number of total processors in global communicator



    ! Assume that nx = ny
    npp = nx/nprocs
    il = pid*npp + 1
    ih = il + npp - 1
    jl = il
    jh = ih

    write(*,*) pid, il,ih, jl,jh

    ! !-----------------------------------------------------------------------------------------------------!
    ! !-----------------------------------------------------------------------------------------------------!
    ! Initialising boundary conditions on temp array

    ! !-----------------------------------------------------------------------------------------------------!
    ! !-----------------------------------------------------------------------------------------------------!
    ! Computation

    ! Calculation of solution using only jacobi solver
    ! call solninit(an,as,ae,aw,ap,b,T,il,ih,jl,jh)

    ! write(*,1600) T
    ! call jacobisolv(an,as,ae,aw,ap,b,x,y,T)


    ! Calculation of solution using Conjugate Gradient Method
    ! call solninit(an,as,ae,aw,ap,b)
    ! call CGSolve(an,as,ae,aw,ap,b,T)
    
    ! write(*,1600) T
    
    
    ! NOTE: Need to update this to match the number of spatial divisions (nx)
    1100 FORMAT(8(F20.10,1x))
    1600 FORMAT(5(F14.8,1x))
    1400 FORMAT(8(F14.8,1x))
    1200 FORMAT(I2.1,I2.1,6(F8.4,1x))
    1300 FORMAT(I2.1,6(F10.8,1x))


    CALL MPI_FINALIZE(ierr)

end program test