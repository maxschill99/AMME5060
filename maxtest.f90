program test

    ! Using modules
    USE jacobi
    USE variablemodule
    USE solinit
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
    Integer(kind = 4) :: il, ih, jl, jh, npp, iter

    ! Solution solver variables
    Real(kind = 8) :: rcurrent, rc, time
    Real(kind = 8), allocatable :: Told(:,:), Tn(:,:), resmat(:,:)

    il = 0; ih = 0; jl = 0; jh = 0


    ! ! Calling subroutines from modules
    CALL intialise()

    ! !-----------------------------------------------------------------------------------------------------!
    ! !-----------------------------------------------------------------------------------------------------!
    ! INITIALISE MPI
	CALL MPI_INIT(ierr)
	CALL MPI_COMM_RANK(MPI_COMM_WORLD, pid, ierr) ! Getting processor ID number
	CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr) ! Getting number of total processors in global communicator

    ! !-----------------------------------------------------------------------------------------------------!
    ! !-----------------------------------------------------------------------------------------------------!
    ! PARTITIONING -- WILL BE COVERED BY CLARA'S STUFF
    ! Assume that nx = ny
    if (nprocs==1) then
        npp = nx/nprocs
        il = pid*npp + 1
        ih = il + npp - 1
        jl = il
        jh = ih
    end if

    if (nprocs==2) then
        if (pid==0) then
            il = 1
            ih = ny/2
            jl = 1
            jh = nx/2
        else
            il = ny/2 + 1
            ih = ny
            jl = nx/2 + 1
            jh = nx
        end if
    end if
    write(*,*) pid, il,ih, jl,jh

    ! !-----------------------------------------------------------------------------------------------------!
    ! !-----------------------------------------------------------------------------------------------------!
    ! Allocation of variable sizes
    call allocatevars(an,as,ae,aw,ap,b,T,Told,Tn,resmat,il,ih,jl,jh)

    ! Initialising boundary conditions on temp array
    ! Doesnt at the moment with processors more than 1 because the partitioning is wrong
    call solutioninit(an,as,ae,aw,ap,b,T,il,ih,jl,jh)
    write(*,1600) T

    ! !-----------------------------------------------------------------------------------------------------!
    ! !-----------------------------------------------------------------------------------------------------!
    ! Computation

    ! Solution loop variable initialisation
    ! Inititialising old temperature array
    Told(:,:) = 0

    ! Initialising time counter
    time = 0
    iter = 0

    ! Begin the solution loop
    ! do while ((t<ttot).and.(r>rmax))
    do while (time<t_final)


        ! Calculation of solution using only jacobi solver
        call jac(an,as,ae,aw,ap,b,T,il,ih,jl,jh,time)

        ! COMMUNICATION
            

        ! computing residuals
        call respar(aw,ae,an,as,ap,b,T,il,ih,jl,jh,resmat)

        ! Calculate Domain averaged residual for stopping critterion
        rcurrent = SUM(SUM(ABS(res(il:ih,jl:jh)),1),1) / ((nx-2)*(ny-2))

        ! Combining all residuals on each processor and sending to processor 0
        call MPI_REDUCE(rcurrent,rc,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

        if (pid == 0) then
            rc = rc/nprocs
        end if

        call MPI_BCAST(rc,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)


        ! Printing to screen after a certain number of iterations
        if (mod(iter,100).eq.0) then
            write(*,*) '      iter', '      res'
            write(*,*) iter, rcurrent
        end if

        ! updating old and new temperature values
        Told = T
        T = Tn

        ! updating counter
        time = time + dt
        iter = iter + 1

    end do






    ! Calculation of solution using Conjugate Gradient Method
    ! call CGSolve(an,as,ae,aw,ap,b,T)

    ! !-----------------------------------------------------------------------------------------------------!
    ! !-----------------------------------------------------------------------------------------------------!
    !! Getting final total temp array
    ! Gathering all temperature areas to processor 0
    ! Ttemp = T(1:ny,il:ih)
    ! call MPI_GATHER(Ttemp,ny*npp,MPI_DOUBLE_PRECISION,ttot,ny*npp,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    
    ! write(*,1600) T
    
    
    ! NOTE: Need to update this to match the number of spatial divisions (nx)
    1100 FORMAT(8(F20.10,1x))
    1600 FORMAT(5(F14.8,1x))
    1400 FORMAT(8(F14.8,1x))
    1200 FORMAT(I2.1,I2.1,6(F8.4,1x))
    1300 FORMAT(I2.1,6(F10.8,1x))


    CALL MPI_FINALIZE(ierr)

end program test