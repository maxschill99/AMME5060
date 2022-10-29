PROGRAM ASSIGN1

    IMPLICIT NONE
	INCLUDE 'mpif.h'

	! Variable declartion
	INTEGER :: pid, nprocs, ierr, tag, stat1(MPI_STATUS_SIZE), stat2(MPI_STATUS_SIZE), root, req1, req2, req3, req4

	! Counter and matrix variables
    REAL (KIND = 8) :: a, Lx, Ly, dx, dy, x, y, pi
    INTEGER :: n, i, j, nx, ny
    REAL (KIND = 8), ALLOCATABLE :: T(:,:), ttot(:,:), Ttemp(:,:), Tnew(:,:)
    REAL (KIND = 8), ALLOCATABLE :: U(:,:), V(:,:), Utemp(:,:), Vtemp(:,:), Utot(:,:), Vtot(:,:)

    ! Array distribution variables
    INTEGER :: nlocal, il, ih, npp

    ! Solver calculations
    REAL (KIND = 8) :: rmax, k, itermax, iter, res, rcurrent
    REAL (KIND = 8) :: topa, topb, topc, topd, bota
    REAL (KIND = 8) :: Tmax, Tave, Ttave
    REAL (KIND = 8), ALLOCATABLE :: aw(:,:), ae(:,:), an(:,:), as(:,:), ap(:,:)
    REAL (KIND = 8), ALLOCATABLE :: r(:,:), rn(:,:)
    REAL (KIND = 8) :: rcurrentc, resc

    ! Writing out for plotting
    REAL (KIND = 8), ALLOCATABLE :: xdat(:), ydat(:), Tdat(:)

    ! CPU timer
    REAL (KIND = 8) :: time1, time2, time
    
    CALL cpu_time(time1)
		
	CALL MPI_INIT(ierr) 					         ! Initialize MPI Program
	CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)   ! Get Size of Communicator (Number of procs)
	CALL MPI_COMM_RANK(MPI_COMM_WORLD,pid,ierr) 	 ! Get Rank of Communicator (proc id)

    !-----------------------------------------------------------------------------------
    ! POSITIONING OF DATA WHEN PRINTING TO SCREEN CAN BE SEEN AS BELOW

    ! |---> rows
    ! |
    ! Y
    ! columns

    ! 1   X   X   X   X   X   .   .   .   
    ! 2   X   X   X   X   X   .   .   .
    ! 3   X   X   X   X   X   .   .   .
    ! 4   X   X   X   X   X   .   .   .
    ! .   X   X   X   X   X   .   .   . 
    ! .   X   X   X   X   X   .   .   .
    ! .   X   X   X   X   X   .   .   .
    ! n   X   X   X   X   X   .   .   .
    !     1   2   3   4   .   .   .   n   

    ! columns going up and down
    ! rows going from left to right
    
    
    !-----------------------------------------------------------------------------------
    ! INITIAL VALUES AND CHECKS

    a = 80
    Lx = 1
    Ly = 1
    pi = (4.D0*DATAN(1.D0))

    ! number of spatial divisions
    nx = 128
    ny = 128

    ! Discretisation in each direction
	dx = Lx/(nx-1)
    dy = Ly/(ny-1)

    ! Running checks
    IF (pid.eq.0) THEN
        IF ((nx).LT.nprocs) THEN
            WRITE(*,*) 'Number of grid points is less than the number of processors'
            STOP
        ELSE IF (MOD(nx,nprocs).ne.0) THEN
            WRITE(*,*) 'nx is not divisible amongst processors'
            STOP
        END IF 
    END IF

    !-----------------------------------------------------------------------------------
    ! GLOBAL/LOCAL INDICES AND INITIAL TEMP DISTRIBUTION

    npp = nx/nprocs
    il = pid*npp + 1
    ih = il + npp - 1

    ! Initialize temp distribution
    ! ny is the number of points in the y-direction, i.e. number of rows
    ! nx is the number of points in the x-direction, i.e. number of columns
    ALLOCATE(T(0:ny+1,il-1:ih+1))

    ! Initialising total areas on processor 0
    if (pid==0) then
        ALLOCATE(ttot(ny,nx))
        ALLOCATE(Utot(ny,nx))
        ALLOCATE(Vtot(ny,nx))
    end if
  

    ! Calculating initial temp profile and setting boundary conditions
    do i = il-1,ih+1
        x = (i-1)*dx
        T(1,i) = SIN((pi*x)/Lx)

        T(ny,i) = 0 
    
        if (pid == 0) then
            T(:,il) = 0
        elseif (pid == nprocs-1) then
            T(:,ih) = 0
        end if

    end do

    ! Creating a temp array to combine all processor data
    ALLOCATE(Ttemp(ny,il:ih))

    ! RUNNING A CHECK TO SEE IF INTIIAL TEMP IS BEING CALCULATED CORRECTLY
    Ttemp = T(1:ny,il:ih)

    ! Gathering data to processor 1
    CALL MPI_GATHER(Ttemp,ny*npp,MPI_DOUBLE_PRECISION,ttot,ny*npp,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

     
    if (pid==0) then
        ! output data into a file 
        write(*,*) 'INITAL TEMP DISTRIBUTION'
        write(*,1100) ttot
    end if
    

    !-----------------------------------------------------------------------------------
    ! INITIALISING VELOCITY FIELDS

    ALLOCATE(U(1:ny,il:ih))
    ALLOCATE(V(1:ny,il:ih))
    ! Calculating fluid velocity fields
    do i = il,ih
        x = (i-1)*dx

        do j = 1,ny
            y = (ny-j)*dy

            U(j,i) = SIN((pi*x)/Lx)*COS((pi*y)/Ly) + SIN((2*pi*x)/Lx)*COS((2*pi*y)/Ly) + SIN((4*pi*x)/Lx)*COS((4*pi*y)/Ly)
            V(j,i) = -COS((pi*x)/Lx)*SIN((pi*y)/Ly) - COS((2*pi*x)/Lx)*SIN((2*pi*y)/Ly) - COS((4*pi*x)/Lx)*SIN((4*pi*y)/Ly)
        end do
    end do
    
    ! ALLOCATE(Utemp(ny,il:ih))
    ! ALLOCATE(Vtemp(ny,il:ih))
    ! Utemp = U(1:ny,il:ih)
    ! Vtemp = V(1:ny,il:ih)

    ! ! Gathering data to processor 0
    ! CALL MPI_GATHER(Utemp,ny*npp,MPI_DOUBLE_PRECISION,Utot,ny*npp,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    ! CALL MPI_GATHER(Vtemp,ny*npp,MPI_DOUBLE_PRECISION,Vtot,ny*npp,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

    ! if (pid==0) then
    !     write(*,*) 'U velocities'
    !     write(*,1100) Utot
    !     write(*,*) 'V velocities'
    !     write(*,1100) Vtot
    ! end if


    !-----------------------------------------------------------------------------------
    ! SOLVER

    ! Creating the A matrix coefficients
    ALLOCATE(aw(ny,il:ih))
    ALLOCATE(ae(ny,il:ih))
    ALLOCATE(an(ny,il:ih))
    ALLOCATE(as(ny,il:ih))
    ALLOCATE(ap(ny,il:ih))

    ! Computing A matrix coefficients
    do i = il,ih
        do j = 1,ny
            aw(j,i) = -((a*U(j,i))/(2*dx) + 1/dx**2)
            ae(j,i) = ((a*U(j,i))/(2*dx) - 1/dx**2)
            as(j,i) = -((a*V(j,i))/(2*dy) + 1/dy**2)
            an(j,i) = ((a*V(j,i))/(2*dy) - 1/dy**2)
            ap(j,i) = -(2/dx**2 + 2/dy**2)
        end do
    end do

    write(*,*) an(1,1), as(1,1), ae(1,1), aw(1,1), ap(1,1)

    ! Creating a new temp array to test residuals
    ALLOCATE(Tnew(0:ny+1,il-1:ih+1))
    Tnew = T

    ! Communication tag
    tag = 101

    ! Initialising residuals and counter
    rmax = 0.000000001
    itermax = 1000000
    iter = 0
    res = 1
    resc = 1

    ALLOCATE(r(ny,il:ih))
    ALLOCATE(rn(ny,il:ih))
    ! Initialising residual array
    do i = il,ih
        do j = 1,ny
            r(j,i) = (T(j,i-1)*aw(j,i) + T(j,i+1)*ae(j,i) + T(j-1,i)*an(j,i) + T(j+1,i)*as(j,i))/ap(j,i)
        end do
    end do   

    ! Calculate Domain averaged residual for stopping critterion
    rcurrent = SUM(SUM(ABS(r(1:ny,il:ih)),1)) / ((ih-il+1)*(ny))

    write(*,*) rcurrent

    ! STEADY STATE COMPUTATION
    do while ((res > rmax).and.(iter <= itermax))
        !---------------------------------------------------------------------------
        ! SOLVER

        ! RED BLACK SOLVER
        ! CALCULATING RED NODES
        do i = il,ih
            do j = 2,ny-1,2
                Tnew(j,i) = (T(j,i-1)*aw(j,i) + T(j,i+1)*ae(j,i) + T(j-1,i)*an(j,i) + T(j+1,i)*as(j,i))/ap(j,i) 
            end do
        end do

        ! Non-blocking processor communication
        if (pid.NE.nprocs-1) then ! not on far right
            call MPI_IRECV(Tnew(:,ih+1),ny,MPI_DOUBLE_PRECISION,pid+1,tag,MPI_COMM_WORLD,req1,ierr) ! receive from the right
        end if
        if (pid.NE.0) then ! not on far left
            call MPI_ISEND(Tnew(:,il),ny,MPI_DOUBLE_PRECISION,pid-1,tag,MPI_COMM_WORLD,req2,ierr) ! send to the left
        end if
        if (pid.NE.0) then ! not on far left
            call MPI_IRECV(Tnew(:,il-1),ny,MPI_DOUBLE_PRECISION,pid-1,tag,MPI_COMM_WORLD,req1,ierr) ! receive from the left
        end if
        if (pid.NE.nprocs-1) then ! not on far right
            call MPI_ISEND(Tnew(:,ih),ny,MPI_DOUBLE_PRECISION,pid+1,tag,MPI_COMM_WORLD,req2,ierr)	! send to the right
        end if

        ! Wait until all processors have finished communicating
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        ! CALCULATING BLACK NODES
        do i = il,ih
            do j = 3,ny-1,2
                Tnew(j,i) = (T(j,i-1)*aw(j,i) + T(j,i+1)*ae(j,i) + T(j-1,i)*an(j,i) + T(j+1,i)*as(j,i))/ap(j,i)
            end do
        end do

        ! Non-blocking processor communication
        if (pid.NE.nprocs-1) then ! not on far right
            call MPI_IRECV(Tnew(:,ih+1),ny,MPI_DOUBLE_PRECISION,pid+1,tag,MPI_COMM_WORLD,req1,ierr) ! receive from the right
        end if
        if (pid.NE.0) then ! not on far left
            call MPI_ISEND(Tnew(:,il),ny,MPI_DOUBLE_PRECISION,pid-1,tag,MPI_COMM_WORLD,req2,ierr) ! send to the left
        end if
        if (pid.NE.0) then ! not on far left
            call MPI_IRECV(Tnew(:,il-1),ny,MPI_DOUBLE_PRECISION,pid-1,tag,MPI_COMM_WORLD,req1,ierr) ! receive from the left
        end if
        if (pid.NE.nprocs-1) then ! not on far right
            call MPI_ISEND(Tnew(:,ih),ny,MPI_DOUBLE_PRECISION,pid+1,tag,MPI_COMM_WORLD,req2,ierr)	! send to the right
        end if

        ! Wait until all processors have finished communicating
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        !---------------------------------------------------------------------------
        ! RESIDUALS
        ! Updating temperature array, residual and iteration counter

        do i = il,ih
            do j = 1,ny
                rn(j,i) = (Tnew(j,i-1)*aw(j,i) + Tnew(j,i+1)*ae(j,i) + Tnew(j-1,i)*an(j,i) + Tnew(j+1,i)*as(j,i))/ap(j,i)
            end do
        end do

        rcurrentc = SUM(ABS(Tnew(1:ny,il:ih) - T(1:ny,il:ih))) / ((ih-il+1)*(ny))
        rcurrent = SUM(ABS(rn - r))/ ((ih-il+1)*ny)

        ! Combining all residuals on each processor and sending to processor 0
        call MPI_REDUCE(rcurrent,res,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(rcurrentc,resc,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

        if (pid == 0) then
            res = res/nprocs
            resc = resc/nprocs
        end if

        call MPI_BCAST(res,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(resc,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

        iter = iter + 1
        T = Tnew
        r = rn

    end do

    !-----------------------------------------------------------------------------------
    ! SYNCHRONIZE AND GATHER I/O

    ! Gathering all temperature areas to processor 0
    Ttemp = T(1:ny,il:ih)
    call MPI_GATHER(Ttemp,ny*npp,MPI_DOUBLE_PRECISION,ttot,ny*npp,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

    ! Computing average temperature across each processors domain
    Tave = sum(Ttemp)/(ny*(ih-il+1))
    ! Sending all average temperatures to processor 0 and summing
    call MPI_REDUCE(Tave, Ttave, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    ! Computing the global average domain temperature based on the number of processors
    Ttave = Ttave/nprocs

    if (pid==0) then
        write(*,*) 'FINAL TEMP DISTRIBUTION'
        ! write(*,1100) ttot
        Tmax = maxval(ttot)
        write(*,*) 'Number of x and y grids =', nx, ny
        write(*,*) 'Maximum number of iterations =', itermax
        write(*,*) 'Residual limit =', rmax
        write(*,*) 'Maximum Temperature =', Tmax
        write(*,*) 'Average Domain Temperature =', Ttave
        write(*,*) 'Number of iterations =', iter
        write(*,*) 'Residual =', res
        write(*,*) 'Residual on processor 0 =', rcurrent
        write(*,*) 'Residual check =', resc
        write(*,*) 'Residual check on processor 0 =', rcurrentc
    end if

    !-----------------------------------------------------------------------------------
    ! WRITING OUT TO TECPLOT FILE

    if (pid == 0) then

        ALLOCATE(xdat(nx))
        ALLOCATE(ydat(ny))

        do i = 1,nx
            x = (i-1)*dx
            xdat(i) = x
        end do
            
        do j = 1,ny
            y = (ny-j)*dy
            ydat(j) = y
        end do
    
        CALL tecplot_2D(1,nx,ny,xdat,ydat,Ttot)

    end if
	
    !-----------------------------------------------------------------------------------
    ! NOTE: Need to update this to match the number of spatial divisions (ny)
    1100 FORMAT(8(F20.10,1x))
    1600 FORMAT(4(F14.8,1x))
    1400 FORMAT(8(F14.8,1x))
    1200 FORMAT(I2.1,I2.1,6(F8.4,1x))
    1300 FORMAT(I2.1,6(F10.8,1x))

    ! FINISHING CPU TIMER CALCS
    ! stop timer
	call cpu_time(time2)

    ! write(*,*) pid, time2-time1

    ! Sending all cpu times to processor to sum for total time
    call MPI_REDUCE(time2-time1, time, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
	
	! write computational time to processor 0
	if (pid == 0) then 
		write(*,*)"Computational time =", time
        ! outputting data for analysis
        call data_out(2,Tmax,Ttave,iter,res,time)
	end if 

    CALL MPI_FINALIZE(ierr)

END PROGRAM ASSIGN1



! SUBROUTINES
!----------------------------------------------------------------------------- 
subroutine tecplot_2D ( iunit, nx, ny, x, y, T )

  IMPLICIT NONE

  Integer ( kind = 4 ) iunit,nx,ny,i,j,ierr
  Real ( kind = 8 ) T(ny,nx)
  Real ( kind = 8 ) x(nx)
  Real ( kind = 8 ) y(ny)

  Character(80), parameter ::  file_name = 'TecPlot2Dnew.tec'

  open ( unit = iunit, file = file_name, form = 'formatted', access = 'sequential', status = 'replace', iostat = ierr )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) '  Error opening file : tecplot_2D '
    stop
  end if
   
  write ( iunit, '(a)' ) 'Title="' // trim ( 'Temperature Data' ) // '"'
  write ( iunit, '(a)' ) 'Variables=' // trim ( '"X","Y","T"' )

  write ( iunit, '(a)' ) ' '
  write ( iunit, '(a,i6,a,i6,a)' ) 'Zone I=', ny, ', J=', nx, ', F=POINT'
 
  do i = 1, nx
    do j = 1, ny
      write ( iunit, '(2f10.3,g15.6)' ) x(i), y(j), T(j,i)
    end do
  end do
  
  close ( unit = iunit )

end subroutine tecplot_2D
!----------------------------------------------------------------------------- 
subroutine data_out ( iunit, Tmax, Tave, Iter, Res, time )

  IMPLICIT NONE

  Integer ( kind = 4 ) iunit, ierr
  Real ( kind = 8 ) Tmax, Tave, Iter, Res, time

  Character(80), parameter ::  file_name = 'data_out.tec'

  open ( unit = iunit, file = file_name, form = 'formatted', access = 'sequential', status = 'replace', iostat = ierr )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) '  Error opening file : tecplot_2D '
    stop
  end if
   
  write ( iunit, '(a)' ) 'Title="' // trim ( 'Output Data' ) // '"'
  write ( iunit, '(a)' ) 'Variables=' // trim ( 'Tmax, Tave, Iter, Res, time' )

  write ( iunit, '(5e30.10)' ) Tmax, Tave, Iter, Res, time
  
  close ( unit = iunit )

end subroutine data_out