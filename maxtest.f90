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
    ! Real(kind = 8), allocatable :: T(:,:)
    ! Real(kind = 8), allocatable :: an(:,:), as(:,:), ae(:,:), aw(:,:), ap(:,:), b(:,:)
    ! Integer(kind = 8) :: il, ih, jl, jh, npp, iter

    ! ! Solution solver variables
    ! Real(kind = 8) :: rcurrent, rc, time
    ! Real(kind = 8), allocatable :: Told(:,:), Tn(:,:), resmat(:,:)

    ! ! Gathering variables for final solution
    ! Real(kind = 8), allocatable :: Ttemp(:,:), Ttot(:,:)
    ! Real(kind = 8) :: numcount
    Integer:: req1, req2


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
            ih = ny
            jl = 1
            jh = nx/2 + 1
            ! il = 1
            ! ih = 500
            ! jl = 1
            ! jh = 251

            resil = 1
            resih = ny
            resjl = 1
            resjh = nx/2
            ! resil = 1
            ! resih = 500
            ! resjl = 1
            ! resjh = 250
        else
            il = 1
            ih = ny
            jl = nx/2
            jh = nx
            ! il = 1
            ! ih = 500
            ! jl = 250
            ! jh = 500

            resil = 1
            resih = ny
            resjl = nx/2+1
            resjh = nx
            ! resil = 1
            ! resih = 500
            ! resjl = 251
            ! resjh = 500
        end if
    end if
    ! write(*,*) pid, il,ih, jl,jh

    ! !-----------------------------------------------------------------------------------------------------!
    ! !-----------------------------------------------------------------------------------------------------!
    ! Allocation of variable sizes
    call allocatevars(an,as,ae,aw,ap,b,T,Told,Tn,resmat,il,ih,jl,jh)

    ! Initialising boundary conditions on temp array
    ! Doesnt at the moment with processors more than 1 because the partitioning is wrong
    call solutioninit(an,as,ae,aw,ap,b,T,il,ih,jl,jh)
    ! write(*,1600) T
   
    rc = 1

    ! !-----------------------------------------------------------------------------------------------------!
    ! !-----------------------------------------------------------------------------------------------------!
    ! ! Computation

    ! Solution loop variable initialisation
    ! Inititialising old temperature array
    Told(:,:) = 0
    Tn(:,:) = 0

    ! Initialising time counter
    time = 0
    iter = 0

    ! Choosing the solver --> jac, redblack, conj
    solvertype = 'jac'

    SELECT CASE (solvertype)

        ! jacobi solver
        CASE ("jac")
            ! Begin the solution loop
            do while ((time<t_final).and.(rc>res_max))
            ! do while (time<50)

                ! Calculation of solution using only jacobi solver
                call jac(an,as,ae,aw,ap,b,T,il,ih,jl,jh)

                ! COMMUNICATION - COPY FROM CLARA's PARTITIONING CODE TEST (need to update correct indices in comms)
                ! write(*,*) 'new loop'
                ! write(*,1600) T
                ! write(*,*) time, iter

                ! Non-blocking processor communication
                if (pid.NE.nprocs-1) then ! not on far right
                    call MPI_IRECV(T(:,ih+1),ny,MPI_DOUBLE_PRECISION,pid+1,tag,MPI_COMM_WORLD,req1,ierr) ! receive from the right
                end if
                if (pid.NE.0) then ! not on far left
                    call MPI_ISEND(T(:,il),ny,MPI_DOUBLE_PRECISION,pid-1,tag,MPI_COMM_WORLD,req2,ierr) ! send to the left
                end if
                if (pid.NE.0) then ! not on far left
                    call MPI_IRECV(T(:,il-1),ny,MPI_DOUBLE_PRECISION,pid-1,tag,MPI_COMM_WORLD,req1,ierr) ! receive from the left
                end if
                if (pid.NE.nprocs-1) then ! not on far right
                    call MPI_ISEND(T(:,ih),ny,MPI_DOUBLE_PRECISION,pid+1,tag,MPI_COMM_WORLD,req2,ierr)	! send to the right
                end if

                ! Wait until all processors have finished communicating
                call MPI_BARRIER(MPI_COMM_WORLD,ierr)
                    

                ! computing residuals
                call respar(aw,ae,an,as,ap,b,T,resil,resih,resjl,resjh,resmat)

                ! Calculate Domain averaged residual for stopping critterion
                rcurrent = SUM(SUM(ABS(resmat(resil:resih,resjl:resjh)),1),1) / ((resih-resil+1)*(resjh-resjl+1))

                ! Combining all residuals on each processor and sending to processor 0
                call MPI_REDUCE(rcurrent,rc,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

                if (pid.eq.0) then
                    rc = rc/nprocs
                end if

                call MPI_BCAST(rc,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

                ! write(*,*) rc


                ! Printing to screen after a certain number of iterations
                if (mod(iter,100).eq.0) then
					write(*,*) '-----------------------------------------'
					write(*,*) 'Time =', time
					write(*,*) 'Iteration =', iter
					write(*,*) 'Residual =', rc
					write(*,*) '-----------------------------------------'
                    write(*,1600) T
                end if


                ! updating counter
                time = time + dt
                iter = iter + 1

            end do

        
        CASE ("redblack")

            ! Begin the solution loop
            ! do while ((t<ttot).and.(r>rmax))
            do while (time<t_final)

                ! calculation of red nodes
                call rednodes(an,as,ae,aw,ap,b,T,il,ih,jl,jh)

                ! COMMUNICATION of red nodes

                ! calculation of black nodes
                call blacknodes(an,as,ae,aw,ap,b,T,il,ih,jl,jh)

                ! COMMUNICATION of black nodes
                    

                ! computing residuals
                call respar(aw,ae,an,as,ap,b,T,il,ih,jl,jh,resmat)

                ! Calculate Domain averaged residual for stopping critterion
                rcurrent = SUM(SUM(ABS(resmat(il:ih,jl:jh)),1),1) / ((nx-2)*(ny-2))

                ! Combining all residuals on each processor and sending to processor 0
                call MPI_REDUCE(rcurrent,rc,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

                if (pid == 0) then
                    rc = rc/nprocs
                end if

                call MPI_BCAST(rc,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)


                ! Printing to screen after a certain number of iterations
                if (mod(iter,100).eq.0) then
                    write(*,*) '      iter', '      res'
                    write(*,*) iter, rc
                end if

                ! updating old and new temperature values
                Told = T
                T = Tn

                ! updating counter
                time = time + dt
                iter = iter + 1

            end do

        CASE("Conj")

            ! ! Calculation of solution using Conjugate Gradient Method
            ! ! call CGSolve(an,as,ae,aw,ap,b,T)

    	CASE DEFAULT 
		WRITE(*,*) "No topology selected or incorrect selection"
		STOP

    END SELECT

    write(*,*) 'Output after solver'
    ! write(*,1600) T
	write(*,*) iter, time, rc
	write(*,*) pid,il,ih,jl,jh
	write(*,*) pid, resil,resih,resjl,resjh

    ! write(*,*) 'residual matrix'
    ! write(*,1600) resmat

    ! ! !-----------------------------------------------------------------------------------------------------!
    ! ! !-----------------------------------------------------------------------------------------------------!
    ! !! Getting final total temp array - SHOULD ALSO BE COVERED BY CLARA
    
    ! Gathering all temperature areas to processor 0
    allocate(Ttemp(il:ih,jl:jh))
    if (pid==0) then
        allocate(Ttot(nx,ny))
        Ttot(:,:) = 0
    end if
    ! allocate(Ttot(1,4))

    ! write(*,*) size(Ttemp), (ih-il+1)*(jh-jl+1), size(Ttot)
    Ttemp = T
    numcount = real((resih-resil+1)*(resjh-resjl+1)) ! Cannot get this working for an arbitrary number of procs atm
    ! numcount = 500*250
    write(*,*) numcount
    call MPI_GATHER(Ttemp,numcount,MPI_DOUBLE_PRECISION,Ttot,numcount, &
                    MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    ! call MPI_GATHER(Ttemp,2,MPI_DOUBLE_PRECISION,Ttot,2, &
    !                 MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

    ! if (pid == 0) then
    !     write(*,*) 'Output after gathering data'
    !     write(*,1600) Ttot
    ! end if

    ! write(*,*) 'im here'
    ! write(*,*) size(Ttot)
    

    if (pid.eq.0) then

        allocate(x(ny))
        allocate(y(nx))
        
        ! x is in the j direction, y is in the i direction
        do i = 1,ny
            y(i) = (i-1)*dy
        end do
        do j = 1,nx
            x(j) = (j-1)*dx
        end do

        ! Writing updated initial distribution to file
        write(file_name, "(A14)") "Tecplotmax.tec"
        call tec_2D ( iunit, nx, ny, x, y, Ttot, file_name )
    end if
    
    ! NOTE: Need to update this to match the number of spatial divisions (nx)
    1100 FORMAT(8(F20.10,1x))
    1600 FORMAT(8(F10.6,1x))
    1400 FORMAT(8(F14.8,1x))
    1200 FORMAT(I2.1,I2.1,6(F8.4,1x))
    1300 FORMAT(I2.1,6(F10.8,1x))


    CALL MPI_FINALIZE(ierr)

end program test



    ! SUBROUTINES
!----------------------------------------------------------------------------- 
subroutine tec_2D ( iunit, nx, ny, x, y, T, file_name )

  IMPLICIT NONE

  Integer ( kind = 4 ) iunit,nx,ny,i,j,ierr
  Real ( kind = 8 ) T(nx,ny)
  Real ( kind = 8 ) x(ny)
  Real ( kind = 8 ) y(nx)
  Character(len=1024) :: file_name

!   Character(80), parameter ::  file_name = 'TecPlot2Dnew.tec'

  open ( unit = iunit, file = file_name, form = 'formatted', access = 'sequential', status = 'replace', iostat = ierr )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) '  Error opening file : tecplot_2D '
    stop
  end if
   
  write ( iunit, '(a)' ) 'Title="' // trim ( 'Temperature Data' ) // '"'
  write ( iunit, '(a)' ) 'Variables=' // trim ( '"X","Y","T"' )

  write ( iunit, '(a)' ) ' '
  write ( iunit, '(a,i6,a,i6,a)' ) 'Zone I=', ny, ', J=', nx, ', F=POINT'
 
  do j = 1, nx
    do i = 1, ny
      write ( iunit, '(2f10.3,g15.6)' ) y(i), x(j), T(i,j)
    end do
  end do
  
  close ( unit = iunit )

end subroutine tec_2D