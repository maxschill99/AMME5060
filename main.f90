! MAJOR PROJECT MAIN CODE OUTLINE
PROGRAM MAIN
! CALLING MODULES


include 'mpif.h'


! INITIALISE VARIABLES 
	! Iniitialise variables
	! variable module - grid size, grid spacing

! INITIALISE PARTITIONING
	! module
	! outputs - local x, local y, low + high indxx, low + high indxy, low + high nodex, low + high nodey, all neighbour pid vals, idx east1/east2/west1/west2, number of send/recv points toe neighbours (east1/east2/west1/west2) 
	

! INITIALISE TEMP DISTRIBUTION 
	! boundary conditions most pizza like


! SOLVER
	! outer loop: time stepping
		! solve for time step n+1 and while r<err
			! multigrid module with jacobi/gauss seidel 
				! inverted V
				! number of levels dictated by grid points and divisible factor used 
			! this will involve looping through temperature arrays
		
		! end solve for timestep n+1
		
		! If statement check that t=somevalue
		! WRITE TO ONE FILE 
		
	! end time stepping
	
	! Residual module
	! Solver module - jacobi/gauss seidel
	! Multigrid module
	! Call communication
	
	





END PROGRAM MAIN
