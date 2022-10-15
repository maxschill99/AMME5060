! MAJOR PROJECT MAIN CODE OUTLINE
PROGRAM MAIN
! CALLING MODULES
	USE variablemodule 		! Contains variable allocation all problem-specific variable values
	USE nodemodule			! Obtains nodes and indices of the domain considered by each processor
	USE partitionmodule

	include 'mpif.h'

! INITIALISE VARIABLES 
	! Variable declaration and assignment
	CALL intialise()

! -------------------------------------------------------------------------------	
!| CODE CASES																	 |
	topology = "graph" 	!														 |
!|	! &&& "cart" "none" "slabs" could also be other cases we could select for    |
! -------------------------------------------------------------------------------

! 	INITIALISE MPI
	CALL MPI_INIT(ierr)
	CALL MPI_COMM_RANK(MPI_COMM_WORLD, pid, ierr) ! Getting processor ID number
	CALL MPI_COMM_SIZE(MPI_COMM_WORLD, Nprocs, ierr) ! Getting number of total processors in global communicator

! INITIALISE PARTITIONING
	! module
	! outputs - local x, local y, low + high indxx, low + high indxy, low + high nodex, low + high nodey, all neighbour pid vals, idx east1/east2/west1/west2, number of send/recv points toe neighbours (east1/east2/west1/west2) 
	
	SELECT CASE (topology)
	
		CASE ("graph")
		
			graph_partition()
			
			! Not reordering because in this configuration MPI doesn't appear to want to reorder them anyway.
			! (I ran the code with reorder=true and compared before/after for different numbers of processors)
			! So setting false to ensure future code still works on the /off/ chance it did reorder them.
			CALL MPI_GRAPH_CREATE(MPI_COMM_WORLD, Nprocs, indexes, edges, .false., COMM_GRAPH, ierr)

			! MPI_Graph_neighbors_count retrieves the number of neighbours for a given rank in the communicator
			CALL MPI_GRAPH_NEIGHBORS_COUNT(COMM_GRAPH, pid, neighbours_count, ierr)
			WRITE(*,*) 'My pid is', pid, '. I have', neighbours_count, 'neighbours.'
			! MPI_Graph_neighbors - Returns the neighbors of a node associated with a graph topology.
			
			ALLOCATE( neighbours_array(neighbours_count) )
			CALL MPI_GRAPH_NEIGHBORS(COMM_GRAPH, pid, neighbours_count, neighbours_array, ierr)
			write(*,*) 'i am pid', pid, 'my neighbours are', neighbours_array
			
		CASE DEFAULT 
		  WRITE(*,*) "No topology selected or incorrect selection"
		  STOP
	END SELECT
	
	
! CREATING COMMUNICATOR

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
