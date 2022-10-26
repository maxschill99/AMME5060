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
!|	! Other options: "cart" and "slabs" 							 			 |
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
		
			CALL graph_partition()
			
			! Not reordering because in this configuration MPI doesn't appear to want to reorder them anyway.
			! (I ran the code with reorder=true and compared before/after for different numbers of processors)
			! So setting false to ensure future code still works on the /off/ chance it did reorder them.
			CALL MPI_GRAPH_CREATE(MPI_COMM_WORLD, Nprocs, indexes, edges, .false., COMM_TOPO, ierr)

			! MPI_Graph_neighbors_count retrieves the number of neighbours for a given rank in the communicator
			CALL MPI_GRAPH_NEIGHBORS_COUNT(COMM_TOPO, pid, neighbours_count, ierr)
			WRITE(*,*) 'My pid is', pid, '. I have', neighbours_count, 'neighbours.'
			! MPI_Graph_neighbors - Returns the neighbors of a node associated with a graph topology.
			
			ALLOCATE( neighbours_array(neighbours_count) )
			CALL MPI_GRAPH_NEIGHBORS(COMM_TOPO, pid, neighbours_count, neighbours_array, ierr)
			write(*,*) 'i am pid', pid, 'my neighbours are', neighbours_array
			
		CASE ("cart")
		
			CALL cart_partition()
		
			! Create a new cartesian communicator based on the above analysis
			CALL MPI_CART_CREATE(MPI_COMM_WORLD, ndims, dims, periods, .true., COMM_TOPO, ierr)
			CALL MPI_COMM_RANK(COMM_TOPO, pid, ierr) ! New processor ID number
			CALL MPI_CART_COORDS(COMM_TOPO, pid, ndims, coords, ierr)
			! coords is an array with pid location in (x, y) but y counts from top down and counts start from zero
			! Send/Recv id's for each processor horizontally and vertically (if applicable)
			CALL MPI_CART_SHIFT(COMM_TOPO, 0, 1, west1, east1, ierr)
			CALL MPI_CART_SHIFT(COMM_TOPO, 1, 1, north, south, ierr)
			
			east2 = MPI_PROC_NULL
			west2 = MPI_PROC_NULL
			
			! Nodes and indices in x
			CALL get_nodes(nx, dims(1), coords(1), node_low_x, node_high_x)
			CALL get_indices(dims(1), coords(1), node_low_x, node_high_x, ind_low_x, ind_high_x, ncalcpoints_x)
			! Nodes and indices in y
			CALL get_nodes(ny, dims(2), coords(2), node_low_y, node_high_y)
			CALL get_indices(dims(2), coords(2), node_low_y, node_high_y, ind_low_y, ind_high_y, ncalcpoints_y)
			
			! Sends/Recvs that happen			
			ind_low_east1 		= ind_low_y+1
			ind_high_east1 		= ind_high_y-1
			ncalcpoints_y_east1 = ncalcpoints_y
			ind_low_west1		= ind_low_y+1
			ind_high_west1		= ind_high_y-1
			ncalcpoints_y_west1 = ncalcpoints_y
			! Sends/Recvs that go to MPI_PROC_NULL since not graph topology
			ind_low_east2 		= ind_low_east1
			ind_high_east2 		= ind_high_east1
			ncalcpoints_y_east2 = ncalcpoints_y_east1	
			ind_low_west2		= ind_low_west1
			ind_high_west2		= ind_high_west1
			ncalcpoints_y_west2 = ncalcpoints_y_west1
			
	
		CASE ("slabs")
		
			CALL slab_partition()
		
			! Create a new cartesian communicator based on the above analysis
			CALL MPI_CART_CREATE(MPI_COMM_WORLD, ndims, dims, periods, .true., COMM_TOPO, ierr)
			CALL MPI_COMM_RANK(COMM_TOPO, pid, ierr) ! New processor ID number
			CALL MPI_CART_COORDS(COMM_TOPO, pid, ndims, coords, ierr)
			! coords is an array with pid location in (x, y) but y counts from top down and counts start from zero
			! Send/Recv id's for each processor horizontally and vertically (if applicable)
			CALL MPI_CART_SHIFT(COMM_TOPO, 0, 1, west1, east1, ierr)
			CALL MPI_CART_SHIFT(COMM_TOPO, 1, 1, north, south, ierr)
			
			east2 = MPI_PROC_NULL
			west2 = MPI_PROC_NULL
			
			! Nodes and indices in x
			CALL get_nodes(nx, dims(1), coords(1), node_low_x, node_high_x)
			CALL get_indices(dims(1), coords(1), node_low_x, node_high_x, ind_low_x, ind_high_x, ncalcpoints_x)
			! Nodes and indices in y
			CALL get_nodes(ny, dims(2), coords(2), node_low_y, node_high_y)
			CALL get_indices(dims(2), coords(2), node_low_y, node_high_y, ind_low_y, ind_high_y, ncalcpoints_y)
			
			! Sends/Recvs that happen			
			ind_low_east1 		= ind_low_y+1
			ind_high_east1 		= ind_high_y-1
			ncalcpoints_y_east1 = ncalcpoints_y
			ind_low_west1		= ind_low_y+1
			ind_high_west1		= ind_high_y-1
			ncalcpoints_y_west1 = ncalcpoints_y
			! Sends/Recvs that go to MPI_PROC_NULL since not graph topology
			ind_low_east2 		= ind_low_east1
			ind_high_east2 		= ind_high_east1
			ncalcpoints_y_east2 = ncalcpoints_y_east1	
			ind_low_west2		= ind_low_west1
			ind_high_west2		= ind_high_west1
			ncalcpoints_y_west2 = ncalcpoints_y_west1
	
		CASE DEFAULT 
		  WRITE(*,*) "No topology selected or incorrect selection"
		  STOP
	END SELECT
	
	! New vector type for sending and receiving the rows north and south (data not contiguous along a row)
	CALL MPI_TYPE_VECTOR(ncalcpoints_x, 1, ncalcpoints_y+2, MPI_DOUBLE_PRECISION, NS_ROW_SENDRECV , ierr)
	CALL MPI_TYPE_COMMIT(NS_ROW_SENDRECV, ierr)
	
	
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


	! ALLOCATE( status_array(MPI_STATUS_SIZE, 12) )
	! ALLOCATE( request_array(12) )
	
	! ! ------- RECEIVE FROM THE RIGHT/EAST, SEND TO THE RIGHT/EAST
	! ! Receive from the east1, and put it into the space for the ghost node on the east side
	! CALL MPI_IRECV( Tk_1(ind_low_east1:ind_high_east1, ind_high_x), ncalcpoints_y_east1, MPI_DOUBLE_PRECISION, &
					! east1, tag, COMM_TOPO, request_array(1), ierr)
	! ! Send to the east1
	! CALL MPI_ISEND( Tk_1(ind_low_east1:ind_high_east1, ind_high_x-1), ncalcpoints_y_east1, MPI_DOUBLE_PRECISION, &
					! east1, tag, COMM_TOPO, request_array(2), ierr)
					
	! ! Receive from the east2, and put it into the space for the ghost node on the east side
	! CALL MPI_IRECV( Tk_1(ind_low_east2:ind_high_east2, ind_high_x), ncalcpoints_y_east2, MPI_DOUBLE_PRECISION, &
					! east2, tag, COMM_TOPO, request_array(3), ierr)
	! ! Send to the east2
	! CALL MPI_ISEND( Tk_1(ind_low_east2:ind_high_east2, ind_high_x-1), ncalcpoints_y_east2, MPI_DOUBLE_PRECISION, &
					! east2, tag, COMM_TOPO, request_array(4), ierr)


	! ! ------- RECEIVE FROM THE LEFT/WEST, SEND TO THE LEFT/WEST
	! ! Receive from the west1, put it into the space for the ghost node on the west side 
	! CALL MPI_IRECV( Tk_1(ind_low_west1:ind_high_west1, ind_low_x), ncalcpoints_y_west1, MPI_DOUBLE_PRECISION, & 
					! west1, tag, COMM_TOPO, request_array(5), ierr)
	! ! Send to the west1
	! CALL MPI_ISEND( Tk_1(ind_low_west1:ind_high_west1, ind_low_x+1), ncalcpoints_y_west1, MPI_DOUBLE_PRECISION, &
					! west1, tag, COMM_TOPO, request_array(6), ierr)
					
	! ! Receive from the west2, put it into the space for the ghost node on the west side 
	! CALL MPI_IRECV( Tk_1(ind_low_west2:ind_high_west2, ind_low_x), ncalcpoints_y_west2, MPI_DOUBLE_PRECISION, & 
					! west2, tag, COMM_TOPO, request_array(7), ierr)
	! ! Send to the west2
	! CALL MPI_ISEND( Tk_1(ind_low_west2:ind_high_west2, ind_low_x+1), ncalcpoints_y_west2, MPI_DOUBLE_PRECISION, &
					! west2, tag, COMM_TOPO, request_array(8), ierr)
					
	! ! ------- RECEIVE FROM THE NORTH, SEND TO THE NORTH
	! ! Receive from the north, put it into the space for the ghost node 
	! CALL MPI_IRECV( Tk_1(ind_low_y, ind_low_x+1), 1, NS_ROW_SENDRECV, &
					! north, tag, COMM_TOPO, request_array(9), ierr)
	! ! Send to the north
	! CALL MPI_ISEND( Tk_1(ind_low_y+1, ind_low_x+1), 1, NS_ROW_SENDRECV, &
					! north, tag, COMM_TOPO, request_array(10), ierr)
	! ! ------- RECEIVE FROM THE SOUTH, SEND TO THE SOUTH
	! ! Receive from the south, put it into the space for the ghost node 
	! CALL MPI_IRECV( Tk_1(ind_high_y, ind_low_x+1), 1, NS_ROW_SENDRECV, &
					! south, tag, COMM_TOPO, request_array(11), ierr)
	! ! Send to the south
	! CALL MPI_ISEND( Tk_1(ind_high_y-1, ind_low_x+1), 1, NS_ROW_SENDRECV, &
					! south, tag, COMM_TOPO, request_array(12), ierr)

	! ! Wait for data sends to complete before black points start referencing red points
	! CALL MPI_WAITALL(12, request_array, status_array, ierr)
