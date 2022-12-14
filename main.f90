! MAJOR PROJECT MAIN CODE OUTLINE
PROGRAM MAIN
! CALLING MODULES
	USE variablemodule 		! Contains variable allocation all problem-specific variable values
	USE nodemodule			! Obtains nodes and indices of the domain considered by each processor
	USE partitionmodule		! Contains subroutines to partition domain and obtain domain info for each processor
	USE outputmodule		! Creates tec file
	USE jacobi				! Solves 2d diffusion using jacobi
	USE residuals			! Computes residuals
	USE solinit				! Initialises the solution

	include 'mpif.h'

! INITIALISE VARIABLES 
	! Variable declaration and assignment
	CALL intialise()

! ----------------------------------------------	
!| CODE CASES									|
	topology = "slabs" 	!						|
!|		Options: "graph" "cart" "slabs" 		|
    solvertype = "jac"	!						|
!| 		Options: "jac", "redblack", "conj"		|
! ----------------------------------------------

! 	INITIALISE MPI
	CALL MPI_INIT(ierr)
	CALL MPI_COMM_RANK(MPI_COMM_WORLD, pid, ierr) ! Getting processor ID number
	CALL MPI_COMM_SIZE(MPI_COMM_WORLD, Nprocs, ierr) ! Getting number of total processors in global communicator

! Begin timer
t1 = MPI_WTIME()

	! -----------------------------------------------------------------------------------------------------
	! TOPOLOGY
	! -----------------------------------------------------------------------------------------------------

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
		
			!  Obtain how many processors in each dimension (set up `dims` and `periods`)
			CALL cart_partition()
		
			! Create a new cartesian communicator based on the above analysis
			CALL MPI_CART_CREATE(MPI_COMM_WORLD, ndims, dims, periods, .true., COMM_TOPO, ierr)
			CALL MPI_COMM_RANK(COMM_TOPO, pid, ierr) ! New processor ID number
			CALL MPI_CART_COORDS(COMM_TOPO, pid, ndims, coords, ierr)
			! coords is an array with pid location in (x, y) but y counts from top down and counts start from zero
			! Send/Recv id's for each processor horizontally and vertically (if applicable)
			CALL MPI_CART_SHIFT(COMM_TOPO, 0, 1, west1, east1, ierr)
			CALL MPI_CART_SHIFT(COMM_TOPO, 1, 1, north, south, ierr)
			
			! Define indices of arrays and nodes considered, and data associated with east2/west2
			CALL cart_nodesindices()
			
		CASE ("slabs")
		
			! Set up `dims` and `periods`
			CALL slab_partition()
		
			! Create a new cartesian communicator based on the above analysis
			CALL MPI_CART_CREATE(MPI_COMM_WORLD, ndims, dims, periods, .true., COMM_TOPO, ierr)
			CALL MPI_COMM_RANK(COMM_TOPO, pid, ierr) ! New processor ID number
			CALL MPI_CART_COORDS(COMM_TOPO, pid, ndims, coords, ierr)
			! coords is an array with pid location in (x, y) but y counts from top down and counts start from zero
			! Send/Recv id's for each processor horizontally and vertically (if applicable)
			CALL MPI_CART_SHIFT(COMM_TOPO, 0, 1, west1, east1, ierr)
			CALL MPI_CART_SHIFT(COMM_TOPO, 1, 1, north, south, ierr)
			
			! Define indices of arrays and nodes considered, and data associated with east2/west2
			CALL slab_nodesindices()
	
		CASE DEFAULT 
		  WRITE(*,*) "No topology selected or incorrect selection"
		  STOP
	END SELECT
	
	! New vector type for sending and receiving the rows north and south (data not contiguous along a row)
	CALL MPI_TYPE_VECTOR(ncalcpoints_x, 1, ncalcpoints_y+2, MPI_DOUBLE_PRECISION, NS_ROW_SENDRECV , ierr)
	CALL MPI_TYPE_COMMIT(NS_ROW_SENDRECV, ierr)
	
	
	! -----------------------------------------------------------------------------------------------------
	! INITIALISATION
	! -----------------------------------------------------------------------------------------------------

	! Indices for computation
	il = ind_low_y; ih = ind_high_y; jl = ind_low_x; jh = ind_high_x ! Indices for temperature calculations
	resil = node_low_y; resih = node_high_y; resjl = node_low_x; resjh = node_high_x ! Nodes for residual calculations

	! INITIALISE TEMP DISTRIBUTION 
	! Allocation of variable sizes
	call allocatevars(an,as,ae,aw,ap,b,T,Told,Tn,res,il,ih,jl,jh)

	! Initialising boundary conditions on temp array
	call solutioninit(an,as,ae,aw,ap,b,T,il,ih,jl,jh)

	! Initialising new temp array and setting its value to T
	Tn(:,:) = 0
	Tn = T
	Told(:,:) = 0

	! INITIAL RESIDUAL
	! Initialising residual matrix
	allocate(resmat(resil:resih,resjl:resjh))
	resmat(:,:) = 0.0
	rc = 1

	! SOLUTION ACCURACY CHECK
	! CFL = ((1/(dx**2) + 1/(dy**2))*alpha*dt)
	CFL = (1 - (4*dt*alpha)/dx**2)
	uncondstab = (dt*alpha)/dx**2
	if (pid.eq.0) then
		write(*,*) 'Explicit Stability: ', CFL, 'Implicit Accuracy: ', uncondstab
	end if
	if (CFL.GT.0.5) then
		write(*,*) 'Explicit Stability not met, needs to be less than 0.5'
		STOP
	elseif (CFL.LT.0) then
		write(*,*) 'Explicit Stability less than 0'
		STOP
	end if

	!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Initialising time counter and iteration counter
	time = 0
	iter = 0

	! -----------------------------------------------------------------------------------------------------
	! SEND/RECV SET UP
	! -----------------------------------------------------------------------------------------------------

	ALLOCATE( status_array(MPI_STATUS_SIZE, 12) )
	ALLOCATE( request_array(12) )
	ALLOCATE( status_array_gather(MPI_STATUS_SIZE, Nprocs+1) )
	ALLOCATE( request_array_gather(Nprocs+1) )

	! All processors first create subarrays. These are the final temperature array cleansed of the extra ghost nodes
	subarray_Nrows = node_high_y-node_low_y + 1
	subarray_Ncols = node_high_x-node_low_x + 1
	
	CALL MPI_Type_create_subarray(2, [ncalcpoints_y+2,ncalcpoints_x+2], [subarray_Nrows,subarray_Ncols], &
	[node_low_y-ind_low_y, node_low_x-ind_low_x], MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, SENDINGSUBBARAY, ierr)
	CALL MPI_TYPE_COMMIT(SENDINGSUBBARAY, ierr)
	
	! If pid 0, then make space to obtain relevant information about each subarray, and the final temp array for the domain 
	IF (pid .EQ. 0) THEN
		ALLOCATE( Tfinal(1:ny,1:nx) )
		ALLOCATE( subarray_rows_array(1:Nprocs) )
		ALLOCATE( subarray_cols_array(1:Nprocs) )
		ALLOCATE( subarray_row_start(1:Nprocs) )
		ALLOCATE( subarray_col_start(1:Nprocs) )
	END IF
	
	! ---- RECIEVE SUBARRAY TYPE ----
	! Gather subarray information from all processors
		! Subarray row sizes
		CALL MPI_GATHER( subarray_Nrows, 1, MPI_INTEGER, subarray_rows_array, 1, MPI_INTEGER, 0, COMM_TOPO, ierr )
		! Subarray column sizes
		CALL MPI_GATHER( subarray_Ncols, 1, MPI_INTEGER, subarray_cols_array, 1, MPI_INTEGER, 0, COMM_TOPO, ierr )
		! Start row location in the final matrix (minus 1 because start location starts at zero)
		CALL MPI_GATHER( node_low_y-1, 	 1, MPI_INTEGER, subarray_row_start,  1, MPI_INTEGER, 0, COMM_TOPO, ierr )
		! Start column location in the final matrix (minus 1 because start location starts at zero)
		CALL MPI_GATHER( node_low_x-1, 	 1, MPI_INTEGER, subarray_col_start,  1, MPI_INTEGER, 0, COMM_TOPO, ierr )


	! -----------------------------------------------------------------------------------------------------
	! SOLVING
	! -----------------------------------------------------------------------------------------------------
    SELECT CASE (solvertype)

        ! JACOBI
        CASE ("jac")
            ! Begin the solution loop
			! do while (rc>res_max)
			do while ((rc>res_max).and.(time<t_final))
				!-------------------------------------------------------------------!
				Told = T

				! Calculation of solution using only jacobi solver				
				call jac(an,as,ae,aw,ap,b,T,Told,il,ih,jl,jh)

				!-------------------------------------------------------------------!
				! COMMUNICATION				
				! ------- RECEIVE FROM THE RIGHT/EAST, SEND TO THE RIGHT/EAST
				! Receive from the east1, and put it into the space for the ghost node on the east side
				CALL MPI_IRECV( T(ind_low_east1:ind_high_east1, ind_high_x), ncalcpoints_y_east1, MPI_DOUBLE_PRECISION, &
								east1, tag, COMM_TOPO, request_array(1), ierr)
				! Send to the east1
				CALL MPI_ISEND( T(ind_low_east1:ind_high_east1, ind_high_x-1), ncalcpoints_y_east1, MPI_DOUBLE_PRECISION, &
								east1, tag, COMM_TOPO, request_array(2), ierr)
								
				! Receive from the east2, and put it into the space for the ghost node on the east side
				CALL MPI_IRECV( T(ind_low_east2:ind_high_east2, ind_high_x), ncalcpoints_y_east2, MPI_DOUBLE_PRECISION, &
								east2, tag, COMM_TOPO, request_array(3), ierr)
				! Send to the east2
				CALL MPI_ISEND( T(ind_low_east2:ind_high_east2, ind_high_x-1), ncalcpoints_y_east2, MPI_DOUBLE_PRECISION, &
								east2, tag, COMM_TOPO, request_array(4), ierr)


				! ------- RECEIVE FROM THE LEFT/WEST, SEND TO THE LEFT/WEST
				! Receive from the west1, put it into the space for the ghost node on the west side 
				CALL MPI_IRECV( T(ind_low_west1:ind_high_west1, ind_low_x), ncalcpoints_y_west1, MPI_DOUBLE_PRECISION, & 
								west1, tag, COMM_TOPO, request_array(5), ierr)
				! Send to the west1
				CALL MPI_ISEND( T(ind_low_west1:ind_high_west1, ind_low_x+1), ncalcpoints_y_west1, MPI_DOUBLE_PRECISION, &
								west1, tag, COMM_TOPO, request_array(6), ierr)
								
				! Receive from the west2, put it into the space for the ghost node on the west side 
				CALL MPI_IRECV( T(ind_low_west2:ind_high_west2, ind_low_x), ncalcpoints_y_west2, MPI_DOUBLE_PRECISION, & 
								west2, tag, COMM_TOPO, request_array(7), ierr)
				! Send to the west2
				CALL MPI_ISEND( T(ind_low_west2:ind_high_west2, ind_low_x+1), ncalcpoints_y_west2, MPI_DOUBLE_PRECISION, &
								west2, tag, COMM_TOPO, request_array(8), ierr)
								
				! ------- RECEIVE FROM THE NORTH, SEND TO THE NORTH
				! Receive from the north, put it into the space for the ghost node 
				CALL MPI_IRECV( T(ind_low_y, ind_low_x+1), 1, NS_ROW_SENDRECV, &
								north, tag, COMM_TOPO, request_array(9), ierr)
				! Send to the north
				CALL MPI_ISEND( T(ind_low_y+1, ind_low_x+1), 1, NS_ROW_SENDRECV, &
								north, tag, COMM_TOPO, request_array(10), ierr)
				! ------- RECEIVE FROM THE SOUTH, SEND TO THE SOUTH
				! Receive from the south, put it into the space for the ghost node 
				CALL MPI_IRECV( T(ind_high_y, ind_low_x+1), 1, NS_ROW_SENDRECV, &
								south, tag, COMM_TOPO, request_array(11), ierr)
				! Send to the south
				CALL MPI_ISEND( T(ind_high_y-1, ind_low_x+1), 1, NS_ROW_SENDRECV, &
								south, tag, COMM_TOPO, request_array(12), ierr)

				! Wait for data sends to complete before black points start referencing red points
				CALL MPI_WAITALL(12, request_array, status_array, ierr)

				!-------------------------------------------------------------------!
				! RESIDUALS
				! computing residuals
				rcurrent = SUM(sum(ABS(Told(resil:resih,resjl:resjh) - T(resil:resih,resjl:resjh)),1),1) &
									/ ((resih-resil+1)*(resjh-resjl+1))

				! Summing processor residuals to get global resiudal and broadcasting to get average
				call MPI_ALLREDUCE(rcurrent,rc,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_TOPO,ierr)
				rc = rc/nprocs

				!-------------------------------------------------------------------!
				! Printing to screen after a certain amount of time
				if ((MOD(iter,1000).eq.0)) then

					! TECPLOT

					! These subarrays are now sent to PID 0
					CALL MPI_ISEND( T, 1, SENDINGSUBBARAY, 0, tag2, COMM_TOPO, request_array_gather(1), ierr)	
							
					! Pid 0 receiving data and putting into final file
					IF (pid .EQ. 0) THEN

						DO i = 0, Nprocs-1
						
							! Creating receive subarray type bespoke to each processor
							CALL MPI_Type_create_subarray(2, [ny, nx], [subarray_rows_array(i+1), subarray_cols_array(i+1)], &
							[subarray_row_start(i+1),subarray_col_start(i+1)], MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, RECVSUBBARAY, ierr)
							CALL MPI_TYPE_COMMIT(RECVSUBBARAY, ierr)
						
							! Receiving with this new receiving subarray type
							CALL MPI_IRECV( Tfinal, 1, RECVSUBBARAY, &
								i, tag2, COMM_TOPO, request_array_gather(i+2), ierr)
						END DO

						CALL MPI_WAITALL(Nprocs+1, request_array_gather, status_array_gather, ierr)

						! --- PUTTING INTO FILE ---
						! x vector of whole domain 
						xtot 	= 0. + dx * [(i, i=0,(nx-1))] ! Implied DO loop used
						! y vector of whole domain 
						ytot 	= 0. + dy * [(i, i=0,(ny-1))] ! Implied DO loop used

						! Writing updated solution to file at each new iteration
						write(file_name, "(A9,I5,A4)") "TecPlot2D",int(time),".tec"
						CALL tecplot_2D ( iunit, nx, ny, xtot, ytot, Tfinal,  file_name )

						write(*,*) '-----------------------------------------'
						write(*,*) 'Time =', time
						write(*,*) 'Iteration =', iter
						write(*,*) 'Residual =', rc
						write(*,*) 'Average Temperature =', sum(sum((Tfinal),1),1)/(nx*ny)
						write(*,*) '-----------------------------------------'
						
					END IF
					
				end if

				!-------------------------------------------------------------------!
				! updating counter
				time = time + dt
				iter = iter + 1

			end do

		! REDBLACK
        CASE ("redblack")

            ! Begin the solution loop
			do while ((rc>res_max).and.(time<t_final))
				!-------------------------------------------------------------------!
				Told = T

				! calculation of red nodes
				call rednodes(an,as,ae,aw,ap,b,T,Told,il,ih,jl,jh)

				!-------------------------------------------------------------------!
				! COMMUNICATION of red nodes			
				! ------- RECEIVE FROM THE RIGHT/EAST, SEND TO THE RIGHT/EAST
				! Receive from the east1, and put it into the space for the ghost node on the east side
				CALL MPI_IRECV( T(ind_low_east1:ind_high_east1, ind_high_x), ncalcpoints_y_east1, MPI_DOUBLE_PRECISION, &
								east1, tag, COMM_TOPO, request_array(1), ierr)
				! Send to the east1
				CALL MPI_ISEND( T(ind_low_east1:ind_high_east1, ind_high_x-1), ncalcpoints_y_east1, MPI_DOUBLE_PRECISION, &
								east1, tag, COMM_TOPO, request_array(2), ierr)
								
				! Receive from the east2, and put it into the space for the ghost node on the east side
				CALL MPI_IRECV( T(ind_low_east2:ind_high_east2, ind_high_x), ncalcpoints_y_east2, MPI_DOUBLE_PRECISION, &
								east2, tag, COMM_TOPO, request_array(3), ierr)
				! Send to the east2
				CALL MPI_ISEND( T(ind_low_east2:ind_high_east2, ind_high_x-1), ncalcpoints_y_east2, MPI_DOUBLE_PRECISION, &
								east2, tag, COMM_TOPO, request_array(4), ierr)


				! ------- RECEIVE FROM THE LEFT/WEST, SEND TO THE LEFT/WEST
				! Receive from the west1, put it into the space for the ghost node on the west side 
				CALL MPI_IRECV( T(ind_low_west1:ind_high_west1, ind_low_x), ncalcpoints_y_west1, MPI_DOUBLE_PRECISION, & 
								west1, tag, COMM_TOPO, request_array(5), ierr)
				! Send to the west1
				CALL MPI_ISEND( T(ind_low_west1:ind_high_west1, ind_low_x+1), ncalcpoints_y_west1, MPI_DOUBLE_PRECISION, &
								west1, tag, COMM_TOPO, request_array(6), ierr)
								
				! Receive from the west2, put it into the space for the ghost node on the west side 
				CALL MPI_IRECV( T(ind_low_west2:ind_high_west2, ind_low_x), ncalcpoints_y_west2, MPI_DOUBLE_PRECISION, & 
								west2, tag, COMM_TOPO, request_array(7), ierr)
				! Send to the west2
				CALL MPI_ISEND( T(ind_low_west2:ind_high_west2, ind_low_x+1), ncalcpoints_y_west2, MPI_DOUBLE_PRECISION, &
								west2, tag, COMM_TOPO, request_array(8), ierr)
								
				! ------- RECEIVE FROM THE NORTH, SEND TO THE NORTH
				! Receive from the north, put it into the space for the ghost node 
				CALL MPI_IRECV( T(ind_low_y, ind_low_x+1), 1, NS_ROW_SENDRECV, &
								north, tag, COMM_TOPO, request_array(9), ierr)
				! Send to the north
				CALL MPI_ISEND( T(ind_low_y+1, ind_low_x+1), 1, NS_ROW_SENDRECV, &
								north, tag, COMM_TOPO, request_array(10), ierr)
				! ------- RECEIVE FROM THE SOUTH, SEND TO THE SOUTH
				! Receive from the south, put it into the space for the ghost node 
				CALL MPI_IRECV( T(ind_high_y, ind_low_x+1), 1, NS_ROW_SENDRECV, &
								south, tag, COMM_TOPO, request_array(11), ierr)
				! Send to the south
				CALL MPI_ISEND( T(ind_high_y-1, ind_low_x+1), 1, NS_ROW_SENDRECV, &
								south, tag, COMM_TOPO, request_array(12), ierr)

				! Wait for data sends to complete before black points start referencing red points
				CALL MPI_WAITALL(12, request_array, status_array, ierr)



				Told = T
				!-------------------------------------------------------------------!
				! calculation of black nodes
				call blacknodes(an,as,ae,aw,ap,b,T,Told,il,ih,jl,jh)


				!-------------------------------------------------------------------!
				! COMMUNICATION of black nodes
				! ------- RECEIVE FROM THE RIGHT/EAST, SEND TO THE RIGHT/EAST
				! Receive from the east1, and put it into the space for the ghost node on the east side
				CALL MPI_IRECV( T(ind_low_east1:ind_high_east1, ind_high_x), ncalcpoints_y_east1, MPI_DOUBLE_PRECISION, &
								east1, tag, COMM_TOPO, request_array(1), ierr)
				! Send to the east1
				CALL MPI_ISEND( T(ind_low_east1:ind_high_east1, ind_high_x-1), ncalcpoints_y_east1, MPI_DOUBLE_PRECISION, &
								east1, tag, COMM_TOPO, request_array(2), ierr)
								
				! Receive from the east2, and put it into the space for the ghost node on the east side
				CALL MPI_IRECV( T(ind_low_east2:ind_high_east2, ind_high_x), ncalcpoints_y_east2, MPI_DOUBLE_PRECISION, &
								east2, tag, COMM_TOPO, request_array(3), ierr)
				! Send to the east2
				CALL MPI_ISEND( T(ind_low_east2:ind_high_east2, ind_high_x-1), ncalcpoints_y_east2, MPI_DOUBLE_PRECISION, &
								east2, tag, COMM_TOPO, request_array(4), ierr)


				! ------- RECEIVE FROM THE LEFT/WEST, SEND TO THE LEFT/WEST
				! Receive from the west1, put it into the space for the ghost node on the west side 
				CALL MPI_IRECV( T(ind_low_west1:ind_high_west1, ind_low_x), ncalcpoints_y_west1, MPI_DOUBLE_PRECISION, & 
								west1, tag, COMM_TOPO, request_array(5), ierr)
				! Send to the west1
				CALL MPI_ISEND( T(ind_low_west1:ind_high_west1, ind_low_x+1), ncalcpoints_y_west1, MPI_DOUBLE_PRECISION, &
								west1, tag, COMM_TOPO, request_array(6), ierr)
								
				! Receive from the west2, put it into the space for the ghost node on the west side 
				CALL MPI_IRECV( T(ind_low_west2:ind_high_west2, ind_low_x), ncalcpoints_y_west2, MPI_DOUBLE_PRECISION, & 
								west2, tag, COMM_TOPO, request_array(7), ierr)
				! Send to the west2
				CALL MPI_ISEND( T(ind_low_west2:ind_high_west2, ind_low_x+1), ncalcpoints_y_west2, MPI_DOUBLE_PRECISION, &
								west2, tag, COMM_TOPO, request_array(8), ierr)
								
				! ------- RECEIVE FROM THE NORTH, SEND TO THE NORTH
				! Receive from the north, put it into the space for the ghost node 
				CALL MPI_IRECV( T(ind_low_y, ind_low_x+1), 1, NS_ROW_SENDRECV, &
								north, tag, COMM_TOPO, request_array(9), ierr)
				! Send to the north
				CALL MPI_ISEND( T(ind_low_y+1, ind_low_x+1), 1, NS_ROW_SENDRECV, &
								north, tag, COMM_TOPO, request_array(10), ierr)
				! ------- RECEIVE FROM THE SOUTH, SEND TO THE SOUTH
				! Receive from the south, put it into the space for the ghost node 
				CALL MPI_IRECV( T(ind_high_y, ind_low_x+1), 1, NS_ROW_SENDRECV, &
								south, tag, COMM_TOPO, request_array(11), ierr)
				! Send to the south
				CALL MPI_ISEND( T(ind_high_y-1, ind_low_x+1), 1, NS_ROW_SENDRECV, &
								south, tag, COMM_TOPO, request_array(12), ierr)

				! Wait for data sends to complete before black points start referencing red points
				CALL MPI_WAITALL(12, request_array, status_array, ierr)
		

				!-------------------------------------------------------------------!
				! RESIDUALS
				! computing residuals
				rcurrent = SUM(sum(ABS(Told(resil:resih,resjl:resjh) - T(resil:resih,resjl:resjh)),1),1) &
									/ ((resih-resil+1)*(resjh-resjl+1))

				! Summing processor residuals to get global resiudal and broadcasting to get average
				call MPI_ALLREDUCE(rcurrent,rc,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_TOPO,ierr)
				rc = rc/nprocs

				!-------------------------------------------------------------------!
				! Printing to screen after a certain amount of time

				if ((MOD(iter,1000).eq.0)) then
					! TECPLOT

					! These subarrays are now sent to PID 0
					CALL MPI_ISEND( T, 1, SENDINGSUBBARAY, 0, tag2, COMM_TOPO, request_array_gather(1), ierr)	
							
					! Pid 0 receiving data and putting into final file
					IF (pid .EQ. 0) THEN

						DO i = 0, Nprocs-1
						
							! Creating receive subarray type bespoke to each processor
							CALL MPI_Type_create_subarray(2, [ny, nx], [subarray_rows_array(i+1), subarray_cols_array(i+1)], &
							[subarray_row_start(i+1),subarray_col_start(i+1)], MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, RECVSUBBARAY, ierr)
							
							CALL MPI_TYPE_COMMIT(RECVSUBBARAY, ierr)
						
							! Receiving with this new receiving subarray type
							CALL MPI_IRECV( Tfinal, 1, RECVSUBBARAY, &
								i, tag2, COMM_TOPO, request_array_gather(i+2), ierr)
						END DO

						CALL MPI_WAITALL(Nprocs+1, request_array_gather, status_array_gather, ierr)
						
						
						! --- PUTTING INTO FILE ---
						! x vector of whole domain 
						xtot 	= 0. + dx * [(i, i=0,(nx-1))] ! Implied DO loop used
						! y vector of whole domain 
						ytot 	= 0. + dy * [(i, i=0,(ny-1))] ! Implied DO loop used

						! Writing updated solution to file at each new iteration
						write(file_name, "(A9,I5,A4)") "TecPlot2D",int(time),".tec"
						CALL tecplot_2D ( iunit, nx, ny, xtot, ytot, Tfinal,  file_name )
					
						write(*,*) '-----------------------------------------'
						write(*,*) 'Time =', time
						write(*,*) 'Iteration =', iter
						write(*,*) 'Residual =', rc
						write(*,*) 'Average Temperature =', sum(sum((Tfinal),1),1)/(nx*ny)
						write(*,*) '-----------------------------------------'
						
					END IF
					
				end if


				! updating counter
				time = time + dt
				iter = iter + 1
			end do

		! CONJUGATE GRADIENT
        CASE("conj")

			! Initialising solution constants and vectors
			! Initialising delta and delta0
			delta = 1.0
			delta_o = 1.0
			! Initialising dot product variable
			dp = 0.0
			! Initialising aconst and beta
			aconst = 1.0
			beta = 1.0

			! Initialising A matrix constants
			do j = jl,jh	
				do i = il,ih
					ap(i,j) = 1
					an(i,j) = ((dt*alpha)/(dx*dx))
					as(i,j) = ((dt*alpha)/(dx*dx))
					ae(i,j) = ((dt*alpha)/(dx*dx))
					aw(i,j) = ((dt*alpha)/(dx*dx))
					b(i,j) = (1- (4*dt*alpha)/(dx*dx))
				end do
			end do

			allocate(dmat(il:ih,jl:jh))
			allocate(dmatT(il:ih,jl:jh))
			allocate(resmatT(il:ih,jl:jh))
			allocate(qmat(il:ih,jl:jh))
			dmat(:,:) = 0
			dmatT(:,:) = 0
			resmatT(:,:) = 0
			qmat(:,:) = 0

			! calculating initial grid residual 
		    do j = jl+1,jh-1
				do i = il+1,ih-1
					call respar(aw,ae,an,as,ap,b,T,il,ih,jl,jh,resmat)
				end do
			end do 


		!-------------------------------------------------------------------!
			! COMMUNICATION 
			! ------- RECEIVE FROM THE RIGHT/EAST, SEND TO THE RIGHT/EAST
			! Receive from the east1, and put it into the space for the ghost node on the east side
			CALL MPI_IRECV( resmat(ind_low_east1:ind_high_east1, ind_high_x), ncalcpoints_y_east1, MPI_DOUBLE_PRECISION, &
							east1, tag, COMM_TOPO, request_array(1), ierr)
			! Send to the east1
			CALL MPI_ISEND( resmat(ind_low_east1:ind_high_east1, ind_high_x-1), ncalcpoints_y_east1, MPI_DOUBLE_PRECISION, &
							east1, tag, COMM_TOPO, request_array(2), ierr)
							
			! Receive from the east2, and put it into the space for the ghost node on the east side
			CALL MPI_IRECV( resmat(ind_low_east2:ind_high_east2, ind_high_x), ncalcpoints_y_east2, MPI_DOUBLE_PRECISION, &
							east2, tag, COMM_TOPO, request_array(3), ierr)
			! Send to the east2
			CALL MPI_ISEND( resmat(ind_low_east2:ind_high_east2, ind_high_x-1), ncalcpoints_y_east2, MPI_DOUBLE_PRECISION, &
							east2, tag, COMM_TOPO, request_array(4), ierr)


			! ------- RECEIVE FROM THE LEFT/WEST, SEND TO THE LEFT/WEST
			! Receive from the west1, put it into the space for the ghost node on the west side 
			CALL MPI_IRECV( resmat(ind_low_west1:ind_high_west1, ind_low_x), ncalcpoints_y_west1, MPI_DOUBLE_PRECISION, & 
							west1, tag, COMM_TOPO, request_array(5), ierr)
			! Send to the west1
			CALL MPI_ISEND( resmat(ind_low_west1:ind_high_west1, ind_low_x+1), ncalcpoints_y_west1, MPI_DOUBLE_PRECISION, &
							west1, tag, COMM_TOPO, request_array(6), ierr)
							
			! Receive from the west2, put it into the space for the ghost node on the west side 
			CALL MPI_IRECV( resmat(ind_low_west2:ind_high_west2, ind_low_x), ncalcpoints_y_west2, MPI_DOUBLE_PRECISION, & 
							west2, tag, COMM_TOPO, request_array(7), ierr)
			! Send to the west2
			CALL MPI_ISEND( resmat(ind_low_west2:ind_high_west2, ind_low_x+1), ncalcpoints_y_west2, MPI_DOUBLE_PRECISION, &
							west2, tag, COMM_TOPO, request_array(8), ierr)
							
			! ------- RECEIVE FROM THE NORTH, SEND TO THE NORTH
			! Receive from the north, put it into the space for the ghost node 
			CALL MPI_IRECV( resmat(ind_low_y, ind_low_x+1), 1, NS_ROW_SENDRECV, &
							north, tag, COMM_TOPO, request_array(9), ierr)
			! Send to the north
			CALL MPI_ISEND( resmat(ind_low_y+1, ind_low_x+1), 1, NS_ROW_SENDRECV, &
							north, tag, COMM_TOPO, request_array(10), ierr)
			! ------- RECEIVE FROM THE SOUTH, SEND TO THE SOUTH
			! Receive from the south, put it into the space for the ghost node 
			CALL MPI_IRECV( resmat(ind_high_y, ind_low_x+1), 1, NS_ROW_SENDRECV, &
							south, tag, COMM_TOPO, request_array(11), ierr)
			! Send to the south
			CALL MPI_ISEND( resmat(ind_high_y-1, ind_low_x+1), 1, NS_ROW_SENDRECV, &
							south, tag, COMM_TOPO, request_array(12), ierr)

			! Wait for data sends to complete before black points start referencing red points
			CALL MPI_WAITALL(12, request_array, status_array, ierr)

			! Assigning the value of residual matrix to d
			dmat = resmat

			resmatT = transpose(resmat)
			! Compute dot product in order to determine delta
			dp = 0
			do j = jl,jh
				do i = il,ih
					dp = dp + resmatT(i,j)*resmat(i,j)
				end do
			end do

			delta_o = delta

			do while (rc>res_max)

				! Computring the q matrix
				do j = jl+1,jh-1
					do i = il+1,ih-1
						qmat(i,j) = aw(i,j)*dmat(i-1,j) + ae(i,j)*dmat(i+1,j) + an(i,j)*dmat(i,j+1) + as(i,j)*dmat(i,j-1) &
								+ ap(i,j)*dmat(i,j)
					end do
				end do

			!-------------------------------------------------------------------!
				! COMMUNICATION
				! ------- RECEIVE FROM THE RIGHT/EAST, SEND TO THE RIGHT/EAST
				! Receive from the east1, and put it into the space for the ghost node on the east side
				CALL MPI_IRECV( qmat(ind_low_east1:ind_high_east1, ind_high_x), ncalcpoints_y_east1, MPI_DOUBLE_PRECISION, &
								east1, tag, COMM_TOPO, request_array(1), ierr)
				! Send to the east1
				CALL MPI_ISEND( qmat(ind_low_east1:ind_high_east1, ind_high_x-1), ncalcpoints_y_east1, MPI_DOUBLE_PRECISION, &
								east1, tag, COMM_TOPO, request_array(2), ierr)
								
				! Receive from the east2, and put it into the space for the ghost node on the east side
				CALL MPI_IRECV( qmat(ind_low_east2:ind_high_east2, ind_high_x), ncalcpoints_y_east2, MPI_DOUBLE_PRECISION, &
								east2, tag, COMM_TOPO, request_array(3), ierr)
				! Send to the east2
				CALL MPI_ISEND( qmat(ind_low_east2:ind_high_east2, ind_high_x-1), ncalcpoints_y_east2, MPI_DOUBLE_PRECISION, &
								east2, tag, COMM_TOPO, request_array(4), ierr)


				! ------- RECEIVE FROM THE LEFT/WEST, SEND TO THE LEFT/WEST
				! Receive from the west1, put it into the space for the ghost node on the west side 
				CALL MPI_IRECV( qmat(ind_low_west1:ind_high_west1, ind_low_x), ncalcpoints_y_west1, MPI_DOUBLE_PRECISION, & 
								west1, tag, COMM_TOPO, request_array(5), ierr)
				! Send to the west1
				CALL MPI_ISEND( qmat(ind_low_west1:ind_high_west1, ind_low_x+1), ncalcpoints_y_west1, MPI_DOUBLE_PRECISION, &
								west1, tag, COMM_TOPO, request_array(6), ierr)
								
				! Receive from the west2, put it into the space for the ghost node on the west side 
				CALL MPI_IRECV( qmat(ind_low_west2:ind_high_west2, ind_low_x), ncalcpoints_y_west2, MPI_DOUBLE_PRECISION, & 
								west2, tag, COMM_TOPO, request_array(7), ierr)
				! Send to the west2
				CALL MPI_ISEND( qmat(ind_low_west2:ind_high_west2, ind_low_x+1), ncalcpoints_y_west2, MPI_DOUBLE_PRECISION, &
								west2, tag, COMM_TOPO, request_array(8), ierr)
								
				! ------- RECEIVE FROM THE NORTH, SEND TO THE NORTH
				! Receive from the north, put it into the space for the ghost node 
				CALL MPI_IRECV( qmat(ind_low_y, ind_low_x+1), 1, NS_ROW_SENDRECV, &
								north, tag, COMM_TOPO, request_array(9), ierr)
				! Send to the north
				CALL MPI_ISEND( qmat(ind_low_y+1, ind_low_x+1), 1, NS_ROW_SENDRECV, &
								north, tag, COMM_TOPO, request_array(10), ierr)
				! ------- RECEIVE FROM THE SOUTH, SEND TO THE SOUTH
				! Receive from the south, put it into the space for the ghost node 
				CALL MPI_IRECV( qmat(ind_high_y, ind_low_x+1), 1, NS_ROW_SENDRECV, &
								south, tag, COMM_TOPO, request_array(11), ierr)
				! Send to the south
				CALL MPI_ISEND( qmat(ind_high_y-1, ind_low_x+1), 1, NS_ROW_SENDRECV, &
								south, tag, COMM_TOPO, request_array(12), ierr)

				! Wait for data sends to complete before black points start referencing red points
				CALL MPI_WAITALL(12, request_array, status_array, ierr)

				! Computing alpha
				dp = 0
				dmatT = transpose(dmat)
				do j = jl,jh
					do i = il,ih
						dp = dp + dmatT(i,j)*qmat(i,j)
					end do
				end do
				aconst = delta/dp
				
				Told = T
				! Updating solution
				do j = jl+1,jh-1
					do i = il+1,ih-1
						T(i,j) = T(i,j) + aconst*dmat(i,j)
					end do
				end do
!-------------------------------------------------------------------!
				! COMMUNICATION
				! ------- RECEIVE FROM THE RIGHT/EAST, SEND TO THE RIGHT/EAST
				! Receive from the east1, and put it into the space for the ghost node on the east side
				CALL MPI_IRECV( T(ind_low_east1:ind_high_east1, ind_high_x), ncalcpoints_y_east1, MPI_DOUBLE_PRECISION, &
								east1, tag, COMM_TOPO, request_array(1), ierr)
				! Send to the east1
				CALL MPI_ISEND( T(ind_low_east1:ind_high_east1, ind_high_x-1), ncalcpoints_y_east1, MPI_DOUBLE_PRECISION, &
								east1, tag, COMM_TOPO, request_array(2), ierr)
								
				! Receive from the east2, and put it into the space for the ghost node on the east side
				CALL MPI_IRECV( T(ind_low_east2:ind_high_east2, ind_high_x), ncalcpoints_y_east2, MPI_DOUBLE_PRECISION, &
								east2, tag, COMM_TOPO, request_array(3), ierr)
				! Send to the east2
				CALL MPI_ISEND( T(ind_low_east2:ind_high_east2, ind_high_x-1), ncalcpoints_y_east2, MPI_DOUBLE_PRECISION, &
								east2, tag, COMM_TOPO, request_array(4), ierr)


				! ------- RECEIVE FROM THE LEFT/WEST, SEND TO THE LEFT/WEST
				! Receive from the west1, put it into the space for the ghost node on the west side 
				CALL MPI_IRECV( T(ind_low_west1:ind_high_west1, ind_low_x), ncalcpoints_y_west1, MPI_DOUBLE_PRECISION, & 
								west1, tag, COMM_TOPO, request_array(5), ierr)
				! Send to the west1
				CALL MPI_ISEND( T(ind_low_west1:ind_high_west1, ind_low_x+1), ncalcpoints_y_west1, MPI_DOUBLE_PRECISION, &
								west1, tag, COMM_TOPO, request_array(6), ierr)
								
				! Receive from the west2, put it into the space for the ghost node on the west side 
				CALL MPI_IRECV( T(ind_low_west2:ind_high_west2, ind_low_x), ncalcpoints_y_west2, MPI_DOUBLE_PRECISION, & 
								west2, tag, COMM_TOPO, request_array(7), ierr)
				! Send to the west2
				CALL MPI_ISEND( T(ind_low_west2:ind_high_west2, ind_low_x+1), ncalcpoints_y_west2, MPI_DOUBLE_PRECISION, &
								west2, tag, COMM_TOPO, request_array(8), ierr)
								
				! ------- RECEIVE FROM THE NORTH, SEND TO THE NORTH
				! Receive from the north, put it into the space for the ghost node 
				CALL MPI_IRECV( T(ind_low_y, ind_low_x+1), 1, NS_ROW_SENDRECV, &
								north, tag, COMM_TOPO, request_array(9), ierr)
				! Send to the north
				CALL MPI_ISEND( T(ind_low_y+1, ind_low_x+1), 1, NS_ROW_SENDRECV, &
								north, tag, COMM_TOPO, request_array(10), ierr)
				! ------- RECEIVE FROM THE SOUTH, SEND TO THE SOUTH
				! Receive from the south, put it into the space for the ghost node 
				CALL MPI_IRECV( T(ind_high_y, ind_low_x+1), 1, NS_ROW_SENDRECV, &
								south, tag, COMM_TOPO, request_array(11), ierr)
				! Send to the south
				CALL MPI_ISEND( T(ind_high_y-1, ind_low_x+1), 1, NS_ROW_SENDRECV, &
								south, tag, COMM_TOPO, request_array(12), ierr)

				! Wait for data sends to complete before black points start referencing red points
				CALL MPI_WAITALL(12, request_array, status_array, ierr)


				! Updating residual
				if ((MOD(iter,50).eq.0)) then
					do j = jl+1,jh-1
						do i = il+1,ih-1
							call respar(aw,ae,an,as,ap,b,T,il,ih,jl,jh,resmat)
						end do
					end do
				else 
					do j = jl+1,jh-1
						do i = il+1,ih-1
							resmat(i,j) = resmat(i,j) - aconst*qmat(i,j)
						end do
					end do
				end if


				rcurrent = SUM(SUM(ABS(resmat(il:ih,jl:jh)),1),1) / ((ih-il+1)*(jh-jl+1))
				! ! Summing processor residuals to get global resiudal and broadcasting to get average
				call MPI_ALLREDUCE(rcurrent,rc,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_TOPO,ierr)
				rc = rc/nprocs	

!-------------------------------------------------------------------!
				! COMMUNICATION 
				! ------- RECEIVE FROM THE RIGHT/EAST, SEND TO THE RIGHT/EAST
				! Receive from the east1, and put it into the space for the ghost node on the east side
				CALL MPI_IRECV( resmat(ind_low_east1:ind_high_east1, ind_high_x), ncalcpoints_y_east1, MPI_DOUBLE_PRECISION, &
								east1, tag, COMM_TOPO, request_array(1), ierr)
				! Send to the east1
				CALL MPI_ISEND( resmat(ind_low_east1:ind_high_east1, ind_high_x-1), ncalcpoints_y_east1, MPI_DOUBLE_PRECISION, &
								east1, tag, COMM_TOPO, request_array(2), ierr)
								
				! Receive from the east2, and put it into the space for the ghost node on the east side
				CALL MPI_IRECV( resmat(ind_low_east2:ind_high_east2, ind_high_x), ncalcpoints_y_east2, MPI_DOUBLE_PRECISION, &
								east2, tag, COMM_TOPO, request_array(3), ierr)
				! Send to the east2
				CALL MPI_ISEND( resmat(ind_low_east2:ind_high_east2, ind_high_x-1), ncalcpoints_y_east2, MPI_DOUBLE_PRECISION, &
								east2, tag, COMM_TOPO, request_array(4), ierr)


				! ------- RECEIVE FROM THE LEFT/WEST, SEND TO THE LEFT/WEST
				! Receive from the west1, put it into the space for the ghost node on the west side 
				CALL MPI_IRECV( resmat(ind_low_west1:ind_high_west1, ind_low_x), ncalcpoints_y_west1, MPI_DOUBLE_PRECISION, & 
								west1, tag, COMM_TOPO, request_array(5), ierr)
				! Send to the west1
				CALL MPI_ISEND( resmat(ind_low_west1:ind_high_west1, ind_low_x+1), ncalcpoints_y_west1, MPI_DOUBLE_PRECISION, &
								west1, tag, COMM_TOPO, request_array(6), ierr)
								
				! Receive from the west2, put it into the space for the ghost node on the west side 
				CALL MPI_IRECV( resmat(ind_low_west2:ind_high_west2, ind_low_x), ncalcpoints_y_west2, MPI_DOUBLE_PRECISION, & 
								west2, tag, COMM_TOPO, request_array(7), ierr)
				! Send to the west2
				CALL MPI_ISEND( resmat(ind_low_west2:ind_high_west2, ind_low_x+1), ncalcpoints_y_west2, MPI_DOUBLE_PRECISION, &
								west2, tag, COMM_TOPO, request_array(8), ierr)
								
				! ------- RECEIVE FROM THE NORTH, SEND TO THE NORTH
				! Receive from the north, put it into the space for the ghost node 
				CALL MPI_IRECV( resmat(ind_low_y, ind_low_x+1), 1, NS_ROW_SENDRECV, &
								north, tag, COMM_TOPO, request_array(9), ierr)
				! Send to the north
				CALL MPI_ISEND( resmat(ind_low_y+1, ind_low_x+1), 1, NS_ROW_SENDRECV, &
								north, tag, COMM_TOPO, request_array(10), ierr)
				! ------- RECEIVE FROM THE SOUTH, SEND TO THE SOUTH
				! Receive from the south, put it into the space for the ghost node 
				CALL MPI_IRECV( resmat(ind_high_y, ind_low_x+1), 1, NS_ROW_SENDRECV, &
								south, tag, COMM_TOPO, request_array(11), ierr)
				! Send to the south
				CALL MPI_ISEND( resmat(ind_high_y-1, ind_low_x+1), 1, NS_ROW_SENDRECV, &
								south, tag, COMM_TOPO, request_array(12), ierr)

				! Wait for data sends to complete before black points start referencing red points
				CALL MPI_WAITALL(12, request_array, status_array, ierr)

				! Updating delta0
				delta_o = delta

				! Updating delta
				dp = 0
				resmatT = transpose(resmat)
				do j = jl+1,jh-1
					do i = il+1,ih-1
						dp = dp + resmatT(i,j)*resmat(i,j)
					end do
				end do
				delta = dp

				! Updating beta
				beta = delta/delta_o

				! Updating d
				do j = jl+1,jh-1
					do i = il+1,ih-1
						dmat(i,j) = resmat(i,j) + beta*dmat(i,j)
					end do
				end do

!-------------------------------------------------------------------!
				! COMMUNICATION 
				! ------- RECEIVE FROM THE RIGHT/EAST, SEND TO THE RIGHT/EAST
				! Receive from the east1, and put it into the space for the ghost node on the east side
				CALL MPI_IRECV( dmat(ind_low_east1:ind_high_east1, ind_high_x), ncalcpoints_y_east1, MPI_DOUBLE_PRECISION, &
								east1, tag, COMM_TOPO, request_array(1), ierr)
				! Send to the east1
				CALL MPI_ISEND( dmat(ind_low_east1:ind_high_east1, ind_high_x-1), ncalcpoints_y_east1, MPI_DOUBLE_PRECISION, &
								east1, tag, COMM_TOPO, request_array(2), ierr)
								
				! Receive from the east2, and put it into the space for the ghost node on the east side
				CALL MPI_IRECV( dmat(ind_low_east2:ind_high_east2, ind_high_x), ncalcpoints_y_east2, MPI_DOUBLE_PRECISION, &
								east2, tag, COMM_TOPO, request_array(3), ierr)
				! Send to the east2
				CALL MPI_ISEND( dmat(ind_low_east2:ind_high_east2, ind_high_x-1), ncalcpoints_y_east2, MPI_DOUBLE_PRECISION, &
								east2, tag, COMM_TOPO, request_array(4), ierr)


				! ------- RECEIVE FROM THE LEFT/WEST, SEND TO THE LEFT/WEST
				! Receive from the west1, put it into the space for the ghost node on the west side 
				CALL MPI_IRECV( dmat(ind_low_west1:ind_high_west1, ind_low_x), ncalcpoints_y_west1, MPI_DOUBLE_PRECISION, & 
								west1, tag, COMM_TOPO, request_array(5), ierr)
				! Send to the west1
				CALL MPI_ISEND( dmat(ind_low_west1:ind_high_west1, ind_low_x+1), ncalcpoints_y_west1, MPI_DOUBLE_PRECISION, &
								west1, tag, COMM_TOPO, request_array(6), ierr)
								
				! Receive from the west2, put it into the space for the ghost node on the west side 
				CALL MPI_IRECV( dmat(ind_low_west2:ind_high_west2, ind_low_x), ncalcpoints_y_west2, MPI_DOUBLE_PRECISION, & 
								west2, tag, COMM_TOPO, request_array(7), ierr)
				! Send to the west2
				CALL MPI_ISEND( dmat(ind_low_west2:ind_high_west2, ind_low_x+1), ncalcpoints_y_west2, MPI_DOUBLE_PRECISION, &
								west2, tag, COMM_TOPO, request_array(8), ierr)
								
				! ------- RECEIVE FROM THE NORTH, SEND TO THE NORTH
				! Receive from the north, put it into the space for the ghost node 
				CALL MPI_IRECV( dmat(ind_low_y, ind_low_x+1), 1, NS_ROW_SENDRECV, &
								north, tag, COMM_TOPO, request_array(9), ierr)
				! Send to the north
				CALL MPI_ISEND( dmat(ind_low_y+1, ind_low_x+1), 1, NS_ROW_SENDRECV, &
								north, tag, COMM_TOPO, request_array(10), ierr)
				! ------- RECEIVE FROM THE SOUTH, SEND TO THE SOUTH
				! Receive from the south, put it into the space for the ghost node 
				CALL MPI_IRECV( dmat(ind_high_y, ind_low_x+1), 1, NS_ROW_SENDRECV, &
								south, tag, COMM_TOPO, request_array(11), ierr)
				! Send to the south
				CALL MPI_ISEND( dmat(ind_high_y-1, ind_low_x+1), 1, NS_ROW_SENDRECV, &
								south, tag, COMM_TOPO, request_array(12), ierr)

				! Wait for data sends to complete before black points start referencing red points
				CALL MPI_WAITALL(12, request_array, status_array, ierr)

				!-------------------------------------------------------------------!
				! Printing to screen after a certain amount of time
				! if ((time-int(time))<dt) then
				if ((MOD(iter,100).eq.0)) then
					! TECPLOT

					! These subarrays are now sent to PID 0
					CALL MPI_ISEND( T, 1, SENDINGSUBBARAY, 0, tag2, COMM_TOPO, request_array_gather(1), ierr)	
							
					! Pid 0 receiving data and putting into final file
					IF (pid .EQ. 0) THEN

						DO i = 0, Nprocs-1
						
							! Creating receive subarray type bespoke to each processor
							CALL MPI_Type_create_subarray(2, [ny, nx], [subarray_rows_array(i+1), subarray_cols_array(i+1)], &
							[subarray_row_start(i+1),subarray_col_start(i+1)], MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, RECVSUBBARAY, ierr)
							
							CALL MPI_TYPE_COMMIT(RECVSUBBARAY, ierr)
						
							! Receiving with this new receiving subarray type
							CALL MPI_IRECV( Tfinal, 1, RECVSUBBARAY, &
								i, tag2, COMM_TOPO, request_array_gather(i+2), ierr)
						END DO

						CALL MPI_WAITALL(Nprocs+1, request_array_gather, status_array_gather, ierr)
						
						! --- PUTTING INTO FILE ---
						! x vector of whole domain (could have also done a mpi_gatherv)
						xtot 	= 0. + dx * [(i, i=0,(nx-1))] ! Implied DO loop used
						! y vector of whole domain (could have also done a mpi_gatherv)
						ytot 	= 0. + dy * [(i, i=0,(ny-1))] ! Implied DO loop used

						! Writing updated solution to file at each new iteration
						write(file_name, "(A9,I5,A4)") "TecPlot2D",int(time),".tec"
						CALL tecplot_2D ( iunit, nx, ny, xtot, ytot, Tfinal,  file_name )
				
						write(*,*) '-----------------------------------------'
						write(*,*) 'Time =', time
						write(*,*) 'Iteration =', iter
						write(*,*) 'Residual =', rc
						write(*,*) 'Average Temperature =', sum(sum((Tfinal),1),1)/(nx*ny)
						write(*,*) '-----------------------------------------'
						
					end if
					
				end if

				! Updating timer and iteration counter
				time = time + dt
				iter = iter + 1

			end do

    	CASE DEFAULT 
		WRITE(*,*) "No solver selected or incorrect selection"
		STOP

    END SELECT


	if (pid.eq.0) then
		write(*,*) 'Output after solver'
		write(*,*) '-----------------------------------------'
		write(*,*) 'Time =', time
		write(*,*) 'Iteration =', iter
		write(*,*) 'Residual =', rc
		write(*,*) 'CFL =', CFL, 'Iplicit Accuracy =', uncondstab
		t2 = MPI_WTIME()
		write(*,*) 'Computational time =', t2-t1, 'seconds'
		write(*,*) 'Average Temperature =', sum(sum((Tfinal),1),1)/(nx*ny)
		write(*,*) '-----------------------------------------'
	end if


CALL MPI_FINALIZE(ierr)

END PROGRAM MAIN

