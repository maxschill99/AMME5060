! Assignment 1
! SID: 	470378812
! Name: Clara Morris
! Due: 	14/9/2022

! Compiling Instructions:
! Run the following code to compile:
! 		mpif90 -o finalcode outputmod.f90 nodemod.f90 solvmod.f90 varmod.f90 morris.f90 
! 		NOTE: 	You will first have to create the object files.
! 				On my computer I can achieve this by running the above line of code twice. 
! 				The first time throws an error (but creates the object files), and the second time will compile everything.
! Run the following code to execute the program on four processors:
! 		mpirun -np 4 ./finalcode

! Program Description:
! This program solves for the temperature in a simple two dimensional advection diffusion problem
! Second order centered finite difference discretisation in space


PROGRAM Assignment1

	! *************************************************
	!!! --- SET UP ---

	! Linking modules
	USE variablemodule 		! Contains variable allocation all problem-specific variable values 
	USE nodemodule			! Obtains nodes and indices of the domain considered by each processor
	USE solvingmodule		! Finds residual and solves for Tk_1 from Tk from each processor
	USE outputmodule		! Creates tec file (base code was given but has been modified)
	
	! Ensuring all names in the program must be explicitly declared
	IMPLICIT NONE
	
	! Allowing MPI to be used 
	INCLUDE 'mpif.h'
		
	! Variable declaration and assignment, MPI initialising
	CALL intialise()
	
	
	! *************************************************
	!!! --- INITIALISE MPI ---
	CALL MPI_INIT(ierr) ! Initialises everything 
	CALL MPI_COMM_RANK(MPI_COMM_WORLD, pid_global, ierr) ! Getting processor ID number
	CALL MPI_COMM_SIZE(MPI_COMM_WORLD, Nprocs, ierr) ! Getting number of total processors in global communicator
	
	! Start timing the program with MPI wall timer
	t1 = MPI_WTIME()
	
	
	! *************************************************
	!!! --- DOMAIN DECOMPOSITION WITH CARTEISAN TOPOLOGY ---
	
	! To minimise communication, the domain is split up in 2 dimensions if possible.
	! First, it needs to be checked if splitting up in 2 dimensions is even possible.
	! If the number of processors is prime, the domain will have to be split into slabs.

    ! Check what the greatest divisor of the number of processors is:
	bestdiv = 1 			! 1 guaranteed
	bestgap = Nprocs-1		! Difference between divisors. Ideally this would be zero (perfect square).
	
	! Start from 2 (A divisor of 1 is already a guarantee) and just go up to half the number of procs,
	! or nearest integer close to half way if Nprocs is odd.
	DO divisor = 2, FLOOR(DBLE(Nprocs)/2.) 
		IF (MOD(Nprocs,divisor) .EQ. 0) THEN ! If it's a divisor
		
			IF (ABS(divisor - Nprocs/divisor) .LE. bestgap) THEN 
				bestgap = ABS(divisor - Nprocs/divisor)
				bestdiv = divisor ! Looking for 2 divisors closest together (take biggest from those 2 if Nprocs not square)
			END IF
			
		END IF
	END DO
	
	ndims = 2 ! Partition in 2 dimensions
	ALLOCATE( dims(ndims), periods(ndims), coords(ndims))
	
	IF (bestdiv .EQ. 1) THEN ! The number is prime and the domain cannot be partitioned in 2 dimensions
		dims(1) = Nprocs ! Number of procs in this one dimension is simply all of them (number of procs in x direction)
		dims(2) = 1 ! Number of procs in y direction is 1
	ELSE
		dims(2) = Nprocs/bestdiv ! number of processors in y direction
		dims(1) = bestdiv ! This number will always be greater than or equal to dim(2) (more columns). Number of procs in x direction
	END IF
	
	! The domain is not periodic
	periods = .false.
	
	! Create a new cartesian communicator based on the above analysis
	CALL MPI_CART_CREATE(MPI_COMM_WORLD, ndims, dims, periods, .true., COMM_CART, ierr)
	
	! New processor ID number
	CALL MPI_COMM_RANK(COMM_CART, pid, ierr)

	!!! --- GATHERING BEARINGS ---
	! Get coordinates for each processor
	CALL MPI_CART_COORDS(COMM_CART, pid, ndims, coords, ierr)
	! coords is an array with pid location in (x, y) but y counts from top down and counts start from zero
	
	! Send/Recv id's for each processor horizontally and vertically (if applicable)
	CALL MPI_CART_SHIFT(COMM_CART, 0, 1, left, right, ierr)
	CALL MPI_CART_SHIFT(COMM_CART, 1, 1, north, south, ierr)
	
	
	! *************************************************
	!!! --- NODES CONSIDERED BY EACH PROCESSOR ---
	
	! Nodes and indices in x
	!get_nodes(Ntpoints, Nprocs, coord_id, node_low, node_high)
	CALL get_nodes(nx, dims(1), coords(1), node_low_x, node_high_x)
	CALL get_indices(dims(1), coords(1), node_low_x, node_high_x, ind_low_x, ind_high_x, ncalcpoints_x)
	
	! Nodes and indices in y
	CALL get_nodes(ny, dims(2), coords(2), node_low_y, node_high_y)
	CALL get_indices(dims(2), coords(2), node_low_y, node_high_y, ind_low_y, ind_high_y, ncalcpoints_y)
	
	
	! *************************************************
	!!! --- SOLVING PREPARATION ---
	
	! Allocate space in matrices for each processor based on their indices
	! Indices in y correspond to row location, while indices in x correspond to column location
	ALLOCATE( x(ind_low_x:ind_high_x), y(ind_low_y:ind_high_y) )
	ALLOCATE( Tk_1(ind_low_y:ind_high_y, ind_low_x:ind_high_x), Tk(ind_low_y:ind_high_y, ind_low_x:ind_high_x) )
	ALLOCATE( U(ind_low_y:ind_high_y, ind_low_x:ind_high_x), V(ind_low_y:ind_high_y, ind_low_x:ind_high_x) )

	! Defining x and y-coordinates considered by each processor (including ghost nodes)
	! In total, values are linearly spaced from 0 to Lx in increments of dx
	x 	= (ind_low_x-1)*dx + dx * [(i, i=0,(ind_high_x-ind_low_x))] ! Implied DO loop used
		
	! In total, values are linearly spaced from 0 to Ly in increments of dy
	! y vector created from 1 to 0 so that at y(1) (top of domain), the axis system is respected
	y 	= (ny-ind_low_y)*dy  -  dy * [(i, i=0,(ind_high_y-ind_low_y))] ! Implied DO loop used
	
	
	! Velocities U and V
	DO j = ind_low_x, ind_high_x
		! then loop through by row 
		DO i = ind_low_y,ind_high_y
			U(i,j) =      ( SIN((pi/Lx)*x(j)) * COS((pi/Ly)*y(i)) ) + ( SIN((2*pi/Lx)*x(j)) * COS((2*pi/Ly)*y(i)) ) &
						+ ( SIN((4*pi/Lx)*x(j)) * COS((4*pi/Ly)*y(i)) )
						
			V(i,j) = (-1)*( COS((pi/Lx)*x(j)) * SIN((pi/Ly)*y(i)) ) - ( COS((2*pi/Lx)*x(j)) * SIN((2*pi/Ly)*y(i)) ) &
						- ( COS((4*pi/Lx)*x(j)) * SIN((4*pi/Ly)*y(i)) )
		END DO
	END DO
	
	!!! --- BOUNDARY CONDITIONS ---
	! Initial temperature distribution guess. Will be somewhere less than 1 and close to 0 based on BCs
	! Just guess zero everywhere
	Tk  	= 0.0
	
	! Boundary conditions along edges
	IF (ind_low_y .EQ. 1) THEN ! Top edge of domain, heated
		Tk(1, :) 	= SIN( (pi/Lx) * x)
	END IF
	IF (ind_high_y .EQ. ny) THEN ! Bottom edge of domain, 0 deg
		Tk(ny, :) 	=  0.
	END IF
	IF (ind_low_x .EQ. 1) THEN ! Left edge of domain, 0 deg
		Tk(:, 1) 	= 0.
	END IF
	IF (ind_high_x .EQ. nx) THEN ! Right edge of domain, 0 deg
		Tk(:, nx) 	= 0.
	END IF
	
	! Setting up the k+1 iteration matrix to be the same as the k matrix.
	Tk_1 = Tk
	
	! This ensures:
	! 	a) the BC's are enforced across all iterations
	! 	b) the whole field is zero unless at top edge. This takes care of the corner "ghost ghost" nodes 
	!	   that will never get used or updated by send/recv. They are just set and kept at zero.
	
	
	! --- SETTING UP SPECIFIC DATA TYPE ---
	! New vector type for sending and receiving the rows north and south (data not contiguous along a row)
	CALL MPI_TYPE_VECTOR(ncalcpoints_x, 1, ncalcpoints_y+2, MPI_DOUBLE_PRECISION, NS_ROW_SENDRECV , ierr)
	CALL MPI_TYPE_COMMIT(NS_ROW_SENDRECV, ierr)
	

	! --- SETTING UP RESIDUAL ---
	! Set up a residual
	ALLOCATE ( Res(node_low_y:node_high_y, node_low_x:node_high_x) )
	Res = 0.0 ! BCs known so all the residuals of nodes along edges are zero
	
	CALL get_residual_finalval() ! Obtaining 1/n * sum of abs value of residuals at each point for each processor 
	! Distributing global residual to all processors
	CALL MPI_ALLREDUCE(res_proc_val, res_final_val, 1, MPI_DOUBLE_PRECISION, MPI_SUM, COMM_CART, ierr)
	
	! --- SEND/RECV ARRAYS ---
	! Send/recv status and request arrays
	ALLOCATE( status_array(MPI_STATUS_SIZE, 8) )
	ALLOCATE( request_array(8) )
	
	! *************************************************
	!!! --- SOLVING ---
	

	! While overall residual too big (>1E-3) and max iterations hasnt been exceeded
	DO WHILE ((k .LE. iter_max) .AND. (res_final_val .GT. res_max))
			
			! ----------------------------- "REDS" -----------------------------------
			! Calculate next iteration
			CALL get_Tk1_Red() ! Fills in Tk_1 matrix with updated reds using last iteration
			
			! ------- RECEIVE FROM THE RIGHT, SEND TO THE RIGHT
			! Receive from the right, and put it into the space for the ghost node 
			CALL MPI_IRECV( Tk_1(ind_low_y+1:ind_high_y-1, ind_high_x), ncalcpoints_y, MPI_DOUBLE_PRECISION, &
							right, tag, COMM_CART, request_array(1), ierr)
			! Send to the right
			CALL MPI_ISEND( Tk_1(ind_low_y+1:ind_high_y-1, ind_high_x-1), ncalcpoints_y, MPI_DOUBLE_PRECISION, &
							right, tag, COMM_CART, request_array(2), ierr)

			! ------- RECEIVE FROM THE LEFT, SEND TO THE LEFT
			! Receive from the left, put it into the space for the ghost node 
			CALL MPI_IRECV( Tk_1(ind_low_y+1:ind_high_y-1, ind_low_x), ncalcpoints_y, MPI_DOUBLE_PRECISION, & 
							left, tag, COMM_CART, request_array(3), ierr)
			! Send to the left
			CALL MPI_ISEND( Tk_1(ind_low_y+1:ind_high_y-1, ind_low_x+1), ncalcpoints_y, MPI_DOUBLE_PRECISION, &
							left, tag, COMM_CART, request_array(4), ierr)
	
			! ------- RECEIVE FROM THE NORTH, SEND TO THE NORTH
			! Receive from the north, put it into the space for the ghost node 
			CALL MPI_IRECV( Tk_1(ind_low_y, ind_low_x+1), 1, NS_ROW_SENDRECV, &
							north, tag, COMM_CART, request_array(5), ierr)
			! Send to the north
			CALL MPI_ISEND( Tk_1(ind_low_y+1, ind_low_x+1), 1, NS_ROW_SENDRECV, &
							north, tag, COMM_CART, request_array(6), ierr)

			! ------- RECEIVE FROM THE SOUTH, SEND TO THE SOUTH
			! Receive from the south, put it into the space for the ghost node 
			CALL MPI_IRECV( Tk_1(ind_high_y, ind_low_x+1), 1, NS_ROW_SENDRECV, &
							south, tag, COMM_CART, request_array(7), ierr)
			! Send to the south
			CALL MPI_ISEND( Tk_1(ind_high_y-1, ind_low_x+1), 1, NS_ROW_SENDRECV, &
							south, tag, COMM_CART, request_array(8), ierr)

			! Wait for data sends to complete before black points start referencing red points
			CALL MPI_WAITALL(8, request_array, status_array, ierr)
			
			
			! ----------------------------- "BLACKS" -----------------------------------
			! Calculate next iteration
			CALL get_Tk1_Black() ! Fills in last spots of Tk_1 matrix with the already updated Tk_1 red spots
			
			! Swapping T_n1 ghost nodes via SEND/RECV commands (none of this will apply if Nprocs = 1)
			! ------- RECEIVE FROM THE RIGHT, SEND TO THE RIGHT
			! Receive from the right, and put it into the space for the ghost node 
			CALL MPI_IRECV( Tk_1(ind_low_y+1:ind_high_y-1, ind_high_x), ncalcpoints_y, MPI_DOUBLE_PRECISION, &
							right, tag, COMM_CART, request_array(1), ierr)
			! Send to the right
			CALL MPI_ISEND( Tk_1(ind_low_y+1:ind_high_y-1, ind_high_x-1), ncalcpoints_y, MPI_DOUBLE_PRECISION, &
							right, tag, COMM_CART, request_array(2), ierr)

			! ------- RECEIVE FROM THE LEFT, SEND TO THE LEFT
			! Receive from the left, put it into the space for the ghost node 
			CALL MPI_IRECV( Tk_1(ind_low_y+1:ind_high_y-1, ind_low_x), ncalcpoints_y, MPI_DOUBLE_PRECISION, & 
							left, tag, COMM_CART, request_array(3), ierr)
			! Send to the left
			CALL MPI_ISEND( Tk_1(ind_low_y+1:ind_high_y-1, ind_low_x+1), ncalcpoints_y, MPI_DOUBLE_PRECISION, &
							left, tag, COMM_CART, request_array(4), ierr)
	
			! ------- RECEIVE FROM THE NORTH, SEND TO THE NORTH
			! Receive from the north, put it into the space for the ghost node 
			CALL MPI_IRECV( Tk_1(ind_low_y, ind_low_x+1), 1, NS_ROW_SENDRECV, &
							north, tag, COMM_CART, request_array(5), ierr)
			! Send to the north
			CALL MPI_ISEND( Tk_1(ind_low_y+1, ind_low_x+1), 1, NS_ROW_SENDRECV, &
							north, tag, COMM_CART, request_array(6), ierr)

			! ------- RECEIVE FROM THE SOUTH, SEND TO THE SOUTH
			! Receive from the south, put it into the space for the ghost node 
			CALL MPI_IRECV( Tk_1(ind_high_y, ind_low_x+1), 1, NS_ROW_SENDRECV, &
							south, tag, COMM_CART, request_array(7), ierr)
			! Send to the south
			CALL MPI_ISEND( Tk_1(ind_high_y-1, ind_low_x+1), 1, NS_ROW_SENDRECV, &
							south, tag, COMM_CART, request_array(8), ierr)

			! Can't change the buffer (by starting next loop) until the sends are complete.
			! Can't assign Tk to new iteration until the proper data has come in. Need to wait.
			CALL MPI_WAITALL(8, request_array, status_array, ierr)


			! Now Tk_1 is populated with new values. Update Tk and Tk1
			Tk = Tk_1
			
			! ----------------------------- RESIDUALS -----------------------------------
			! For columns node_low_x+1, ind_high_x-1 and rows ind_low_y+1,ind_high_y-1
			CALL get_residual_finalval() ! Obtaining sum of abs value of residuals at each point for each processor 
			! gets out res_proc_val
			
			CALL MPI_ALLREDUCE(res_proc_val, res_final_val, 1, MPI_DOUBLE_PRECISION, MPI_SUM, COMM_CART, ierr)		
		
			! Have one processor write to screen every so often so that user can see residual come down and knows simulation is working
			IF (pid == 0) THEN
				IF (mod(k, 400) == 0) THEN
				WRITE(*,*) "Iteration number: ", k, "Residual: ", res_final_val
				END IF
			END IF
			
			! Update iteration
			k = k+1
	
	! End while loop
	END DO
	
	! Write to the screen if the simulation didn't converge to ideal value to let user know
	IF ((k .GE. iter_max) .AND. (pid == 0)) THEN
		WRITE(*,'(A52, I8, A10)') "Residual target not reached. Maximum iterations of ",iter_max, " exceeded."
	END IF
	
	
	! *************************************************
	! --- Collating Final Temperature Arrays ---
	! Tk arrays from node values need to be sent to processor zero who needs to put them in the right place.
	! Sending must consider that ghost nodes need to be ignored, and that looks different for all processors
	! Receving will also need to consider where the matrix goes in the final domain (will not be contiguous)
	
	! Unfortunately, GATHERv will only work for subarrays of equal size, so individual send/recvs must be used.
	! Because the size of the arrays may differ processor-to-processor (as any number of points in x and y are allowed
	! which means extra nodes may be given to some processors), varying receive-subarray-types have to be created.
	! GATHERv does not like this. This seems to be the price that is paid for having code that works for the general case.
	
	! ---- MAXIMUM AND AVERAGE TEMPERATURE ----
	! Printing maximum and average temperature to 8 decimal places:
	
	! Max temperature. Can take max over ghost nodes too because it doesn't matter - after global max.
	CALL MPI_REDUCE(maxval(Tk), &
	T_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, COMM_CART, ierr)
	
	! Average temperature. Only calculate average over nodes considered by each processor. Ignore ghost nodes.
	CALL MPI_REDUCE( (1/(DBLE(nx)*DBLE(ny))) * sum( Tk(node_low_y:node_high_y, node_low_x:node_high_x) ), &
 	T_avg, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, COMM_CART, ierr)
	
	IF (pid .EQ. 0) THEN 
		WRITE(*,'(A35, f10.8, A16)') "Maximum temperature across domain: ", T_max, " degrees Celcius"
		WRITE(*,'(A35, f10.8, A16)') "Average temperature across domain: ", T_avg, " degrees Celcius"
	END IF
	
	! ---- SEND DATA TYPE CREATION ----
	! The final send/recvs will be dependent on number of processors unlike before. Re-allocate relevant arrays to reflect this
	DEALLOCATE( status_array  )
	DEALLOCATE( request_array )
	ALLOCATE( status_array(MPI_STATUS_SIZE, Nprocs+1) )
	ALLOCATE( request_array(Nprocs+1) )
	
	! All processors first create subarrays. These are the final temperature array cleansed of the extra ghost nodes
	subarray_Nrows = node_high_y-node_low_y + 1
	subarray_Ncols = node_high_x-node_low_x + 1
	
	CALL MPI_Type_create_subarray(2, [ncalcpoints_y+2,ncalcpoints_x+2], [subarray_Nrows,subarray_Ncols], &
	[node_low_y-ind_low_y, node_low_x-ind_low_x], MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, SENDINGSUBBARAY, ierr)
	CALL MPI_TYPE_COMMIT(SENDINGSUBBARAY, ierr)
	
	! These subarrays are now sent to PID 0
	CALL MPI_ISEND( Tk, 1, SENDINGSUBBARAY, 0, tag2, COMM_CART, request_array(1), ierr)
		
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
		CALL MPI_GATHER( subarray_Nrows, 1, MPI_INTEGER, subarray_rows_array, 1, MPI_INTEGER, 0, COMM_CART, ierr )
		! Subarray column sizes
		CALL MPI_GATHER( subarray_Ncols, 1, MPI_INTEGER, subarray_cols_array, 1, MPI_INTEGER, 0, COMM_CART, ierr )
		! Start row location in the final matrix (minus 1 because start location starts at zero)
		CALL MPI_GATHER( node_low_y-1, 	 1, MPI_INTEGER, subarray_row_start,  1, MPI_INTEGER, 0, COMM_CART, ierr )
		! Start column location in the final matrix (minus 1 because start location starts at zero)
		CALL MPI_GATHER( node_low_x-1, 	 1, MPI_INTEGER, subarray_col_start,  1, MPI_INTEGER, 0, COMM_CART, ierr )
		
			
	! Pid 0 receiving data and putting into final file
	IF (pid .EQ. 0) THEN

		DO i = 0, Nprocs-1
		
			! Creating receive subarray type bespoke to each processor
			CALL MPI_Type_create_subarray(2, [ny, nx], [subarray_rows_array(i+1), subarray_cols_array(i+1)], &
			[subarray_row_start(i+1),subarray_col_start(i+1)], MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, RECVSUBBARAY, ierr)
			
			CALL MPI_TYPE_COMMIT(RECVSUBBARAY, ierr)
		
			! Receiving with this new receiving subarray type
			CALL MPI_IRECV( Tfinal, 1, RECVSUBBARAY, &
				i, tag2, COMM_CART, request_array(i+2), ierr)
		END DO

		CALL MPI_WAITALL(Nprocs+1, request_array, status_array, ierr)
		
		
		! --- PUTTING INTO FILE ---
		! x vector of whole domain (could have also done a mpi_gatherv)
		xtot 	= 0. + dx * [(i, i=0,(nx-1))] ! Implied DO loop used
		! y vector of whole domain (could have also done a mpi_gatherv)
		ytot 	= 1. - dy * [(i, i=0,(ny-1))] ! Implied DO loop used

		CALL tecplot_2D ( iunit, nx, ny, xtot, ytot, Tfinal )
		
		! Finally, the program has finshed.
		t2 = MPI_WTIME()
		WRITE(*,*) "Time elapsed for simulation: ", t2-t1, "seconds"
		
	END IF

	
	CALL MPI_FINALIZE(ierr) ! Terminate MPI Program
	
END PROGRAM Assignment1