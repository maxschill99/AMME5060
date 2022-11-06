! Partitions Square Domain into rectangles of equal area and minimised perimeter

MODULE partitionmodule

	USE variablemodule
	USE nodemodule

	! Ensuring all names in the program must be explicitly declared
	IMPLICIT NONE
	
	! --------------------------------------------
	CONTAINS
	! --------------------------------------------

	SUBROUTINE graph_partition()
		IMPLICIT NONE
	
		! Constants p and n
		! p = number of rectangles the square is to be divided into
		! n = nearest square less than or equal to p
		p = Nprocs 
		n = FLOOR( SQRT(DBLE(Nprocs)) )
		
		! Finding r and s according to partitioning algorithm
		IF ( p .LE. n*(n+1) ) THEN 
			r = n*(n+1) - p
			s = p - n**2
		ELSE
			r = (n+1)**2 - p
			s = p - n*(n+1)
		END IF
		
		! Writing to screen how the domain will be partitioned 
		! if (pid == 0) then
			! write(*,*) r, 'columns of ', n, 'rectangles with height (rows)', 1/DBLE(n), 'width (cols)', DBLE(n)/DBLE(p)
			! write(*,*) s, 'columns of ', n+1, 'rectangles with height (rows)', 1 /(DBLE(n) + 1), 'width (cols)', ( DBLE(n)+1) / DBLE(p)
		! end if 

		! **********************************************************************
		! ******************  G R A P H   T O P O L O G Y  *********************
		! **********************************************************************
		
		! Total number of neighbours based on partitioning variables given above
		size_edges = (n-1)*r + (n-1)*r + (r-1)*n + n*r + n + n*s + n*s + (n+1)*s -1 + n + (n+1)*(s-1)
		
		! Allocating arrays that are required to define the graph topology communicator
		ALLOCATE( edges(size_edges) )
		ALLOCATE( indexes(Nprocs) )
		
		! integers that increment in the following loops
		alloc_pid = -1				! pid in consideration inside the loop
		cumulative_neighbours = 0	! growing sum of number of neighbours (needed to define the graph topology communicator later)
		
		! -----
		! LOOPING THROUGH ZONE A (r columns of n rectangles)
		DO k = 1, r
			DO m = 1, n
			
				! Initialising all neighbour directions to nothing. 
				! This makes all the following `IF` statements tidier and helps later with general send/receives.
				! 2 east/west directions required because some processors will overlap with more than 1 neighbour (but never more than 2)
				znorth = -2 ! MPI_PROC_NULL
				zsouth = -2 ! MPI_PROC_NULL
				zeast1 = -2 ! MPI_PROC_NULL
				zeast2 = -2 ! MPI_PROC_NULL
				zwest1 = -2 ! MPI_PROC_NULL
				zwest2 = -2 ! MPI_PROC_NULL
		
				alloc_pid = alloc_pid + 1 ! Starts at zero.
			
				IF (m .NE. 1) THEN
					znorth = alloc_pid - 1
					cumulative_neighbours = cumulative_neighbours+1
					edges(cumulative_neighbours) = znorth
				END IF 
				IF (m .NE. n) THEN
					zsouth = alloc_pid + 1
					cumulative_neighbours = cumulative_neighbours+1
					edges(cumulative_neighbours) = zsouth
				END IF
				IF (k .NE. 1) THEN
					zwest1 = alloc_pid - n 
					cumulative_neighbours = cumulative_neighbours+1
					edges(cumulative_neighbours) = zwest1
				END IF
					
				!   (zwest2 = MPI_PROC_NULL)
				
				IF ((Nprocs .GT. 1) .AND. (k .NE. s+r)) THEN
					zeast1 = alloc_pid + n
					cumulative_neighbours = cumulative_neighbours+1
					edges(cumulative_neighbours) = zeast1
				END IF
				IF ((Nprocs .GT. 1) .AND. (k .EQ. r) .AND. (s .NE. 0)) THEN
					zeast2 = alloc_pid + n + 1
					cumulative_neighbours = cumulative_neighbours+1
					edges(cumulative_neighbours) = zeast2
				END IF
				
				! If the pid reading matches the alloc_pid in consideration, it stores its own neighbours for later communication with them
				IF (pid .EQ. alloc_pid) THEN
					north = znorth
					south = zsouth
					east1 = zeast1
					east2 = zeast2
					west1 = zwest1
					west2 = zwest2
					zone  = 0 		! 0 = zone A
					zone_col = k-1 	! Starting from zero
					zone_row = m-1 	! Starting from zero
				END IF
				
				! Adding to indexes array needed to define the graph topology communicator later
				indexes(alloc_pid+1) = cumulative_neighbours
				
				
			END DO
		END DO
		
		! -----
		! LOOPING THROUGH ZONE B (s columns of n+1 rectangles)
		DO d = 1, s
			DO e = 1, (n+1)
			
				! Initialising all neighbour directions to nothing
				znorth = -2 ! MPI_PROC_NULL
				zsouth = -2 ! MPI_PROC_NULL
				zeast1 = -2 ! MPI_PROC_NULL
				zeast2 = -2 ! MPI_PROC_NULL
				zwest1 = -2 ! MPI_PROC_NULL
				zwest2 = -2 ! MPI_PROC_NULL
				
				alloc_pid = alloc_pid + 1 ! Keep numbering going after zone A!
				
				IF (e .NE. 1) THEN	
					znorth = alloc_pid - 1
					cumulative_neighbours = cumulative_neighbours+1
					edges(cumulative_neighbours) = znorth
				END IF
				IF (e .NE. n+1) THEN 
					zsouth = alloc_pid + 1
					cumulative_neighbours = cumulative_neighbours+1
					edges(cumulative_neighbours) = zsouth
				END IF
				IF ( (d.NE.1).OR.(e.NE.1).AND.(r .NE. 0) ) THEN	
					zwest1 = alloc_pid - (n+1)
					cumulative_neighbours = cumulative_neighbours+1
					edges(cumulative_neighbours) = zwest1
				END IF
				IF ( (d.EQ.1).AND.(e.LT. n+1).AND.(r .NE. 0))	THEN
					zwest2 = alloc_pid - n
					cumulative_neighbours = cumulative_neighbours+1
					edges(cumulative_neighbours) = zwest2
				END IF
				IF (d .NE. s) THEN
					zeast1 = alloc_pid + n+1
					cumulative_neighbours = cumulative_neighbours+1
					edges(cumulative_neighbours) = zeast1
				END IF

				! 	(zeast2 = MPI_PROC_NULL)
					
				! If the pid reading matches the alloc_pid in consideration, it stores its own neighbours for later communication with them
				IF (pid .EQ. alloc_pid) THEN
					north = znorth
					south = zsouth
					east1 = zeast1
					east2 = zeast2
					west1 = zwest1
					west2 = zwest2
					zone  = 1		! 1 = zone B
					zone_col = d-1 	! Starting from zero
					zone_row = e-1 	! Starting from zero
				END IF
				
				! Adding to indexes array needed to define the graph topology communicator later
				indexes(alloc_pid+1) = cumulative_neighbours
				
				! IF (pid .EQ. alloc_pid) THEN
					! WRITE(*,*) "i am pid", pid, north, south, east1, east2, west1, west2
				! end if

			END DO
		END DO
		
		! **********************************************************************
		! ************  CONSIDERED NODE NUMBERS W/O GHOST NODE  ****************
		! **********************************************************************

		! ------------------- X (COLUMNS)
		
		! In both zone A and zone B, column widths are multiples of 1/p according to the partitioning algorithm
		! Checking remainder when dividing number of points in x by p
		rem = MOD(nx, p)
		
		IF (rem .EQ. 0) THEN
			! Everything works out nicely and the points in x can be divided perfectly between the two zones
			nx_zA = nx * r *     n/p ! Number of nodes in zone A
			nx_zB = nx * s * (n+1)/p ! Number of nodes in zone B
			
		ELSE ! The number of points can't be nicely divided for each column.
				
			! Of the leftover points found by `rem`, these should be distributed amongst the columns as evenly as possible
			rem_2 = MOD(rem, r+s) ! remainder dividing spare points by number of columns across whole domain (r columns in zone A, s columns in zone B)
			
			! Extra points given to each zone, with any final spare points just given to zone A
			IF (r .EQ. 0) then ! (Zone A doesnt exist, give extra point to zone B)
				gift_to_zone_A = FLOOR( rem/( DBLE(r)+DBLE(s) ) )*r 
				gift_to_zone_B = FLOOR( rem/( DBLE(r)+DBLE(s) ) )*s + rem_2
			ELSE
				gift_to_zone_A = FLOOR( rem/( DBLE(r)+DBLE(s) ) )*r + rem_2
				gift_to_zone_B = FLOOR( rem/( DBLE(r)+DBLE(s) ) )*s
			END IF
			
			! Number of points in x in zones A and B
			nx_zA = (nx-rem) * r *     n/p + gift_to_zone_A
			nx_zB = (nx-rem) * s * (n+1)/p + gift_to_zone_B
			
		END IF
		
		
		! With the number of points in x for both zones known, now each processor can use their knowledge of what column they are in
		! to figure out what nodes in the domain they will be considering
		IF ( zone .EQ. 0) THEN 		! Zone A
			CALL get_nodes(nx_zA, r, zone_col, node_low_x, node_high_x)
		ELSE IF ( zone .EQ. 1) THEN ! Zone B
			CALL get_nodes(nx_zB, s, zone_col, node_low_x, node_high_x)
			! shift over since get_nodes expects points to start from 1, but zone B is to the right of zone A
			node_low_x  = node_low_x  + nx_zA
			node_high_x = node_high_x + nx_zA
		END IF
		
		! For clarity: each pid comes out of the `get_nodes` function above with different node numbers because of its `zone_col` value

		! ------------------- Y (ROWS)
		! For `zone A` in y-direction there are `n` processors.
		! For `zone B` in y-direction there are `n+1` processors.
		IF ( zone .EQ. 0) THEN 
			CALL get_nodes(ny, n,   zone_row, node_low_y, node_high_y)
		ELSE IF ( zone .EQ. 1) THEN
			CALL get_nodes(ny, n+1, zone_row, node_low_y, node_high_y)
		END IF
		
		
		! **********************************************************************
		! *****************  INDICES WITH GHOST NODES NOW  *********************
		! **********************************************************************
		! Each processor knows what nodes its considering, but now it needs to include the ghost nodes too.
		! The indices obtained below define the temperature array for each pid
		
		! ------------------- X (COLUMNS)
		! Include both zones when calculating indices in x (total columns = r+s) so that whole domain span in x considered.
		IF ( zone .EQ. 0) THEN 
			CALL get_indices(r+s, zone_col, node_low_x, node_high_x, ind_low_x, ind_high_x, ncalcpoints_x)
		ELSE IF ( zone .EQ. 1) THEN
			! Changing zone_col coordinate id away from zero. Adding r to it instead. 
			! That way `get_indices` recognises that the LHS of zone B is in the middle, not at a boundary. 
			CALL get_indices(r+s, zone_col+r, node_low_x, node_high_x, ind_low_x, ind_high_x, ncalcpoints_x)
		END IF
		
		! ------------------- Y (ROWS)
		IF ( zone .EQ. 0) THEN 
			CALL get_indices(n, zone_row, node_low_y, node_high_y, ind_low_y, ind_high_y, ncalcpoints_y)
		ELSE IF ( zone .EQ. 1) THEN
			CALL get_indices(n+1, zone_row, node_low_y, node_high_y, ind_low_y, ind_high_y, ncalcpoints_y)
		END IF


		! ********************************************************************
		! ****** INDICES NEEDED TO SEND/RECV IN EAST/WEST DIRECTIONS *********
		! ********************************************************************
		
		! To send to the east and west, some processors will have two neighbours in one direction.
		! All processors (regardless of with/without 2 neighbours in east or west) need to know what portion of their
		! temperature array to send to their east/west neighbour(s), i.e. the indices corresponding to that data.
		
		! Note: All processors go through the following steps regardless of whether they have each type of neighbour.
		! Even though this may result in nonsensical data, it doesn't matter because they have MPI_PROC_NULL target rank in these sends and recvs
		
		! ------------------- EAST 1 DIRECTION
		ind_low_east1 = ind_low_y + 1
		IF (east2 .NE. -2) THEN ! east1 will not be a complete value. 
			! Sending y value based on y-span of the (shorter) east1 neighbour in zone B.
			! For this, we need to call the `get_nodes` function and "pretend to be the east 1 neighbour" to obtain
			! their index information
			CALL get_nodes(ny, n+1, zone_row, ignore_val, ind_high_east1) ! east1 neighbour has the same zone_row.
		ELSE ! Everyone else lines up entirely with their east1 neighbour.
			ind_high_east1 = ind_high_y-1
		END IF
		ncalcpoints_y_east1 = ind_high_east1 - ind_low_east1 + 1 ! For send/recv `count`
		
		! ------------------- EAST 2 DIRECTION
		ind_low_east2 = ind_high_east1+1 ! picks up just after where east1 index left off
		ind_high_east2 = ind_high_y-1
		ncalcpoints_y_east2 = ind_high_east2 - ind_low_east2 + 1! For send/recv `count`
		
		! ------------------- WEST 1 DIRECTION
		ind_low_west1 = ind_low_y + 1	
		IF ( (zone .EQ. 1) .AND. (zone_row .GT. 0) .AND. (zone_row .NE. n) .AND. (zone_col .EQ. 0) .AND. (r .NE. 0)) THEN 
			! Sending y value based on y-span of the (longer) west1 neighbour in zone A.
			! west1 neighbour has zone_row-1.
			CALL get_nodes(ny, n, zone_row-1, ignore_val, ind_high_west1)
		ELSE ! Everyone else lines up entirely with their west1 neighbour.
			ind_high_west1 = ind_high_y-1
		END IF
		ncalcpoints_y_west1 = ind_high_west1 - ind_low_west1 + 1 ! For send/recv `count`
		
		! ------------------- WEST 2 DIRECTION
		IF ( (zone .EQ. 1) .AND. (zone_row .EQ. 0) ) THEN ! Special case for the only proc with no west1 neighbour
			ind_low_west2 = ind_low_y+1
		ELSE
			ind_low_west2 = ind_high_west1+1
		END IF
		ind_high_west2 = ind_high_y-1
		ncalcpoints_y_west2 = ind_high_west2 - ind_low_west2 + 1 ! For send/recv `count`
	
	END SUBROUTINE graph_partition
	
	! ---------------------------------------------------------------------------------------
	!	   A_A
	!	  (-.-)
	!	   |-|
	!	  /   \						cat-esian decomposition subroutine
	!	 |     |   __
	!	 |  || |  |  \__
	!	  \_||_/_/
	! ---------------------------------------------------------------------------------------
	
	SUBROUTINE cart_partition()
		IMPLICIT NONE
		
		! The domain is split up in 2 dimensions if possible.
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
		
		! USE SLABS IF YOU WANT TO TEST PRIME NUMBERS OF PROCESSORS
		IF (bestdiv .EQ. 1) THEN ! The number is prime and the domain cannot be partitioned in 2 dimensions
			WRITE(*,*) 'ERROR: You have picked cartesian 2D but the number of processors is prime. Program stopping.'
			STOP
		ELSE
			dims(2) = Nprocs/bestdiv ! number of processors in y direction
			dims(1) = bestdiv ! This number will always be greater than or equal to dim(2) (more columns). Number of procs in x direction
		END IF
		
		! The domain is not periodic
		periods = .false.
		
	
	END SUBROUTINE cart_partition
	
	! ----------------------
	
	SUBROUTINE cart_nodesindices()
	
		IMPLICIT NONE
		
		east2 = -2 ! MPI_PROC_NULL
		west2 = -2 ! MPI_PROC_NULL
			
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
	
	END SUBROUTINE cart_nodesindices
	
	! --------------------------------------------------------------------------------------- 
	!     _____     ____
	!    /      \  |  o | 
	!   |        |/ ___\|                slab decomposition subroutine
	!   |_________/     
	!   |_|_| |_|_|
	!	  
	! ---------------------------------------------------------------------------------------
	
	
	SUBROUTINE slab_partition()
		IMPLICIT NONE
		
		! The domain is split up in 1 dimension
		ndims = 2 
		ALLOCATE( dims(ndims), periods(ndims), coords(ndims))
		
		dims(1) = Nprocs ! Number of procs in this one dimension is simply all of them (number of procs in x direction)
		dims(2) = 1 ! Number of procs in y direction is 1
		
		! The domain is not periodic
		periods = .false.
		
	END SUBROUTINE slab_partition
	
	SUBROUTINE slab_nodesindices()
	
		IMPLICIT NONE
		
		east2 = -2 ! MPI_PROC_NULL
		west2 = -2 ! MPI_PROC_NULL
			
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
	
	END SUBROUTINE slab_nodesindices

END MODULE partitionmodule