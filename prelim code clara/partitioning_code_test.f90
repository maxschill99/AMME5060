! Test code

PROGRAM test
	
	USE nodemodule
	
	IMPLICIT NONE
	
	INCLUDE "mpif.h"
	
	! Variable Declaration
    integer :: ierr, Nprocs, pid, nx, ny
	
	REAL (KIND = 8) :: Lx, Ly, dx, dy
	
	integer :: tag
	
	INTEGER, ALLOCATABLE :: status_array(:,:), request_array(:)

	REAL (KIND = 8) :: NS_ROW_SENDRECV

	integer :: i, p, n, r, s, k , m , d, e, alloc_pid, cumulative_neighbours, zone, zone_col, zone_row
	real (kind=8) :: r_col_factor, r_row_factor, s_col_factor, s_row_factor
	
	integer :: COMM_GRAPH, neighbours_count, Nedges, size_edges
	integer :: north, south, east1, east2, west1, west2, znorth, zsouth, zeast1, zeast2, zwest1, zwest2
	integer :: ncalcpoints_y_east1, ncalcpoints_y_east2, ncalcpoints_y_west1, ncalcpoints_y_west2
	integer :: ind_low_east1, ind_low_east2, ind_low_west1, ind_low_west2
	integer :: ind_high_east1, ind_high_east2, ind_high_west1, ind_high_west2, ignore_val
	
	INTEGER, ALLOCATABLE :: indexes(:), edges(:), graphget_edges(:), graphget_indexes(:), neighbours_array(:)
	INTEGER, ALLOCATABLE :: expected_edges(:), expected_indexes(:)
	
	integer :: node_low_x, node_low_y, node_high_x, node_high_y, gift_to_zone_A, gift_to_zone_B, rem, rem_2, nx_zA, nx_zB
	INTEGER :: ind_low_x, ind_high_x, ind_low_y, ind_high_y, ncalcpoints_x, ncalcpoints_y
	REAL (KIND = 8), ALLOCATABLE :: x(:), y(:), Tk_1(:,:)

	!!! --- VARIABLE ASSIGNMENT ---
	Lx		= 1. 				! x-domain length
	Ly		= 1. 				! y-domain length
		
	!! - Modelling-specific Parameters
	nx 		= 13 !2**3+1			! Number of points in domain in x direction (evenly spaced)
	ny 		= 13 !2**3+1			! Number of points in domain in y direction (evenly spaced)
	dx 		= Lx/(nx-1)			! Spatial step based on desired number of points [m]
	dy 		= Ly/(ny-1)			! Spatial step based on desired number of points [m]
	
	!! - Send/Recv -
	tag		= 1


	CALL MPI_INIT(ierr) ! Initialises everything 
	CALL MPI_COMM_RANK(MPI_COMM_WORLD, pid, ierr) ! Getting processor ID number
	CALL MPI_COMM_SIZE(MPI_COMM_WORLD, Nprocs, ierr) ! Getting number of total processors in global communicator
	
	! ----------------------------------------------------------
	! RUNNING THE OPTIMISED PARTITIONING OF A SQUARE ALGORITHM
	! ----------------------------------------------------------
	p = Nprocs
	n = FLOOR( SQRT(DBLE(Nprocs)) )
	
	IF ( p .LE. n*(n+1) ) THEN 
		r = n*(n+1) - p
		s = p - n**2
	ELSE
		r = (n+1)**2 - p
		s = p - n*(n+1)
	END IF

	
	! r rows of n rectangles. Side lengths
	r_row_factor = 1/DBLE(n)
	r_col_factor = DBLE(n)/DBLE(p)
	
	! s rows of n+1 rectangles
	s_row_factor = 1 /(DBLE(n) + 1)
	s_col_factor = ( DBLE(n)+1) / DBLE(p)
	
	if (pid == 0) then
		write(*,*) r, 'columns of ', n, 'rectangles with height (rows)', r_row_factor, 'width (cols)', r_col_factor
		write(*,*) s, 'columns of ', n+1, 'rectangles with height (rows)', s_row_factor, 'width (cols)', s_col_factor
	end if 

	
	! ***********************************************************
	! ************  G R A P H   T O P O L O G Y  ****************
	! ***********************************************************
	
	! Total number of neighbours. Sum of north south east west neighbours for all processors
	size_edges = (n-1)*r + (n-1)*r + (r-1)*n + n*r + n + n*s + n*s + (n+1)*s -1 + n + (n+1)*(s-1)
	ALLOCATE( edges(size_edges) )
	ALLOCATE( indexes(Nprocs) )
	
	! integers that increment in the loops
	alloc_pid = -1
	cumulative_neighbours = 0
	
	
	! Loop through zone A
	DO k = 1, r
		DO m = 1, n
		
			! Initialising all neighbour directions to nothing. 
			! This makes all the following `IF` statements tidier and helps later with general send/receives.
			! 2 east/west directions required because some processors will overlap with more than 1 neighbour (but never more than 2)
			znorth = MPI_PROC_NULL
			zsouth = MPI_PROC_NULL
			zeast1 = MPI_PROC_NULL
			zeast2 = MPI_PROC_NULL
			zwest1 = MPI_PROC_NULL
			zwest2 = MPI_PROC_NULL
	
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
				
			!   zwest2 = MPI_PROC_NULL
			
			IF (Nprocs .GT. 1) THEN
				zeast1 = alloc_pid + n
				cumulative_neighbours = cumulative_neighbours+1
				edges(cumulative_neighbours) = zeast1
			END IF
			IF ((Nprocs .GT. 1) .AND. (k .EQ. r)) THEN
				zeast2 = alloc_pid + n + 1
				cumulative_neighbours = cumulative_neighbours+1
				edges(cumulative_neighbours) = zeast2
			END IF
			
			! While we're here, each processor stores its own neighbours for later communication with them
			IF (pid .EQ. alloc_pid) THEN
				north = znorth
				south = zsouth
				east1 = zeast1
				east2 = zeast2
				west1 = zwest1
				west2 = zwest2
				zone  = 0
				zone_col = k-1 ! Start from zero
				zone_row = m-1 ! Start from zero
			END IF
			
			indexes(alloc_pid+1) = cumulative_neighbours
			
		END DO
	END DO
	
	
	! Loop through zone B
	DO d = 1, s
		DO e = 1, (n+1)
		
			! Initialising all neighbour directions to nothing
			znorth = MPI_PROC_NULL
			zsouth = MPI_PROC_NULL
			zeast1 = MPI_PROC_NULL
			zeast2 = MPI_PROC_NULL
			zwest1 = MPI_PROC_NULL
			zwest2 = MPI_PROC_NULL
			
			alloc_pid = alloc_pid + 1 ! Keep numbering going after zone A
						
			! Make the indexes array. Make the edges array. 
			
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
			IF ( (d.NE.1).OR.(e.NE.1) ) THEN	
				zwest1 = alloc_pid - (n+1)
				cumulative_neighbours = cumulative_neighbours+1
				edges(cumulative_neighbours) = zwest1
			END IF
			IF ( (d.EQ.1).AND.(e.LT. n+1) )	THEN
				zwest2 = alloc_pid - n
				cumulative_neighbours = cumulative_neighbours+1
				edges(cumulative_neighbours) = zwest2
			END IF
			IF (d .NE. s) THEN
				zeast1 = alloc_pid + n+1
				cumulative_neighbours = cumulative_neighbours+1
				edges(cumulative_neighbours) = zeast1
			END IF

			! 	zeast2 = MPI_PROC_NULL
				
			! While we're here, each processor stores its own neighbours for later communication with them
			IF (pid .EQ. alloc_pid) THEN
				north = znorth
				south = zsouth
				east1 = zeast1
				east2 = zeast2
				west1 = zwest1
				west2 = zwest2
				zone  = 1
				zone_col = d-1 ! Start from zero
				zone_row = e-1 ! Start from zero
			END IF
			
			indexes(alloc_pid+1) = cumulative_neighbours

		END DO
	END DO
	
	! ********************************************************************
	! ------------------- CONSIDERED NODE NUMBERS W/O GHOST NODES
	! ********************************************************************

	! ------------------- X (COLUMNS)
	rem = MOD(nx, p)
	
	IF (rem .EQ. 0) THEN
		! Number of nodes in zone A and zone B
		nx_zA = nx * r *     n/p 
		nx_zB = nx * s * (n+1)/p
	ELSE ! The number of points can't be nicely divided for each processor.
	
		! Find out what the extra points are, and give them as equally as possible to zones A and B
		rem_2 = MOD(rem, r+s)
		gift_to_zone_A = FLOOR( rem/( DBLE(r)+DBLE(s) ) )*r + rem_2
		gift_to_zone_B = FLOOR( rem/( DBLE(r)+DBLE(s) ) )*s
		
		! Number of points in x in zones A and B
		nx_zA = (nx-rem) * r *     n/p + gift_to_zone_A
		nx_zB = (nx-rem) * s * (n+1)/p + gift_to_zone_B
		
	END IF
	IF ( zone .EQ. 0) THEN 
		CALL get_nodes(nx_zA, r, zone_col, node_low_x, node_high_x)
	ELSE IF ( zone .EQ. 1) THEN
		CALL get_nodes(nx_zB, s, zone_col, node_low_x, node_high_x)
		! shift over since get_nodes expects points to start from 1, but zone B is to the right of zone A
		node_low_x  = node_low_x  + nx_zA
		node_high_x = node_high_x + nx_zA
	END IF
	
	! For `zone A` in y-direction there are `n` processors.
	! For `zone B` in y-direction there are `n+1` processors.
	

	! ------------------- Y (ROWS)
	IF ( zone .EQ. 0) THEN 
		CALL get_nodes(ny, n,   zone_row, node_low_y, node_high_y)
	ELSE IF ( zone .EQ. 1) THEN
		CALL get_nodes(ny, n+1, zone_row, node_low_y, node_high_y)
	END IF
	
	! ********************************************************************
	! ------------------- INDICES WITH GHOST NODES NOW
	! ********************************************************************
	! Each processor knows what nodes its considering, but now it needs to include the ghost nodes too.
	
	
	! CALL get_indices(Nprocs, zone, coord_id, node_low, node_high, ind_low, ind_high, N_points_looped_thru)
	
	! ------------------- Y (ROWS)
	IF ( zone .EQ. 0) THEN 
		! CALL get_nodes(ny, n,   zone_row, node_low_y, node_high_y)
		CALL get_indices(n, zone_row, node_low_y, node_high_y, ind_low_y, ind_high_y, ncalcpoints_y)
	ELSE IF ( zone .EQ. 1) THEN
		CALL get_indices(n+1, zone_row, node_low_y, node_high_y, ind_low_y, ind_high_y, ncalcpoints_y)
	END IF
	
	! ------------------- X (COLUMNS)
	! Include both zones when calculating indices in x so that whole domain span in x considered.
	IF ( zone .EQ. 0) THEN 
		CALL get_indices(r+s, zone_col, node_low_x, node_high_x, ind_low_x, ind_high_x, ncalcpoints_x)
	ELSE IF ( zone .EQ. 1) THEN
		! Changing zone_col coordinate id away from zero. Adding r to it instead. That way it recognises that zone B is in the middle.
		CALL get_indices(r+s, zone_col+r, node_low_x, node_high_x, ind_low_x, ind_high_x, ncalcpoints_x)
	END IF
	
	! ********************************************************************
	! ------------------- SEND/RECV INDICES IN EAST 1 DIRECTION
	! ********************************************************************
	ind_low_east1 = ind_low_y + 1
	IF (east2 .NE. -2) THEN ! east1 will not be a complete value. 
		! Sending y value based on y-span of the (shorter) east1 neighbour in zone B.
		! Could have done an MPI send/recv here for this value, or calculate it by itself
		! east1 neighbour has the same zone_row.
		CALL get_nodes(ny, n+1, zone_row, ignore_val, ind_high_east1)
	ELSE ! Everyone else lines up entirely with their east1 neighbour.
		ind_high_east1 = ind_high_y-1
	END IF
	ncalcpoints_y_east1 = ind_high_east1 - ind_low_east1 + 1
	
	! ********************************************************************
	! ------------------- SEND/RECV INDICES IN EAST 2 DIRECTION
	! ********************************************************************
	ind_low_east2 = ind_high_east1+1
	ind_high_east2 = ind_high_y-1
	ncalcpoints_y_east2 = ind_high_east2 - ind_low_east2 + 1
	
	! ********************************************************************
	! ------------------- SEND/RECV INDICES IN WEST 1 DIRECTION
	! ********************************************************************
	ind_low_west1 = ind_low_y + 1	
	IF ( (zone .EQ. 1) .AND. (zone_row .GT. 0) .AND. (zone_row .NE. n) ) THEN 
		! Sending y value based on y-span of the (longer) west1 neighbour in zone A.
		! Could have done an MPI send/recv here for this value, or calculate it by itself
		! west1 neighbour has zone_row-1.
		CALL get_nodes(ny, n, zone_row-1, ignore_val, ind_high_west1)
	ELSE ! Everyone else lines up entirely with their west1 neighbour.
		ind_high_west1 = ind_high_y-1
	END IF
	ncalcpoints_y_west1 = ind_high_west1 - ind_low_west1 + 1
	
	! ********************************************************************
	! ------------------- SEND/RECV INDICES IN WEST 2 DIRECTION
	! ********************************************************************
	IF ( (zone .EQ. 1) .AND. (zone_row .EQ. 0) ) THEN ! Special case for the only proc with no west1 neighbour
		ind_low_west2 = ind_low_y+1
	ELSE
		ind_low_west2 = ind_high_west1+1
	END IF
	ind_high_west2 = ind_high_y-1
	ncalcpoints_y_west2 = ind_high_west2 - ind_low_west2 + 1
	
	
	
	ALLOCATE( x(ind_low_x:ind_high_x), y(ind_low_y:ind_high_y) )

	! Defining x and y-coordinates considered by each processor (including ghost nodes)
	! In total, values are linearly spaced from 0 to Lx in increments of dx
	x 	= (ind_low_x-1)*dx + dx * [(i, i=0,(ind_high_x-ind_low_x))] ! Implied DO loop used
	! In total, values are linearly spaced from 0 to Ly in increments of dy
	! y vector created from 1 to 0 so that at y(1) (top of domain), the axis system is respected
	y 	= (ny-ind_low_y)*dy  -  dy * [(i, i=0,(ind_high_y-ind_low_y))] ! Implied DO loop used
	
	! if (pid == 0) then
		! write(*,*) 'points in x', nx
		! write(*,*) 'north', north, 'south', south, 'east1', east1, 'east2', east2, 'west1', west1, 'west2', west2
		! write(*,*) 'node_low_x', node_low_x, 'node_low_y', node_low_y
		! write(*,*) 'node_high_x', node_high_x, 'node_high_y', node_high_y
		! write(*,*) 'ind_low_x', ind_low_x, 'ind_low_y', ind_low_y
		! write(*,*) 'ind_high_x', ind_high_x, 'ind_high_y', ind_high_y
	! end if
	
	! -----------------------------------------------------------------------------------------------
	!	WHO KNOWS IF THIS SECTION IS NECESSARY
	! -----------------------------------------------------------------------------------------------

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
	
	
	! -----------------------------------------------------------------------------------------------
	! TESTING CODE WITH SEND/RECVS
	
	ALLOCATE( Tk_1(ind_low_y:ind_high_y, ind_low_x:ind_high_x) )
	Tk_1 = 0.
	Tk_1(node_low_y:node_high_y, node_low_x:node_high_x) = DBLE(pid)
	
	
	! ! ! ! --- SETTING UP SPECIFIC DATA TYPE ---
	! New vector type for sending and receiving the rows north and south (data not contiguous along a row)
	CALL MPI_TYPE_VECTOR(ncalcpoints_x, 1, ncalcpoints_y+2, MPI_DOUBLE_PRECISION, NS_ROW_SENDRECV , ierr)
	CALL MPI_TYPE_COMMIT(NS_ROW_SENDRECV, ierr)
	
	
	ALLOCATE( status_array(MPI_STATUS_SIZE, 12) )
	ALLOCATE( request_array(12) )
	
	! write(*,*) 'pid', pid, 'north', north, 'south', south, 'east1', east1, 'east2', east2, 'west1', west1, 'west2', west2
	! write(*,*) 'pid', pid, 'ILE1', ind_low_east1, 'IHE1',  ind_high_east1, 'ncalcpoints_y_east1', ncalcpoints_y_east1
	! write(*,*) 'pid', pid, 'ILE2', ind_low_east2, 'IHE2',  ind_high_east2, 'ncalcpoints_y_east2', ncalcpoints_y_east2
	! write(*,*) 'pid', pid, 'ILW1', ind_low_west1, 'IHW1',  ind_high_west1, 'ncalcpoints_y_west1', ncalcpoints_y_west1
	! write(*,*) 'pid', pid, 'ILW2', ind_low_west2, 'IHW2',  ind_high_west2, 'ncalcpoints_y_west2', ncalcpoints_y_west2

	
	! ------- RECEIVE FROM THE RIGHT/EAST, SEND TO THE RIGHT/EAST
	! Receive from the east1, and put it into the space for the ghost node on the east side
	CALL MPI_IRECV( Tk_1(ind_low_east1:ind_high_east1, ind_high_x), ncalcpoints_y_east1, MPI_DOUBLE_PRECISION, &
					east1, tag, COMM_GRAPH, request_array(1), ierr)
	! Send to the east1
	CALL MPI_ISEND( Tk_1(ind_low_east1:ind_high_east1, ind_high_x-1), ncalcpoints_y_east1, MPI_DOUBLE_PRECISION, &
					east1, tag, COMM_GRAPH, request_array(2), ierr)
					
	! Receive from the east2, and put it into the space for the ghost node on the east side
	CALL MPI_IRECV( Tk_1(ind_low_east2:ind_high_east2, ind_high_x), ncalcpoints_y_east2, MPI_DOUBLE_PRECISION, &
					east2, tag, COMM_GRAPH, request_array(3), ierr)
	! Send to the east2
	CALL MPI_ISEND( Tk_1(ind_low_east2:ind_high_east2, ind_high_x-1), ncalcpoints_y_east2, MPI_DOUBLE_PRECISION, &
					east2, tag, COMM_GRAPH, request_array(4), ierr)


	! ------- RECEIVE FROM THE LEFT/WEST, SEND TO THE LEFT/WEST
	! Receive from the west1, put it into the space for the ghost node on the west side 
	CALL MPI_IRECV( Tk_1(ind_low_west1:ind_high_west1, ind_low_x), ncalcpoints_y_west1, MPI_DOUBLE_PRECISION, & 
					west1, tag, COMM_GRAPH, request_array(5), ierr)
	! Send to the west1
	CALL MPI_ISEND( Tk_1(ind_low_west1:ind_high_west1, ind_low_x+1), ncalcpoints_y_west1, MPI_DOUBLE_PRECISION, &
					west1, tag, COMM_GRAPH, request_array(6), ierr)
					
	! Receive from the west2, put it into the space for the ghost node on the west side 
	CALL MPI_IRECV( Tk_1(ind_low_west2:ind_high_west2, ind_low_x), ncalcpoints_y_west2, MPI_DOUBLE_PRECISION, & 
					west2, tag, COMM_GRAPH, request_array(7), ierr)
	! Send to the west2
	CALL MPI_ISEND( Tk_1(ind_low_west2:ind_high_west2, ind_low_x+1), ncalcpoints_y_west2, MPI_DOUBLE_PRECISION, &
					west2, tag, COMM_GRAPH, request_array(8), ierr)
					
	! ------- RECEIVE FROM THE NORTH, SEND TO THE NORTH
	! Receive from the north, put it into the space for the ghost node 
	CALL MPI_IRECV( Tk_1(ind_low_y, ind_low_x+1), 1, NS_ROW_SENDRECV, &
					north, tag, COMM_GRAPH, request_array(9), ierr)
	! Send to the north
	CALL MPI_ISEND( Tk_1(ind_low_y+1, ind_low_x+1), 1, NS_ROW_SENDRECV, &
					north, tag, COMM_GRAPH, request_array(10), ierr)
	! ------- RECEIVE FROM THE SOUTH, SEND TO THE SOUTH
	! Receive from the south, put it into the space for the ghost node 
	CALL MPI_IRECV( Tk_1(ind_high_y, ind_low_x+1), 1, NS_ROW_SENDRECV, &
					south, tag, COMM_GRAPH, request_array(11), ierr)
	! Send to the south
	CALL MPI_ISEND( Tk_1(ind_high_y-1, ind_low_x+1), 1, NS_ROW_SENDRECV, &
					south, tag, COMM_GRAPH, request_array(12), ierr)

	! Wait for data sends to complete before black points start referencing red points
	CALL MPI_WAITALL(d, request_array, status_array, ierr)

	call sleep(2)
	if (pid == 3) then
		write(*,*) "pid 3 After send/recv"
		call printmatrix(Tk_1, SIZE(Tk_1, DIM = 1) , SIZE(Tk_1, DIM = 2))
	end if 


	CALL MPI_FINALIZE(ierr) ! Terminate MPI Program
	
END PROGRAM test

! b - matrix
subroutine printmatrix(b,n,m)
	integer::n,m
	real (kind=8)::b(n,m) !n = # rows, m = # columns
	do i=1,n; print '(20f6.2)',b(i,1:m); enddo
endsubroutine