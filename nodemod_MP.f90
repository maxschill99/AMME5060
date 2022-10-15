! Determination of nodes and indices for each processor

MODULE nodemodule
	
	! Ensuring all names in the program must be explicitly declared
	IMPLICIT NONE
	
	! --------------------------------------------
	CONTAINS
	! --------------------------------------------
	
	SUBROUTINE get_nodes(Ntpoints, Nprocs, coord_id, node_low, node_high)
	
		IMPLICIT NONE
	
		INTEGER, INTENT(IN) :: Ntpoints, Nprocs, coord_id
		INTEGER, INTENT(OUT) :: node_low, node_high
		INTEGER :: Npoints_rem, Npoints
		
		! Number of points remaining when dividing number of rod points by number of processors
		Npoints_rem = MOD(Ntpoints, Nprocs)
	
		! If it works out nicely that each processor can work with the same number of points
		IF ( Npoints_rem .EQ. 0 ) THEN
			Npoints 	= 	Ntpoints/Nprocs		! Number of points considered for each processor
			! Domain decomposition - node numbers assigned to each processor depending on coord_id
			node_low 	= 	coord_id*Npoints + 1
			node_high 	= 	(coord_id+1)*Npoints
		
		! If it doesn't work out nicely, then everything should be divided as equally as possible.
		ELSE
			! Baseline number of points for each processor will be the nearest (low/floor) integer
			Npoints 		=  FLOOR( DBLE(Ntpoints)/DBLE(Nprocs) )
		
			! The extra points need to still be considered and distributed evenly:
			IF ( coord_id .LT. Npoints_rem ) THEN ! Not .LE. because coord_id #s start at 0
				! Giving the extra points to the first processors
				Npoints		= Npoints + 1
				! Node numbers assigned each processor depending on coord_id
				node_low 	= coord_id*Npoints + 1
				node_high 	= (coord_id+1)*Npoints
			ELSE
				! Node numbers assigned each processor depending on coord_id
				! Accounting for extra nodes handed out to the first processors
				node_low 	= coord_id*Npoints + 1 + Npoints_rem
				node_high 	= (coord_id+1)*Npoints + Npoints_rem
			END IF
		END IF
	
	END SUBROUTINE get_nodes
	
	
	! *******************************************************
	
	
	SUBROUTINE get_indices(Nprocs, coord_id, node_low, node_high, ind_low, ind_high, N_points_looped_thru)
	
		IMPLICIT NONE
	
		INTEGER, INTENT(IN) :: Nprocs, coord_id, node_low, node_high
		INTEGER, INTENT(OUT) :: ind_low, ind_high, N_points_looped_thru
		
		! Determining indices based on: Nodes considered AND extra spaces for ghost nodes
	
		IF ( (coord_id .EQ. 0) .AND. (Nprocs .EQ. 1) ) THEN 
			! No ghost nodes (only one processor is being used)
			ind_low 	= node_low
			ind_high	= node_high 
		
		ELSE IF ( (coord_id .EQ. 0) .AND. (Nprocs .GT. 1) ) THEN ! First processor on the left end of the rod
			! Ghost node is only on the right 
			ind_low 	= node_low
			ind_high	= node_high + 1
		
		ELSE IF ( coord_id .EQ. Nprocs-1)  THEN ! Last processor (Nprocs-1 since first processor # starts at zero)
			! Ghost node is only on the left 
			ind_low 	= node_low-1
			ind_high	= node_high
		
		ELSE ! All processors considering points in the middle of the rod
			! Ghost nodes are on either side
			ind_low 	= node_low-1
			ind_high	= node_high+1
		END IF
		
		! In the iterations, only the non-ghost nodes/ BCs will be looped through.
		! Calculating this now because it helps with the send/recv commands later.
		N_points_looped_thru = (ind_high-1) - (ind_low+1) + 1
	
	END SUBROUTINE get_indices
	
END MODULE nodemodule