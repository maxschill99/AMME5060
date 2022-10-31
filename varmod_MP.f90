! Contains all variables and the intialisation subroutine

MODULE variablemodule

	! Ensuring all names in the program must be explicitly declared
	IMPLICIT NONE
	
	! *************************************************
	!!! --- VARIABLE DECLARATION ---
	
	! Initialisation
	REAL (KIND = 8) :: pi, Lx, Ly, dx, dy, alpha, t_final, dt
	INTEGER :: nx, ny, iunit, Ntsteps
	
	! Problem variables / arrays
	CHARACTER(LEN=10) :: topology, solvertype
	REAL (KIND = 8), ALLOCATABLE :: x(:), y(:)
	
	! Partitioning and node/indices
	INTEGER :: p, n, r, s, k , m , d, e, alloc_pid, cumulative_neighbours, zone, zone_col, zone_row
	INTEGER :: node_low_x, node_low_y, node_high_x, node_high_y, gift_to_zone_A, gift_to_zone_B, rem, rem_2, nx_zA, nx_zB
	INTEGER :: ind_low_x, ind_high_x, ind_low_y, ind_high_y
		! Specifically for cartesian/slabs:
	INTEGER :: bestdiv, bestgap, divisor, ndims, COMM_CART
	INTEGER, ALLOCATABLE :: dims(:), coords(:)
	LOGICAL, ALLOCATABLE :: periods(:)
	
	! MPI Communicator
	INTEGER :: COMM_GRAPH, neighbours_count, Nedges, size_edges
	INTEGER, ALLOCATABLE :: indexes(:), edges(:), graphget_edges(:), graphget_indexes(:), neighbours_array(:)
	
	! Solving
	REAL (KIND = 8) :: res_proc_val, res_final_val
	INTEGER :: iter_max
	
	! MPI Send/Recv and Neighbour Information
	INTEGER :: north, south, east1, east2, west1, west2, znorth, zsouth, zeast1, zeast2, zwest1, zwest2
	INTEGER :: ierr, pid, Nprocs, i, j, tag
	INTEGER, ALLOCATABLE :: status_array(:,:), request_array(:), status_array_gather(:,:), request_array_gather(:)
	REAL (KIND = 8) :: NS_ROW_SENDRECV 
	INTEGER :: ind_low_east1, ind_low_east2, ind_low_west1, ind_low_west2
	INTEGER :: ind_high_east1, ind_high_east2, ind_high_west1, ind_high_west2, ignore_val
	INTEGER :: ncalcpoints_x, ncalcpoints_y, ncalcpoints_y_east1, ncalcpoints_y_east2, ncalcpoints_y_west1, ncalcpoints_y_west2

	! Initialising solver variables
    Real(kind = 8), allocatable :: T(:,:)
    Real(kind = 8), allocatable :: an(:,:), as(:,:), ae(:,:), aw(:,:), ap(:,:), b(:,:)
    Integer(kind = 8) :: il, ih, jl, jh, npp, iter, resil, resih, resjl, resjh
	Real(kind = 8) :: CFL

    ! Solution solver variables
    Real(kind = 8) :: rcurrent, rc, time
    Real(kind = 8), allocatable :: Told(:,:), Tn(:,:), resmat(:,:)

    ! Gathering variables for final solution
    Real(kind = 8), allocatable :: Ttemp(:,:), Ttot(:,:), Tinit(:,:), Tinittot(:,:)
    Real(kind = 8) :: numcount

	! Plotting
	Character(len=1024) ::  file_name

	! -------------------------------------------------------------------------------------------------------------------------------
	! `&&&` is used as a tag in comments to come back to something later. At end of assignment we can ctrlF for &&& and fix those things
	! &&& Variables not in code yet but maybe needed later, stashing them here for now:
	
	! Collecting information to one processor at end:
	INTEGER, ALLOCATABLE :: subarray_rows_array(:), subarray_cols_array(:), subarray_row_start(:), subarray_col_start(:)
	REAL (KIND = 8) :: SENDINGSUBBARAY, RECVSUBBARAY, RESIZED_RECVSUBBARAY
	INTEGER :: subarray_Nrows, subarray_Ncols, tag2
	
	
	! Arrays, solving constants, max/av temp, wall timing start/end
	REAL (KIND = 8) :: c1, c2, c3, c4, c5, t1, t2, T_max, T_avg, res_max
	REAL (KIND = 8), ALLOCATABLE :: Tk(:,:), Tk_1(:,:), Res(:,:), Res_array(:), Tfinal(:,:)
	REAL (KIND = 8), ALLOCATABLE :: xtot(:), ytot(:)
	
	
	! --------------------------------------------
	CONTAINS
	! --------------------------------------------

	SUBROUTINE intialise()
	
	! CHANGE WHAT'S IN HERE, MARK
	! ---------------------------------------------------------------------------------------------------
		nx 		= 15 !500				! Number of points in domain in x direction (evenly spaced)
		ny 		= 15 !500				! Number of points in domain in y direction (evenly spaced)
		res_max = 1e-7					! Maximum residual value
		dt		= t_final/Ntsteps 		! Time step size [s]
	! ---------------------------------------------------------------------------------------------------	
	
		! Universal constants
		pi = ACOS(-1.)
	
		!!! --- VARIABLE ASSIGNMENT ---
		!! - Problem Physical Parameters
		Lx		= 1. 				! x-domain length
		Ly		= 1.				! y-domain length
		alpha 	= 0.128E-6			! Pizza Crust Thermal Diffusivity [m^2/s] (from https://www.tandfonline.com/doi/pdf/10.1081/JFP-120015599)
		! alpha = 1

		!! - File I/O
		iunit 	= 11 				! I/O unit number
	
		!! - Modelling-specific Parameters
		dx 		= Lx/(nx-1)			! Spatial step based on desired number of points [m]
		dy 		= Ly/(ny-1)			! Spatial step based on desired number of points [m]
		Ntsteps = 1000				! Number of time steps
		t_final = 600.				! [s] &&& just a guess for now
		
		
		!! - Send/Recv -
		tag		= 1
		tag2	= 5
		
		!! - Iteration Loops -
		
		iter_max= 200000 			! Maximum number of iterations allowed
		k 		= 1					! Starting loop counter
	
	END SUBROUTINE intialise
	

END MODULE variablemodule