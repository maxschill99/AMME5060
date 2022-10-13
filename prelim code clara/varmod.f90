! Contains all variables and the intialisation subroutine

MODULE variablemodule

	! Ensuring all names in the program must be explicitly declared
	IMPLICIT NONE
	
	! *************************************************
	!!! --- VARIABLE DECLARATION ---
	REAL (KIND = 8) :: Lx, Ly, dx, dy, alpha, c1, c2, c3, c4, c5, t1, t2, T_max, T_avg, res_max
	REAL (KIND = 8) :: pi
	REAL (KIND = 8) :: res_proc_val, res_final_val
	INTEGER :: nx, ny, iunit
	INTEGER :: ierr, pid_global, pid, Nprocs, i, j, tag, tag2, subarray_Nrows, subarray_Ncols, k, iter_max
	INTEGER :: bestdiv, bestgap, divisor, ndims, COMM_CART, left, right, north, south, y_flip_coord
	INTEGER, ALLOCATABLE :: dims(:), coords(:)
	LOGICAL, ALLOCATABLE :: periods(:)
	
	INTEGER :: node_low_x, node_high_x, node_low_y, node_high_y
	INTEGER :: ind_low_x, ind_high_x, ind_low_y, ind_high_y, ncalcpoints_x, ncalcpoints_y
	
	REAL (KIND = 8), ALLOCATABLE :: x(:), y(:), U(:,:), V(:,:), Tk(:,:), Tk_1(:,:), Res(:,:), Res_array(:), Tfinal(:,:)
	REAL (KIND = 8), ALLOCATABLE :: xtot(:), ytot(:)

	REAL (KIND = 8) :: SENDINGSUBBARAY, RECVSUBBARAY, RESIZED_RECVSUBBARAY, NS_ROW_SENDRECV
	INTEGER, ALLOCATABLE :: subarray_rows_array(:), subarray_cols_array(:), subarray_row_start(:), subarray_col_start(:)
	INTEGER, ALLOCATABLE :: status_array(:,:), request_array(:)
	
	
	! --------------------------------------------
	CONTAINS
	! --------------------------------------------

	SUBROUTINE intialise()
	
		! Universal constants
		pi = ACOS(-1.)
	
		!!! --- VARIABLE ASSIGNMENT ---
		!! - Problem Physical Parameters
		Lx		= 1. 				! x-domain length
		Ly		= 1. 				! y-domain length
		alpha 	= 80				! Thermal Conductivity [m^2/s]	&&&&&
	
		!! - File I/O
		iunit 	= 11 				! I/O unit number
	
		!! - Modelling-specific Parameters
		nx 		= 150				! Number of points in domain in x direction (evenly spaced)
		ny 		= 150				! Number of points in domain in y direction (evenly spaced)
		dx 		= Lx/(nx-1)			! Spatial step based on desired number of points [m]
		dy 		= Ly/(ny-1)			! Spatial step based on desired number of points [m]
		
		!! - Send/Recv -
		tag		= 1
		tag2	= 2
		
		!! - Iteration Loops -
		res_max = 0.001				! Maximum residual value
		iter_max= 200000 			! Maximum number of iterations allowed
		k 		= 1					! Starting loop counter
	
	END SUBROUTINE intialise
	

END MODULE variablemodule