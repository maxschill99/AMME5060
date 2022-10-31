! MAJOR PROJECT MAIN CODE OUTLINE
PROGRAM MAIN
! CALLING MODULES
	USE variablemodule 		! Contains variable allocation all problem-specific variable values
	USE nodemodule			! Obtains nodes and indices of the domain considered by each processor
	USE partitionmodule		! Contains subroutines to partition domain and obtain domain info for each processor
	USE outputmodule		! Creates tec file
	USE cjgradient
	USE jacobi
	USE residuals
	USE solinit

	include 'mpif.h'

! INITIALISE VARIABLES 
	! Variable declaration and assignment
	CALL intialise()

! -------------------------------------------------------------------------------	
!| CODE CASES										───▄▀▀▀▄▄▄▄▄▄▄▀▀▀▄───		 |
	topology = "slabs" 	!							───█▒▒░░░░░░░░░▒▒█───		 |
!|		Options: "graph" "cart" "slabs" 			────█░░█░░░░░█░░█────		 |
    solvertype = "jac"	!							─▄▄──█░░░▀█▀░░░█──▄▄─		 |
!| 		Options: "jac", "redblack", "conj"			█░░█─▀▄░░░░░░░▄▀─█░░█		 |
! -------------------------------------------------------------------------------

! 	INITIALISE MPI
	CALL MPI_INIT(ierr)
	CALL MPI_COMM_RANK(MPI_COMM_WORLD, pid, ierr) ! Getting processor ID number
	CALL MPI_COMM_SIZE(MPI_COMM_WORLD, Nprocs, ierr) ! Getting number of total processors in global communicator

! Begin timer
t1 = MPI_WTIME()

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
			
			CALL cart_nodesindices()
			
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
			
			CALL slab_nodesindices()
	
		CASE DEFAULT 
		  WRITE(*,*) "No topology selected or incorrect selection"
		  STOP
	END SELECT
	
	! New vector type for sending and receiving the rows north and south (data not contiguous along a row)
	CALL MPI_TYPE_VECTOR(ncalcpoints_x, 1, ncalcpoints_y+2, MPI_DOUBLE_PRECISION, NS_ROW_SENDRECV , ierr)
	CALL MPI_TYPE_COMMIT(NS_ROW_SENDRECV, ierr)
	

! ! !-----------------------------------------------------------------------------------------------------!
! ! !-----------------------------------------------------------------------------------------------------!
		! ┈┈╱▔▔▔▔▔▔▔▏ 
		! ┈╱ ╭▏╮╭┻┻╮╭┻┻╮╭▏ 
		! ▕╮╰▏╯┃╭╮┃┃╭╮┃╰▏ 
		! ▕╯┈▏┈┗┻┻┛┗┻┻┻╮▏ 
		! ▕╭╮▏╮┈┈┈┈┏━━━╯▏
		! ▕╰╯▏╯╰┳┳┳┳┳┳╯╭▏ 
		! ▕┈╭▏╭╮┃┗┛┗┛┃┈╰▏ 
		! ▕┈╰▏╰╯╰━━━━╯┈┈▏

	! ! !-----------------------------------------------------------------------------------------------------!
	! ! !-----------------------------------------------------------------------------------------------------!
	! INITIALISATION

	! Indices for computation
	il = ind_low_y; ih = ind_high_y; jl = ind_low_x; jh = ind_high_x ! Indices for temperature calculations
	resil = node_low_y; resih = node_high_y; resjl = node_low_x; resjh = node_high_x ! Nodes for residual calculations

	!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! INITIALISE TEMP DISTRIBUTION 
	! boundary conditions most pizza like
	! Allocation of variable sizes
	call allocatevars(an,as,ae,aw,ap,b,T,Told,Tn,res,il,ih,jl,jh)

	!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Initialising boundary conditions on temp array
	call solutioninit(an,as,ae,aw,ap,b,T,il,ih,jl,jh)

	! !!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Initialising new temp array and setting its value to T
	Tn(:,:) = 0
	Tn = T

	!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! INITIAL RESIDUAL
	! Initialising residual matrix
	allocate(resmat(resil:resih,resjl:resjh))
	resmat(:,:) = 0.0
	rc = 1

	!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! STABILITY CRITERION CHECK
	CFL = (1/(dx**2) + 1/(dy**2))*alpha*dt
	if (pid.eq.0) then
		write(*,*) 'Stability criterion: ', CFL
	end if
	if (CFL.GT.0.5) then
		write(*,*) 'Stability criterion not met, needs to be less than 0.5'
		STOP
	elseif (CFL.LT.0) then
		write(*,*) 'Stability criterion less than 0'
		STOP
	end if

	!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Initialising time counter and iteration counter
	time = 0
	iter = 0

	! ! !-----------------------------------------------------------------------------------------------------!
	! ! !-----------------------------------------------------------------------------------------------------!
	!!!!!!!!!!!!!!!!!!!!!!!!!!! COMPUTATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! ░▄▀▀▀▀▄░░▄▄
	! █░░░░░░▀▀░░█░░░░░░▄░▄
	! █░║░░░░██░████████████ 
	! █░░░░░░▄▄░░█░░░░░░▀░▀
	! ░▀▄▄▄▄▀░░▀▀


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
			
	! SOLVING

    SELECT CASE (solvertype)

        ! jacobi solver
        CASE ("jac")
            ! Begin the solution loop
            do while ((time<t_final).and.(rc>res_max))

				!-------------------------------------------------------------------!
                ! Calculation of solution using only jacobi solver				
                call jac(an,as,ae,aw,ap,b,T,il,ih,jl,jh)

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
                call respar(aw,ae,an,as,ap,b,T,resil,resih,resjl,resjh,resmat)

				! Calculate Domain averaged residual for stopping critterion
				rcurrent = SUM(SUM(ABS(resmat(resil:resih,resjl:resjh)),1),1) / ((resih-resil+1)*(resjh-resjl+1))

				! Summing processor residuals to get global resiudal and broadcasting to get average
				call MPI_ALLREDUCE(rcurrent,rc,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_TOPO,ierr)
				rc = rc/nprocs

				!-------------------------------------------------------------------!
				! Printing to screen after a certain amount of time
				if ((time-int(time))==0) then

					! TECPLOT

					! --------------------------------------------------------------------------------------- 
					!     _____     ____
					!    /      \  |  o | 
					!   |        |/ ___\|                TECPLOT
					!   |_________/     
					!   |_|_| |_|_|
					!	  
					! ---------------------------------------------------------------------------------------

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
						
					END IF
					
                end if

				!-------------------------------------------------------------------!
                ! updating counter
                time = time + dt
                iter = iter + 1

            end do


! 			 ;               ,           
!          ,;                 '.         
!         ;:                   :;        
!        ::                     ::       
!        ::                     ::       
!        ':                     :        
!         :.                    :        
!      ;' ::                   ::  '     
!     .'  ';                   ;'  '.    
!    ::    :;                 ;:    ::   
!    ;      :;.             ,;:     ::   
!    :;      :;:           ,;"      ::   
!    ::.      ':;  ..,.;  ;:'     ,.;:   
!     "'"...   '::,::::: ;:   .;.;""'    
!         '"""....;:::::;,;.;"""         
!     .:::.....'"':::::::'",...;::::;.   
!    ;:' '""'"";.,;:::::;.'""""""  ':;   
!   ::'         ;::;:::;::..         :;  
!  ::         ,;:::::::::::;:..       :: 
!  ;'     ,;;:;::::::::::::::;";..    ':.
! ::     ;:"  ::::::"""'::::::  ":     ::
!  :.    ::   ::::::;  :::::::   :     ; 
!   ;    ::   :::::::  :::::::   :    ;  
!    '   ::   ::::::....:::::'  ,:   '   
!     '  ::    :::::::::::::"   ::       
!        ::     ':::::::::"'    ::       
!        ':       """""""'      ::       
!         ::                   ;:        
!         ':;                 ;:"        
!           ';              ,;'          
!             "'           '"            
!               '

        CASE ("redblack")

            ! Begin the solution loop
            do while ((time<t_final).and.(rc>res_max))

				!-------------------------------------------------------------------!
                ! calculation of red nodes
                call rednodes(an,as,ae,aw,ap,b,T,il,ih,jl,jh)

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


				!-------------------------------------------------------------------!
                ! calculation of black nodes
                call blacknodes(an,as,ae,aw,ap,b,T,il,ih,jl,jh)


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
                call respar(aw,ae,an,as,ap,b,T,resil,resih,resjl,resjh,resmat)

				! Calculate Domain averaged residual for stopping critterion
				rcurrent = SUM(SUM(ABS(resmat(resil:resih,resjl:resjh)),1),1) / ((resih-resil+1)*(resjh-resjl+1))

				! Summing processor residuals to get global resiudal and broadcasting to get average
				call MPI_ALLREDUCE(rcurrent,rc,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_TOPO,ierr)
				rc = rc/nprocs

				!-------------------------------------------------------------------!
				! Printing to screen after a certain amount of time
				if ((time-int(time))==0) then
					! TECPLOT

					! --------------------------------------------------------------------------------------- 
					!     _____     ____
					!    /      \  |  o | 
					!   |        |/ ___\|                TECPLOT
					!   |_________/     
					!   |_|_| |_|_|
					!	  
					! ---------------------------------------------------------------------------------------

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
						write(file_name, "(A9,I5,A4)") "TecPlot2D",iter,".tec"
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

        CASE("Conj")

            ! ! Calculation of solution using Conjugate Gradient Method
            ! ! call CGSolve(an,as,ae,aw,ap,b,T)

    	CASE DEFAULT 
		WRITE(*,*) "No solver selected or incorrect selection"
		STOP

    END SELECT


	if (pid.eq.0) then
		write(*,*) 'Output after solver'
		! write(*,1600) Tfinal
		write(*,*) '-----------------------------------------'
		write(*,*) 'Time =', time
		write(*,*) 'Iteration =', iter
		write(*,*) 'Residual =', rc
		write(*,*) 'CFL =', CFL
		t2 = MPI_WTIME()
		write(*,*) 'Computational time =', t2-t1, 'seconds'
		write(*,*) 'Average Temperature =', sum(sum((Tfinal),1),1)/(nx*ny)
		write(*,*) '-----------------------------------------'

		! Writing the final iteration to a tecplot file
		write(file_name, "(A9,I5,A4)") "TecPlot2D",int(time),".tec"
		CALL tecplot_2D ( iunit, nx, ny, xtot, ytot, Tfinal,  file_name )

		! 	! Checking what size the A matrix coefficients are
		! 	write(*,*) 'A matrix coefficients'
		! 	! write(*,*) an(1,1), as(1,1), ae(1,1), aw(1,1), ap(1,1)
		! 	write(*,*) 'an =', an(1,1)
		! 	write(*,*) 'as =', as(1,1)
		! 	write(*,*) 'ae =', ae(1,1)
		! 	write(*,*) 'aw =', aw(1,1)
		! 	write(*,*) 'ap =', ap(1,1)
		! 	write(*,*) 'b =', b(1,1)
	end if
	! Outputting the domain indices
	! write(*,*) nx, ny
	! write(*,*) pid,il,ih,jl,jh
	! write(*,*) pid, resil,resih,resjl,resjh


	
	
	
	! *************************************************
	! --- Collating Final Temperature Arrays ---
	! T arrays from node values need to be sent to processor zero who needs to put them in the right place.
	! Sending must consider that ghost nodes need to be ignored, and that looks different for all processors
	! Receving will also need to consider where the matrix goes in the final domain (will not be contiguous)
	
	! ---- MAXIMUM AND AVERAGE TEMPERATURE ----
	! Printing maximum and average temperature to 8 decimal places:
	
	! ! Max temperature. Can take max over ghost nodes too because it doesn't matter - after global max.
	! CALL MPI_REDUCE(maxval(T), &
	! T_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, COMM_TOPO, ierr)
	
	! ! Average temperature. Only calculate average over nodes considered by each processor. Ignore ghost nodes.
	! CALL MPI_REDUCE( (1/(DBLE(nx)*DBLE(ny))) * sum( T(node_low_y:node_high_y, node_low_x:node_high_x) ), &
 	! T_avg, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, COMM_TOPO, ierr)
	
	! IF (pid .EQ. 0) THEN 
		! WRITE(*,'(A35, f10.8, A16)') "Maximum temperature across domain: ", T_max, " degrees Celcius"
		! WRITE(*,'(A35, f10.8, A16)') "Average temperature across domain: ", T_avg, " degrees Celcius"
	! END IF
	
	! ---- SEND DATA TYPE CREATION ----
	! The final send/recvs will be dependent on number of processors unlike before. Re-allocate relevant arrays to reflect this

CALL MPI_FINALIZE(ierr)

    ! NOTE: Need to update this to match the number of spatial divisions (nx)
    1100 FORMAT(8(F20.10,1x))
    1600 FORMAT(8(F14.8,1x))
    1400 FORMAT(8(F14.8,1x))
    1200 FORMAT(I2.1,I2.1,6(F8.4,1x))
    1300 FORMAT(I2.1,6(F10.8,1x))

END PROGRAM MAIN

! b - matrix
subroutine printmatrix(b,n,m)
	integer::n,m
	real (kind=8)::b(n,m) !n = # rows, m = # columns
	do i=1,n; print '(20f6.2)',b(i,1:m); enddo
endsubroutine