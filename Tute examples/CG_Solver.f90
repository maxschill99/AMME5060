!-Conjugate gradient Solver Template!
PROGRAM CGSolver
	IMPLICIT NONE
	
	INTEGER :: kx,ky,il,ih,jl,jh
	INTEGER :: i,niter_precon
	
	Logical :: precon 
	
	REAL(KIND=8), PARAMETER :: pi = 3.14159265358979323846 
	REAL(KIND=8) :: time1,time2,timetotal
	REAL(KIND=8) :: meshlx,meshly,gam,rcurrent,rmax,dx,dy
	REAL(KIND=8), ALLOCATABLE :: aw(:,:),ae(:,:),ap(:,:),an(:,:),as(:,:),t(:,:),b(:,:),x(:),y(:)

	
	!-----------------------------------------
	! USER INPUT------------------------------
	kx = 100; ky = 100 					! MESHSIZE	
	meshlx = 1.0; meshly = 1.0 	! DOMAIN SIZE
	rmax = 0.0001 							! MAX RESIDUAL 	  
  gam = 100.0 								! HEAT FORCING TERM
	niter_precon = 5 						! NUMBER OF PRECON ITERATIONS
	precon = .TRUE.  						! USE A PRECONDITIONER
	
	
	
	!-----------------------------------------
	! Allocate memory
	CALL AllocateMemory()
	
	
	!-----------------------------------------
	! Initialize Solution Parameters
	WRITE(*,'(a)') "# Initializing Simulation - CG Solver for Steady-State Heat "	
	CALL SolInit()
	WRITE(*,'(a20,i14.1)') "# Nx = ", kx
	WRITE(*,'(a20,i14.1)') "# Ny = ", ky
	WRITE(*,'(a20,f14.1)') "# Lx = ", meshlx
	WRITE(*,'(a20,f14.1)') "# Ly = ", meshly
	WRITE(*,'(a20,E20.10)') "# Max Residual = ", rmax
	IF (precon) THEN
		WRITE(*,'(a,i10.1)') "# Number of Jacobi Precon Iterations =  ", niter_precon
  ENDIF
  WRITE(*,'(a)') "# Solution Initialized "
  
  
  !-----------------------------------------
  ! Solver Start
  
	! Start Timer
	CALL CPU_TIME(time1)
	
	! Solve Using CG method
	CALL CGSolve(aw,ae,an,as,ap,b,t,rmax)
	
	! Finnish Timer
  CALL CPU_TIME(time2)
  WRITE(*,'(a)') "# Simulation Finnished "
	WRITE(*,'(a15,f14.10)') "# Total WTime = ",  time2 - time1
	
	! Write out full Temperature Field to disk for visualization
	
	!-----------------------------
	! STUDENTS TO COMPLETE HERE---
	!-----------------------------
	
Contains

!-------------------------------------------------
! Subroutine to allocate memory
SUBROUTINE AllocateMemory()
	IMPLICIT NONE
	
	ALLOCATE(aw(kx,ky))
	ALLOCATE(ae(kx,ky))
	ALLOCATE(an(kx,ky))
    ALLOCATE(as(kx,ky))
    ALLOCATE(ap(kx,ky))
	ALLOCATE(t(kx,ky))
	ALLOCATE(b(kx,ky))
	ALLOCATE(x(kx))
	ALLOCATE(y(ky))

END SUBROUTINE AllocateMemory

!-------------------------------------------------
! Subroutine to Initialize Solution and establish Variables
SUBROUTINE SolInit()
	IMPLICIT NONE
	INTEGER :: i,j
	
	! Interior Indices
	il = 2; jl = 2;
	ih = kx-1; jh = ky-1
	
	! Grid Size
	dx = meshlx/REAL(kx-1)
	dy = meshly/REAL(ky-1)

  ! x coordinate
	DO i = 1,kx
		x(i) = (i-1)*dx											
	ENDDO
	
	! y coordinate
	DO j = 1,ky
		y(j) = (j-1)*dy											
	ENDDO

	! Initialize solution matrices
	t(:,:) = 0.0
	
	! A and b matrices
	DO j = 1,ky
		DO i = 1,kx
		b(i,j) = gam*x(i)*y(j)*sin(x(i)*pi/meshlx)*sin(y(j)*pi/meshly)  
		aw(i,j) = -1/dx**2										
		ae(i,j) = -1/dx**2		
		an(i,j) = -1/dy**2										
		as(i,j) = -1/dy**2								
		ap(i,j) =  2/dx**2 + 2/dy**2		 
		ENDDO
	ENDDO
	

END SUBROUTINE Solinit

!-------------------------------------------------
! Conjugent Gradient Solver with Jacobi Preconditioner
SUBROUTINE CGSolve(aw,ae,an,as,ap,b,t,rmax)
	IMPLICIT NONE
	REAL(KIND = 8), DIMENSION(kx,ky), INTENT(IN) :: aw,ae,an,as,ap,b
  REAL(KIND = 8), DIMENSION(kx,ky), INTENT(INOUT) :: t
  REAL(KIND = 8) :: rmax,rcurrent
  
  REAL(KIND = 8), DIMENSION(kx,ky) :: d,q,s,r,k  
  REAL(KIND = 8) :: delta,delta_o,alpha,beta,gamma
  
  INTEGER :: i,j,iter
  
  ! Initialize the solution constants and vectors
  delta = 1.0 ; delta_o = 1.0
  alpha = 1.0; beta = 1.0; gamma = 1.0
  d = 0.0; q(:,:) = 0.0; s(:,:) = 0.0
  
  ! Iteration Counter
  iter = 0
  
  ! Initialiaze 1/ap Matrix
  DO j = 1,ky
			DO i = 1,kx
			
				k(i,j) = 1/ap(i,j)
			
			ENDDO
  	ENDDO
  
  ! Initialize Precondtioned search vector
  IF (precon) then
  	DO j = jl,jh
			DO i = il,ih
			
				d(i,j) = r(i,j)*k(i,j)
			
			ENDDO
  	ENDDO
  ENDIF
  
  ! Get Initial Residual
  CALL GetResidual(aw,ae,an,as,ap,b,t,r)
  
  ! Calculate Domain averaged residual for stopping critterion
  rcurrent = SUM(SUM(ABS(r(il:ih,jl:jh)),1),1) / ((kx-2)*(ky-2))
  
	IF (precon) then
  	! Get Preconditoned Matrix 'd'  
  	CALL JacobiPrecon(aw,ae,an,as,ap,r,d,k)
  ELSE
    ! Or else search vector is the residual
  	d = r
  ENDIF
   	
  ! Get delta and delta_o  
  CALL DotProduct(d,r,delta)
  delta_o = delta
  

  ! Start Solution Loop
  DO WHILE (rcurrent.GT.rmax)		
    
  	! Get 'q' Matrix
  	!-----------------------------
		! STUDENTS TO COMPLETE HERE---
		!-----------------------------
  	
  	! Get 'gamma' = d^T*q
  	!-----------------------------
		! STUDENTS TO COMPLETE HERE---
		!-----------------------------
 
  	
  	! Get 'alpha'
  	!-----------------------------
		! STUDENTS TO COMPLETE HERE---
		!-----------------------------
  	
  	! Update 't' 
  	DO j = jl,jh
			DO i = il,ih 
				
				!-----------------------------
				! STUDENTS TO COMPLETE HERE  -
				!-----------------------------
				
			ENDDO
		ENDDO	
  
  	! Update Residual
  	IF (MOD(iter,50).eq.0) THEN
  		!-----------------------------
			! STUDENTS TO COMPLETE HERE  -
			!-----------------------------
  	ELSE
  		DO j = jl,jh
				DO i = il,ih 
					
				!-----------------------------
				! STUDENTS TO COMPLETE HERE  -
				!-----------------------------
				
				ENDDO
			ENDDO		
		ENDIF
		
		IF (precon) THEN
			! Apply a Preconditioner to get 's'
			CALL JacobiPrecon(aw,ae,an,as,ap,r,s,k)
		ELSE
			s = r
		ENDIF		
		
		! Update 'delta_o'
		delta_o = delta
		
		! Update 'delta'
		CALL DotProduct(r,s,delta)
		
		! Get Beta
		beta = delta/delta_o
		
		! Update 'd'
		DO j = jl,jh
			DO i = il,ih
			
				!-----------------------------
				! STUDENTS TO COMPLETE HERE  -
				!-----------------------------
			
			ENDDO
		ENDDO	
			
		! Update iteration counter and residual check	
		rcurrent = SUM(SUM(ABS(r(il:ih,jl:jh)),1),1) / ((kx-2)*(ky-2))	
  	iter = iter + 1
    
  	WRITE(*,'(a20,i20.1)') '# Current Iteration = ', iter
  	WRITE(*,'(a20,E20.6)') '# Current Residual = ', rcurrent
  	WRITE(*,'(a20,f20.10)') '# Current Max Temp = ', MAXVAL(t)

  
  ENDDO
  
END SUBROUTINE CGSOLVE

!-------------------------------------------------
! Subroutine to get Residual Matrix
SUBROUTINE GetResidual(aw,ae,an,as,ap,b,x,r)
	IMPLICIT NONE
	REAL(KIND = 8), DIMENSION(kx,ky), INTENT(IN) :: aw,ae,an,as,ap,b,x
	REAL(KIND = 8), DIMENSION(kx,ky), INTENT(OUT) :: r
	
	INTEGER :: i,j
	
	r(:,:) = 0.0
	
	DO j = jl,jh
		DO i = il,ih
		
				!-----------------------------
				! STUDENTS TO COMPLETE HERE  -
				!-----------------------------	
		
		ENDDO
	ENDDO		
END SUBROUTINE GetResidual


!-------------------------------------------------
! Subroutine to get Preconditoner Matrix 'd' using Jacobi
SUBROUTINE JacobiPrecon(aw,ae,an,as,ap,b,x,k)
	IMPLICIT NONE
	REAL(KIND = 8), DIMENSION(kx,ky), INTENT(IN) :: aw,ae,an,as,ap,b,k
	REAL(KIND = 8), DIMENSION(kx,ky), INTENT(INOUT) :: x

	REAL(KIND = 8), DIMENSION(kx,ky) :: r
	INTEGER :: i,j,ii
	
	DO j = jl,jh
		DO i = il,ih
			
			x(i,j) = b(i,j)*k(i,j)
			
		ENDDO
  ENDDO	

	DO ii = 1,niter_precon
    
		! Get Residual of Current system Ax = b
		CALL GetResidual(aw,ae,an,as,ap,b,x,r)
		
		! Update Solution
		DO j = jl,jh
			DO i = il,ih 
				
				x(i,j) = x(i,j) + r(i,j)*k(i,j)

			ENDDO
		ENDDO	
			
	ENDDO
	
END SUBROUTINE JacobiPrecon
	
!-------------------------------------------------
! Subroutine to get Dot Product of two solution matrices
SUBROUTINE DotProduct(l,m,dp)
	IMPLICIT NONE
	REAL(KIND = 8), DIMENSION(kx,ky), INTENT(IN) :: l,m
	REAL(KIND = 8), INTENT(OUT)	:: dp	
	INTEGER :: i,j
	
	dp = 0.0
	
		! Update Dot Product
		DO j = jl,jh
			DO i = il,ih
				
				!-----------------------------
				! STUDENTS TO COMPLETE HERE  -
				!-----------------------------
			
			ENDDO
		ENDDO	

END SUBROUTINE DotProduct
!-------------------------------------------------
! Subroutine to get Matrix Multiplication of two solution matrices
SUBROUTINE MatrixMultiply(aw,ae,an,as,ap,l,m)
	IMPLICIT NONE
	REAL(KIND = 8), DIMENSION(kx,ky), INTENT(IN) :: aw,ae,an,as,ap,l
	REAL(KIND = 8), DIMENSION(kx,ky), INTENT(OUT) :: m
	
	INTEGER :: i,j
		
		m = 0.0
		
		! Get m = A*l	
		DO j = jl,jh
			DO i = il,ih
			
				!-----------------------------
				! STUDENTS TO COMPLETE HERE  -
				!-----------------------------
			
			ENDDO
	  ENDDO		 
	

END SUBROUTINE MatrixMultiply

END PROGRAM CGSolver
