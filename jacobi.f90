module jacobi
    ! NOTE - need to call subroutines in this order --> allocatemmeory, solninit, solver
    ! NOTE - solution matrix should be initialise in main array
    ! SUBROUTINE GUIDE
        ! AllocateMemory() - allocates the memory size required for each variable
        ! SolnInit - Initialises the solution matrix
        ! Solver - solves the jacobi iteratively

    ! Calling modules
    USE variablemodule
    USE residuals
    ! USE cjgradient

    IMPLICIT NONE

    ! List of subroutines within module
    !-------------------------------------------------------------------------------------------!
    !-------------------------------------------------------------------------------------------!
    CONTAINS

    !-------------------------------------------------------------------------------------------!
    !-------------------------------------------------------------------------------------------!
    !-------------------------------------------------------------------------------------------!
    !-------------------------------------------------------------------------------------------!
    !-------------------------------------------------------------------------------------------!
    !-------------------------------------------------------------------------------------------!
    ! subroutine to allocate variables
    subroutine allocatevar()

        ! Real(kind = 8), INTENT(IN) :: an(:,:), as(:,:), ae(:,:), aw(:,:), ap(:,:), b(:,:)

        Real(kind = 8), allocatable :: an(:,:), as(:,:), ae(:,:), aw(:,:), ap(:,:), b(:,:)

        allocate(an(nx,ny))
        allocate(as(nx,ny))
        allocate(ae(nx,ny))
        allocate(aw(nx,ny))
        allocate(ap(nx,ny))
        allocate(b(nx,ny))
        ! allocate(x(nx))
        ! allocate(y(ny))
        ! allocate(Tin(nx,ny))
        ! allocate(Tout(nx,ny))


    end subroutine allocatevar

    !-------------------------------------------------------------------------------------------!
    !-------------------------------------------------------------------------------------------!
    !-------------------------------------------------------------------------------------------!
    !-------------------------------------------------------------------------------------------!
    !-------------------------------------------------------------------------------------------!
    !-------------------------------------------------------------------------------------------!
    ! subroutine to initialise solution
    SUBROUTINE solninit(an,as,ae,aw,ap,b,T)

        Real(kind = 8), INTENT(OUT) :: an(nx,ny), as(nx,ny), ae(nx,ny), aw(nx,ny), ap(nx,ny), b(nx,ny), T(nx,ny)

        ! Real(kind = 8) :: an(nx,ny), as(nx,ny), ae(nx,ny), aw(nx,ny), ap(nx,ny), b(nx,ny)
        ! Real(kind = 8) :: Tin(nx,ny), Tout(nx,ny)

        T(:,:) = 0

        ! Computing A and B matrices - need matrices for conjugate gradient method
        do j = 1,ny
            do i = 1,nx
                ! simplifying variable notation
                an(i,j) = (2*dt*alpha)/dy**2
                as(i,j) = (2*dt*alpha)/dy**2
                ae(i,j) = (2*dt*alpha)/dx**2
                aw(i,j) = (2*dt*alpha)/dx**2
                ap(i,j) = (-2/dx**2 - 2/dy**2)*2*dt*alpha
                b(i,j) = 0
            end do
        end do

        do i = 1,nx
            x(i) = (i-1)*dx
        end do
        do j = 1,ny
            y(j) = (j-1)*dy
        end do


        ! Setting solver boundary conditions

        do i = 1,nx
            T(1,i) = sin((pi*x(i))/Lx)
        end do
        T(:,1) = 0
        T(:,nx) = 0
        T(ny,:) = 0

    END SUBROUTINE solninit


    !-------------------------------------------------------------------------------------------!
    !-------------------------------------------------------------------------------------------!
    !-------------------------------------------------------------------------------------------!
    !-------------------------------------------------------------------------------------------!
    !-------------------------------------------------------------------------------------------!
    ! subroutine to solve jacobi
    SUBROUTINE jacobisolv(an,as,ae,aw,ap,b,x,y,T)

        IMPLICIT NONE
        ! Getting variables from varmod module
            ! nx - number of x points
            ! ny - number of y points
            ! dx - discretisation in x
            ! dy - discretisation in y
            ! dt - temporal discretisation
            ! dt - temporal discretisation
            ! alpha - thermal diffusivity
            ! t_final - total time

        ! Defining variables
        REAL(kind=8), INTENT(IN) :: an(nx,ny), as(nx,ny), ae(nx,ny), aw(nx,ny), ap(nx,ny), b(nx,ny), x(nx), y(ny)
        REAL(kind=8), INTENT(INOUT) :: T(nx,ny)

        Real(kind = 8) :: time, rcurrent
        Real(kind = 8) :: Tn(nx,ny), Told(nx,ny), res(nx,ny)
        Real(kind = 8) :: Minv(nx,ny)
        Integer :: i,j,iter, ii

        ! Set code case for type of solver, norm - for normal solving or conj - for conjugate gradient
        solving = "norm"

        !------------------------------------------------------------------------------------!
        !------------------------------------------------------------------------------------!

        ! Inititialising old temperature array
        Told(:,:) = 0

        ! Initialising time counter
        time = 0
        iter = 0

 
        ! Iteratively solving the jacobi equation
        ! Solving unsteady 2D Heat diffusion - dT/dt = alpha*(d^2T/dx^2 + d^2T/dy^2)

        SELECT CASE (solving)

            ! Jacobi Pre-conditioner
            CASE ("conj")

                do j = 1,ny
                ! do j = jl,jh
                    do i = 1,nx
                    ! do i = il,ih                       
                        Minv(i,j) = 1/ap(i,j)
                        T(i,j) = b(i,j)*Minv(i,j)                       
                    end do
                end do

                ! do ii = 1,niter_precon
                do ii = 1,5
                
                    ! Get Residual of Current system Ax = b
                    call residcalc(aw,ae,an,as,ap,b,T,res)
                    
                    ! Update Solution
                    do j = 1,ny
                    ! do j = jl,jh
                        do i = 1,nx
                        ! do i = il,ih               
                            T(i,j) = T(i,j) + res(i,j)*Minv(i,j)
                        end do
                    end do
                        
                end do


            !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!
            !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!
            !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!
            ! Jacobi solver
            CASE ("norm")
       
                ! do while ((t<ttot).and.(res>rmax))
                do while (time<t_final)
                    if (time == 0) then
                        do j = 2,(ny-1)
                        ! do j = jl,jh
                            do i = 2,(nx-1)
                            ! do i = il,ih
                                Tn(i,j) = (T(i+1,j)*ae(i,j) + T(i-1,j)*aw(i,j) + T(i,j+1)*an(i,j) &
                                    + T(i,j-1)*as(i,j) + T(i,j)*ap(i,j))/2
                            end do
                        end do
                    else
                        do j = 2,(ny-1)
                        ! do j = jl,jh
                            do i = 2,(nx-1)
                            ! do i = il,ih
                                Tn(i,j) = T(i+1,j)*ae(i,j) + T(i-1,j)*aw(i,j) + T(i,j+1)*an(i,j) &
                                    + T(i,j-1)*as(i,j) + T(i,j)*ap(i,j) - Told(i,j)
                            end do
                        end do
                    end if


                    ! computing residuals
                    call residcalc(aw,ae,an,as,ap,b,T,res)

                    ! Calculate Domain averaged residual for stopping critterion
                    ! rcurrent = SUM(SUM(ABS(r(il:ih,jl:jh)),1),1) / ((kx-2)*(ky-2))
                    rcurrent = SUM(SUM(ABS(res(1:nx,1:ny)),1),1) / ((nx-2)*(ny-2))

                    if (mod(iter,100).eq.0) then
                        write(*,*) '      iter', '      res'
                        write(*,*) iter, rcurrent
                    end if

                    ! updating old and new temperature values
                    Told = T
                    T = Tn

                    ! updating counter
                    time = time + dt
                    iter = iter + 1

                end do

            CASE DEFAULT 
                WRITE(*,*) "No solver selected or incorrect selection"
                STOP
	    END SELECT

        1600 FORMAT(5(F14.8,1x))

    END SUBROUTINE jacobisolv



    SUBROUTINE jacobiprecon(ae,aw,an,as,ap,b,Minv,T)

        IMPLICIT NONE
        Real(kind = 8), INTENT(IN) :: ae(nx,ny), aw(nx,ny), an(nx,ny), as(nx,ny), ap(nx,ny), b(nx,ny), Minv(nx,ny)
        Real(kind = 8), INTENT(INOUT) :: T(nx,ny)

        Real(kind = 8) :: res(nx,ny)
        Integer :: i,j,ii, niter_precon

        ! Setting number of preconditioning iterations
        niter_precon = 5
        
        do j = 1,ny
            do i = 1,nx
                T(i,j) = b(i,j)*Minv(i,j)
            end do
        end do	

        ! do ii = 1,5
        do ii = 1,niter_precon
        
            ! Get Residual of Current system Ax = b
            CALL residcalc(aw,ae,an,as,ap,b,T,res)
            
            ! Update Solution
            do j = 1,ny
                do i = 1,nx 
                    T(i,j) = T(i,j) + res(i,j)*Minv(i,j)
                end do
            end do
                
	    end do

    END SUBROUTINE jacobiprecon


end module jacobi