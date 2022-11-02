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
    ! parallel jacobi solver
    SUBROUTINE jac(an,as,ae,aw,ap,b,T,il,ih,jl,jh)


        ! Defining variables
        Real(kind=8), INTENT(IN), dimension(il:ih,jl:jh) :: an,as,ae,aw,ap,b
        Integer(kind = 8), INTENT(IN) :: il,ih,jl,jh
        Real(kind=8), INTENT(INOUT) :: T(il:ih,jl:jh)

        Integer :: i,j
        Real(kind = 8) :: C, F, hsq, k
        Real(kind = 8), dimension(il:ih,jl:jh) :: bmat

        C = ((dt*alpha)/dx**2)
        F = 1 + 4*C
        hsq = dx**2
        k = dt

        do j = jl+1,jh-1
            do i = il+1,ih-1
                ! EXPLICIT
                ! UNSTEADY
                ! Tn(i,j) = T(i+1,j)*an(i,j) + T(i-1,j)*as(i,j) + T(i,j+1)*ae(i,j) &
                !     + T(i,j-1)*aw(i,j) + T(i,j)*ap(i,j)
                ! Tn(i,j) = T(i,j)*ap(i,j) + alpha*dt*((T(i,j+1)+T(i,j-1))/(dx**2)) + alpha*dt*((T(i+1,j)+T(i-1,j))/(dy**2))  

                ! STEADY
                ! Tn(i,j) = (alpha*((T(i,j+1)+T(i,j-1))/(dx**2)) & 
                !             + alpha*((T(i+1,j)+T(i-1,j))/(dy**2)))/((2*alpha)/dx**2 + (2*alpha)/dy**2)    

                ! IMPLICIT
                ! JACOBI
                bmat(i,j) = T(i,j)/(1 + 4*C)
                Tn(i,j) = bmat(i,j) + C*(Tn(i+1,j) + Tn(i-1,j) + Tn(i,j+1) + Tn(i,j-1))/(1 + 4*C)

                ! CRANK NICHOLSON
                ! bmat(i,j) = (1 - (2*alpha*k)/hsq)*T(i,j) + (alpha*k)/(2*hsq)*(T(i+1,j) + T(i-1,j) + T(i,j+1) + T(i,j-1))
                ! Tn(i,j) = (bmat(i,j) + (alpha*k)/(2*hsq)*(Tn(i+1,j) + Tn(i-1,j) + Tn(i,j+1) + Tn(i,j-1)))/(1 + (2*alpha*k)/hsq)

            end do
        end do

        T = Tn

    END SUBROUTINE jac

    !-------------------------------------------------------------------------------------------!
    !-------------------------------------------------------------------------------------------!
    !-------------------------------------------------------------------------------------------!
    !-------------------------------------------------------------------------------------------!
    !-------------------------------------------------------------------------------------------!
    ! parallel redblack solver
    ! red nodes
    SUBROUTINE rednodes(an,as,ae,aw,ap,b,T,il,ih,jl,jh)

        ! Defining variables
        Real(kind=8), INTENT(IN), dimension(il:ih,jl:jh) :: an,as,ae,aw,ap,b
        Integer(kind = 8), INTENT(IN) :: il,ih,jl,jh
        Real(kind=8), INTENT(INOUT) :: T(il:ih,jl:jh)

        Integer :: i,j
        Real(kind = 8) :: C, F, hsq, k
        Real(kind = 8), dimension(il:ih,jl:jh) :: bmat

        C = ((dt*alpha)/dx**2)
        F = 1 + 4*C
        hsq = dx**2
        k = dt

        do j = jl+1,jh-1
            do i = il+1 + mod(j,2),ih-1,2
                ! EXPLICIT
                ! UNSTEADY
                Tn(i,j) = T(i+1,j)*an(i,j) + T(i-1,j)*as(i,j) + T(i,j+1)*ae(i,j) &
                    + T(i,j-1)*aw(i,j) + T(i,j)*ap(i,j)

                ! ! STEADY
                ! Tn(i,j) = (alpha*((T(i,j+1)+T(i,j-1))/(dx**2)) & 
                !             + alpha*((T(i+1,j)+T(i-1,j))/(dy**2)))/((2*alpha)/dx**2 + (2*alpha)/dy**2) 

                ! IMPLICIT
                ! JACOBI
                bmat(i,j) = T(i,j)/(1 + 4*C)
                Tn(i,j) = bmat(i,j) + C*(Tn(i+1,j) + Tn(i-1,j) + Tn(i,j+1) + Tn(i,j-1))/(1 + 4*C)

                ! CRANK NICHOLSON
                ! bmat(i,j) = (1 - (2*alpha*k)/hsq)*T(i,j) + (alpha*k)/(2*hsq)*(T(i+1,j) + T(i-1,j) + T(i,j+1) + T(i,j-1))
                ! Tn(i,j) = (bmat(i,j) + (alpha*k)/(2*hsq)*(Tn(i+1,j) + Tn(i-1,j) + Tn(i,j+1) + Tn(i,j-1)))/(1 + (2*alpha*k)/hsq)

            end do
        end do

        T = Tn

    END SUBROUTINE rednodes

    ! black nodes
    SUBROUTINE blacknodes(an,as,ae,aw,ap,b,T,il,ih,jl,jh)

        ! Defining variables
        Real(kind=8), INTENT(IN), dimension(il:ih,jl:jh) :: an,as,ae,aw,ap,b
        Integer(kind = 8), INTENT(IN) :: il,ih,jl,jh
        Real(kind=8), INTENT(INOUT) :: T(il:ih,jl:jh)

        Integer :: i,j
        Real(kind = 8) :: C, F, hsq, k
        Real(kind = 8), dimension(il:ih,jl:jh) :: bmat

        C = ((dt*alpha)/dx**2)
        F = 1 + 4*C
        hsq = dx**2
        k = dt

        do j = jl+1,jh-1
            do i = il+2 - mod(j,2),ih-1,2
                ! EXPLICIT
                ! UNSTEADY
                Tn(i,j) = T(i+1,j)*an(i,j) + T(i-1,j)*as(i,j) + T(i,j+1)*ae(i,j) &
                    + T(i,j-1)*aw(i,j) + T(i,j)*ap(i,j)

                ! ! STEADY
                ! Tn(i,j) = (alpha*((T(i,j+1)+T(i,j-1))/(dx**2)) & 
                !             + alpha*((T(i+1,j)+T(i-1,j))/(dy**2)))/((2*alpha)/dx**2 + (2*alpha)/dy**2) 

                ! IMPLICIT
                ! JACOBI
                bmat(i,j) = T(i,j)/(1 + 4*C)
                Tn(i,j) = bmat(i,j) + C*(Tn(i+1,j) + Tn(i-1,j) + Tn(i,j+1) + Tn(i,j-1))/(1 + 4*C)

                ! CRANK NICHOLSON
                ! bmat(i,j) = (1 - (2*alpha*k)/hsq)*T(i,j) + (alpha*k)/(2*hsq)*(T(i+1,j) + T(i-1,j) + T(i,j+1) + T(i,j-1))
                ! Tn(i,j) = (bmat(i,j) + (alpha*k)/(2*hsq)*(Tn(i+1,j) + Tn(i-1,j) + Tn(i,j+1) + Tn(i,j-1)))/(1 + (2*alpha*k)/hsq)

            end do
        end do

        T = Tn

    END SUBROUTINE blacknodes



    !-------------------------------------------------------------------------------------------!
    !-------------------------------------------------------------------------------------------!
    !-------------------------------------------------------------------------------------------!
    !-------------------------------------------------------------------------------------------!
    !-------------------------------------------------------------------------------------------!
    ! jacobi preconditioner for conjugate gradient
    SUBROUTINE jacobiprecon(ae,aw,an,as,ap,b,Minv,T,il,jl,ih,jh)

        IMPLICIT NONE
        Real(kind = 8), INTENT(IN) :: ae(nx,ny), aw(nx,ny), an(nx,ny), as(nx,ny), ap(nx,ny), b(nx,ny), Minv(nx,ny)
        Real(kind = 8), INTENT(INOUT) :: T(nx,ny)
        Integer(kind = 8), INTENT(IN) :: il,jl,ih,jh

        Real(kind = 8) :: res(nx,ny)
        Integer(kind = 8) :: i,j,ii, niter_precon

        ! Setting number of preconditioning iterations
        niter_precon = 5
        
        do j = 1,ny
        ! do j = jl,jh
            do i = 1,nx
            ! do i = il,ih
                T(i,j) = b(i,j)*Minv(i,j)
            end do
        end do	

        ! do ii = 1,5
        do ii = 1,niter_precon
        
            ! Get Residual of Current system Ax = b
            CALL respar(aw,ae,an,as,ap,b,T,il,ih,jl,jh,res)
            
            ! Update Solution
            do j = 1,ny
            ! do j = jl,jh
                do i = 1,nx 
                ! do i = il,ih
                    T(i,j) = T(i,j) + res(i,j)*Minv(i,j)
                end do
            end do
                
	    end do

    END SUBROUTINE jacobiprecon


end module jacobi