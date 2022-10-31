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

        do j = jl+1,jh-1
            do i = il+1,ih-1
                ! UNSTEADY
                Tn(i,j) = T(i+1,j)*an(i,j) + T(i-1,j)*as(i,j) + T(i,j+1)*ae(i,j) &
                    + T(i,j-1)*aw(i,j) + T(i,j)*ap(i,j)
                ! Tn(i,j) = T(i,j)*ap(i,j) + alpha*dt*((T(i,j+1)+T(i,j-1))/(dx**2)) + alpha*dt*((T(i+1,j)+T(i-1,j))/(dy**2))  

                ! STEADY
                ! Tn(i,j) = (alpha*((T(i,j+1)+T(i,j-1))/(dx**2)) & 
                !             + alpha*((T(i+1,j)+T(i-1,j))/(dy**2)))/((2*alpha)/dx**2 + (2*alpha)/dy**2)       
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

        do j = jl+1,jh-1
            do i = il+1 + mod(j,2),ih-1,2
                ! UNSTEADY
                Tn(i,j) = T(i+1,j)*an(i,j) + T(i-1,j)*as(i,j) + T(i,j+1)*ae(i,j) &
                    + T(i,j-1)*aw(i,j) + T(i,j)*ap(i,j)

                ! ! STEADY
                ! Tn(i,j) = (alpha*((T(i,j+1)+T(i,j-1))/(dx**2)) & 
                !             + alpha*((T(i+1,j)+T(i-1,j))/(dy**2)))/((2*alpha)/dx**2 + (2*alpha)/dy**2) 
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

        do j = jl+1,jh-1
            do i = il+2 - mod(j,2),ih-1,2
                Tn(i,j) = T(i+1,j)*an(i,j) + T(i-1,j)*as(i,j) + T(i,j+1)*ae(i,j) &
                    + T(i,j-1)*aw(i,j) + T(i,j)*ap(i,j)
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
    SUBROUTINE jacobiprecon(ae,aw,an,as,ap,b,Minv,T)

        IMPLICIT NONE
        Real(kind = 8), INTENT(IN) :: ae(nx,ny), aw(nx,ny), an(nx,ny), as(nx,ny), ap(nx,ny), b(nx,ny), Minv(nx,ny)
        Real(kind = 8), INTENT(INOUT) :: T(nx,ny)

        Real(kind = 8) :: res(nx,ny)
        Integer(kind = 8) :: i,j,ii, niter_precon, il,jl,ih,jh

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
            CALL residcalc(aw,ae,an,as,ap,b,T,res)
            
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