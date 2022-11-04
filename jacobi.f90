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
    SUBROUTINE jac(an,as,ae,aw,ap,b,T,Told,il,ih,jl,jh)


        ! Defining variables
        Real(kind=8), INTENT(IN), dimension(il:ih,jl:jh) :: an,as,ae,aw,ap,b,Told
        Integer(kind = 8), INTENT(IN) :: il,ih,jl,jh
        Real(kind=8), INTENT(INOUT) :: T(il:ih,jl:jh)

        Integer :: i,j

        do j = jl+1,jh-1
            do i = il+1,ih-1            
                Tn(i,j) = (1- (4*dt*alpha)/(dx*dx))*Told(i,j) + ((dt*alpha)/(dx*dx))*(T(i+1,j) + T(i-1,j) + T(i,j+1) + T(i,j-1))
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
    SUBROUTINE rednodes(an,as,ae,aw,ap,b,T,Told,il,ih,jl,jh)

        ! Defining variables
        Real(kind=8), INTENT(IN), dimension(il:ih,jl:jh) :: an,as,ae,aw,ap,b,Told
        Integer(kind = 8), INTENT(IN) :: il,ih,jl,jh
        Real(kind=8), INTENT(INOUT) :: T(il:ih,jl:jh)

        Integer :: i,j

        do j = jl+1,jh-1
            do i = il+1 + mod(j,2),ih-1,2
                Tn(i,j) = ap(i,j)*Told(i,j) + an(i,j)*(T(i+1,j) + T(i-1,j) + T(i,j+1) + T(i,j-1))
            end do
        end do

        T = Tn

    END SUBROUTINE rednodes

    ! black nodes
    SUBROUTINE blacknodes(an,as,ae,aw,ap,b,T,Told,il,ih,jl,jh)

        ! Defining variables
        Real(kind=8), INTENT(IN), dimension(il:ih,jl:jh) :: an,as,ae,aw,ap,b,Told
        Integer(kind = 8), INTENT(IN) :: il,ih,jl,jh
        Real(kind=8), INTENT(INOUT) :: T(il:ih,jl:jh)

        Integer :: i,j

        do j = jl+1,jh-1
            do i = il+2 - mod(j,2),ih-1,2
                Tn(i,j) = (1- (4*dt*alpha)/(dx*dx))*Told(i,j) + ((dt*alpha)/(dx*dx))*(T(i+1,j) + T(i-1,j) + T(i,j+1) + T(i,j-1))
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
    subroutine jacpre(ae,aw,an,as,ap,b,Minv,T,il,ih,jl,jh)

        IMPLICIT NONE
        Real(kind = 8), INTENT(IN), dimension(il:ih,jl:jh) :: ae, aw, an, as, ap, b, Minv
        Real(kind = 8), INTENT(INOUT), dimension(il:ih,jl:jh) :: T
        Integer(kind = 8), INTENT(IN) :: il,ih,jl,jh

        Integer(kind = 8) :: i,j,ii, niter_precon

        ! Setting number of preconditioning iterations
        niter_precon = 5
        
        do j = jl,jh
            do i = il,ih
                T(i,j) = b(i,j)*Minv(i,j)
            end do
        end do	

        do ii = 1,niter_precon
        
            ! Get Residual of Current system Ax = b
            CALL respar(aw,ae,an,as,ap,b,T,il,ih,jl,jh,resmat)
            
            ! Update Solution
            do j = jl,jh
                do i = il,ih
                    T(i,j) = T(i,j) + resmat(i,j)*Minv(i,j)
                end do
            end do
                
	    end do

    end subroutine


    ! b - matrix
subroutine printmatrix(b,n,m)
	integer::n,m
	real (kind=8)::b(n,m) !n = # rows, m = # columns
	do i=1,n; print '(20f16.8)',b(i,1:m); enddo
endsubroutine


end module jacobi