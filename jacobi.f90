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
    ! USE nodemodule
    ! USE partitionmodule


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
    SUBROUTINE jac(an,as,ae,aw,ap,b,T,il,ih,jl,jh,time)


        ! Defining variables
        Real(kind=8), INTENT(IN), dimension(il:ih,jl:jh) :: an,as,ae,aw,ap,b
        Real(kind=8), INTENT(IN) :: time
        Integer(kind = 8), INTENT(IN) :: il,ih,jl,jh
        Real(kind=8), INTENT(INOUT) :: T(il:ih,jl:jh)

        Real(kind = 8) :: Told(il:ih,jl:jh), Tn(il:ih,jl:jh)
        Integer :: i,j


        if (time == 0) then
            do j = jl,jh
                do i = il,ih
                    Tn(i,j) = (T(i+1,j)*ae(i,j) + T(i-1,j)*aw(i,j) + T(i,j+1)*an(i,j) &
                        + T(i,j-1)*as(i,j) + T(i,j)*ap(i,j))/2
                end do
            end do
        else
            do j = jl,jh
                do i = il,ih
                    Tn(i,j) = T(i+1,j)*ae(i,j) + T(i-1,j)*aw(i,j) + T(i,j+1)*an(i,j) &
                        + T(i,j-1)*as(i,j) + T(i,j)*ap(i,j) - Told(i,j)
                end do
            end do
        end if

    END SUBROUTINE jac

    !-------------------------------------------------------------------------------------------!
    !-------------------------------------------------------------------------------------------!
    !-------------------------------------------------------------------------------------------!
    !-------------------------------------------------------------------------------------------!
    !-------------------------------------------------------------------------------------------!
    ! parallel redblack solver
    ! red nodes
    SUBROUTINE rednodes(an,as,ae,aw,ap,b,T,il,ih,jl,jh,time)

        ! Defining variables
        Real(kind=8), INTENT(IN), dimension(il:ih,jl:jh) :: an,as,ae,aw,ap,b
        Real(kind=8), INTENT(IN) :: time
        Integer(kind = 8), INTENT(IN) :: il,ih,jl,jh
        Real(kind=8), INTENT(INOUT) :: T(il:ih,jl:jh)

        Real(kind=8) :: Told(il:ih,jl:jh), Tn(il:ih,jl:jh)
        Integer :: i,j


        if (time == 0) then
            ! RED NODES CALCULATION
            do j = 2,(ny-1),2
            ! do j = jl,jh,2
                do i = 2,(nx-1),2
                ! do i = il,ih,2
                    Tn(i,j) = (T(i+1,j)*ae(i,j) + T(i-1,j)*aw(i,j) + T(i,j+1)*an(i,j) &
                        + T(i,j-1)*as(i,j) + T(i,j)*ap(i,j))/2
                end do
            end do
        else
            do j = 2,(ny-1),2
            ! do j = jl,jh,2
                do i = 2,(nx-1),2
                ! do i = il,ih,2
                    Tn(i,j) = T(i+1,j)*ae(i,j) + T(i-1,j)*aw(i,j) + T(i,j+1)*an(i,j) &
                        + T(i,j-1)*as(i,j) + T(i,j)*ap(i,j) - Told(i,j)
                end do
            end do
        end if

    END SUBROUTINE rednodes

    ! black nodes
    SUBROUTINE blacknodes(an,as,ae,aw,ap,b,T,il,ih,jl,jh,time)

        ! Defining variables
        Real(kind=8), INTENT(IN), dimension(il:ih,jl:jh) :: an,as,ae,aw,ap,B
        Real(kind=8), INTENT(IN) :: time
        Integer(kind = 8), INTENT(IN) :: il,ih,jl,jh
        Real(kind=8), INTENT(INOUT) :: T(il:ih,jl:jh)

        Real(kind=8) :: Told(il:ih,jl:jh), Tn(il:ih,jl:jh)
        Integer :: i,j


        if (time == 0) then
            ! BLACK NODES CALCULATION
            do j = 3,(ny-1),2
            ! do j = jl+1,jh,2
                do i = 3,(nx-1),2
                ! do i = il+1,ih,2
                    Tn(i,j) = (T(i+1,j)*ae(i,j) + T(i-1,j)*aw(i,j) + T(i,j+1)*an(i,j) &
                        + T(i,j-1)*as(i,j) + T(i,j)*ap(i,j))/2
                end do
            end do
        else
            do j = 3,(ny-1),2
            ! do j = jl+1,jh,2
                do i = 3,(nx-1),2
                ! do i = il+1,ih,2
                    Tn(i,j) = T(i+1,j)*ae(i,j) + T(i-1,j)*aw(i,j) + T(i,j+1)*an(i,j) &
                        + T(i,j-1)*as(i,j) + T(i,j)*ap(i,j) - Told(i,j)
                end do
            end do
        end if

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