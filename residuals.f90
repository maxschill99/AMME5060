module residuals

    ! Calling modules
    USE variablemodule
    ! USE jacobi

    IMPLICIT NONE

    ! List of subroutines within module
    !-------------------------------------------------------------------------------------------!
    CONTAINS
    !-------------------------------------------------------------------------------------------!
    !-------------------------------------------------------------------------------------------!
    !-------------------------------------------------------------------------------------------!
    !-------------------------------------------------------------------------------------------!
    !-------------------------------------------------------------------------------------------!
    !-------------------------------------------------------------------------------------------!

    subroutine residcalc(aw,ae,an,as,ap,b,x,res)
        IMPLICIT NONE

        ! Initialising variables
        real(kind=8), dimension(nx,ny), INTENT(IN) :: aw,ae,an,as,ap,b,x
        real(kind=8), dimension(nx,ny), INTENT(OUT) :: res
        
        INTEGER :: i,j
        
        ! Initialising residual
        res(:,:) = 0.0
        
        do j = 2,(ny-1)
        ! do j = jl,jh
            do i = 2,(nx-1)
            ! do i = il,ih           
                res(i,j) = b(i,j) - ( aw(i,j)*x(i-1,j) + ae(i,j)*x(i+1,j) + an(i,j)*x(i,j+1) &
                                                    + as(i,j)*x(i,j-1) + ap(i,j)*x(i,j) )		           
            end do
        end do		

    end subroutine residcalc


    subroutine parres(aw,ae,an,as,ap,b,T,il,ih,jl,jh,res)

        ! Initialising variables
        Integer, INTENT(IN) :: il,ih,jl,jh
        real(kind=8), dimension(il:ih,jl:jh), INTENT(IN) :: aw,ae,an,as,ap,b,T
        real(kind=8), dimension(il:ih,jl:jh), INTENT(OUT) :: res
        
        ! Initialising residual
        res(:,:) = 0.0
        
        do j = jl,jh
            do i = il,ih           
                res(i,j) = b(i,j) - ( aw(i,j)*T(i-1,j) + ae(i,j)*T(i+1,j) + an(i,j)*T(i,j+1) &
                                                    + as(i,j)*T(i,j-1) + ap(i,j)*T(i,j) )		           
            end do
        end do	

    end subroutine parres


end module residuals