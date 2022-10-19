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


end module residuals