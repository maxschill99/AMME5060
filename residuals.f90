module residuals

    ! Calling modules
    USE variablemodule
    USE jacobi

    IMPLICIT NONE

    ! List of subroutines within module
    !-------------------------------------------------------------------------------------------!
    CONTAINS
    !-------------------------------------------------------------------------------------------!

    subroutine residcalc(aw,ae,an,as,ap,b,x,r)
        IMPLICIT NONE
        REAL(KIND = 8), DIMENSION(kx,ky), INTENT(IN) :: aw,ae,an,as,ap,b,x
        REAL(KIND = 8), DIMENSION(kx,ky), INTENT(OUT) :: r
        
        INTEGER :: i,j
        
        r(:,:) = 0.0
        
        DO j = jl,jh
            DO i = il,ih
            
                r(i,j) = b(i,j) - ( aw(i,j)*x(i-1,j) + ae(i,j)*x(i+1,j) + an(i,j)*x(i,j+1) &
                                                    + as(i,j)*x(i,j-1) + ap(i,j)*x(i,j) )		
            
            ENDDO
        ENDDO		

    end subroutine residcalc


end module residuals