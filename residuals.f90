module residuals

    ! NOTE THAT THE INPUT OF THIS SHOULD BE THE NODES IN THE PROCESSOR DOMAIN EXCLUDING GHOST NODES
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

    ! parallel residual calculation
    subroutine respar(aw,ae,an,as,ap,b,T,il,ih,jl,jh,res)

        ! Initialising variables
        Integer(kind = 8), INTENT(IN) :: il,ih,jl,jh
        real(kind=8), dimension(il:ih,jl:jh), INTENT(IN) :: aw,ae,an,as,ap,b,T
        real(kind=8), dimension(il:ih,jl:jh), INTENT(OUT) :: res
        
        do j = jl+1,jh-1
            do i = il+1,ih-1           
                res(i,j) = b(i,j) - ( as(i,j)*T(i-1,j) + an(i,j)*T(i+1,j) + ae(i,j)*T(i,j+1) &
                                                    + aw(i,j)*T(i,j-1) + ap(i,j)*T(i,j) )		           
            end do
        end do	

    end subroutine respar


end module residuals