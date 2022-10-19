module cjgradient
    ! NOTE - need equation being solved to be symmetric

    ! Calling modules
    USE variablemodule
    USE jacobi
    USE residuals

    IMPLICIT NONE


    ! Variables obtained from variable module
    ! pi
    ! dx,dy
    ! grid size
    ! constants
    ! rcurrent, rmax


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
    
        subroutine CGSolve(an,as,ae,aw,ap,b,x,y,Tin,Tout)

        ! Initialising variables
        Real(kind = 8), INTENT(IN) :: Tin(nx,ny), an(nx,ny), as(nx,ny), ae(nx,ny), aw(nx,ny), ap(nx,ny), b(nx,ny), x(nx), y(ny)
        Real(kind = 8), INTENT(OUT) :: Tout(nx,ny)

        Real(kind = 8), dimension(nx,ny) :: Minv, d, res, delta, delta_o
        Real(kind = 8) :: rcurrent
        Integer :: i,j

        ! Calling solution initialisations
        ! Initialising - ae,aw,as,an,ap,T,x,y
        Call SolInit()

        ! initiliasing residuals
        call residcalc(aw,ae,an,as,ap,b,x,res)
          
        ! Calculate Domain averaged residual for stopping critterion
        ! rcurrent = SUM(SUM(ABS(r(il:ih,jl:jh)),1),1) / ((kx-2)*(ky-2))
        rcurrent = SUM(SUM(ABS(res(1:nx,1:ny)),1),1) / ((nx-2)*(ny-2))

        ! Calculating d matrix
        do j = 1,ny
            do i = 1,nx
                Minv(i,j) = 1/ap(i,j)
                d(i,j) = Minv(i,j)*res(i,j)
            end do
        end do

        ! Calculating delta and delta0
        do j = 1,ny
            do i = 1,nx
                
                





        end subroutine CGsolve



end module cjgradient