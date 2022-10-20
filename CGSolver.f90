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

        Real(kind = 8), dimension(nx,ny) :: Minv, d, res, q, s, T
        Real(kind = 8) :: delta, delta_o, dp, alpha, beta
        Real(kind = 8) :: rcurrent
        Integer :: i,j,niter_precon
        Integer :: iter
        Real(kind = 8) :: time
        
        Logical :: precon 


        ! Pre-conditioner counter and check
        niter_precon = 5
	    precon = .TRUE.

        ! Initialising solution constants and vectors
        ! Initialising delta and delta0
        delta = 1
        delta_o = 1
        ! Initialising dot product variable
        dp = 0
        ! Initialising alpha and beta
        alpha = 0
        beta = 0
        ! Initialising q and s
        q(:,:) = 0
        s(:,:) = 0
        ! Initialising temp array
        T = Tin

        ! Calling solution initialisations
        ! Initialising - ae,aw,as,an,ap,T,x,y
        Call SolInit()

        ! initiliasing residuals
        call residcalc(aw,ae,an,as,ap,b,x,res)
          
        ! Calculate Domain averaged residual for stopping critterion
        ! rcurrent = SUM(SUM(ABS(r(il:ih,jl:jh)),1),1) / ((kx-2)*(ky-2))
        rcurrent = SUM(SUM(ABS(res(1:nx,1:ny)),1),1) / ((nx-2)*(ny-2))

        ! Calculating d matrix (search vector matrix)
        do j = 1,ny
            do i = 1,nx
                Minv(i,j) = 1/ap(i,j)
                d(i,j) = Minv(i,j)*res(i,j)
            end do
        end do

        ! Pre-conditioning matrix
        if (precon) then
            jacobisolv(an,as,ae,aw,ap,b,x,y,Tin,Tout)
        else
            d = res
        end if

        ! Calculating delta and delta0
        do j = 1,ny
            do i = 1,nx
                dp = dp + r(j,i)*d(i,j)
            end do
        end do

        ! Assigning values for delta and delta_o
        delta = dp
        delta_o = delta

        ! Initialisint timer and iteration counter
        time = 0
        iter = 0


        !-------------------------------------------------------------------------!
        !-------------------------------------------------------------------------!
        !-------------------------------------------------------------------------!
        ! Begin solution loop
        do while ((time < t_final).and.(rcurrent < res_max))

            ! Compute q matrix
            do j = 1,ny
                do i = 1,nx
                    q(i,j) = ae(i,j)*d(i-1,j) + aw(i,j)*d(i+1,j) + an(i,j)*d(i,j+1) + as(i,j)*d(i,j-1) &
                                ap(i,j)*d(i,j)
                end do
            end do


            dp = 0
            ! compute alpha - note computing dot product in the numerator
            do j = 1,ny
                do i = 1,nx
                    dp = dp + d(j,i)*q(i,j)
                end do
            end do
            alpha = delta/dp

            ! updating T
            do j = 1,ny
                do i = 1,nx
                    T(i,j) = T(i,j) + alpha*d(i,j)
                end do
            end do


            ! Updating residual
            if ((MOD(iter,50).eq.0)) then
                call residcalc(aw,ae,an,as,ap,b,x,res)
            else
                do j = 1,ny
                    do i = 1,nx
                        res(i,j) = res(i,j) - alpha*q(i,j)
                    end do
                end do
            end if


            ! Pre-conditioning s matrix
            if (precon) then
                jacobisolv(an,as,ae,aw,ap,b,x,y,Tin,Tout)
            else
                s = r
            end


            dp = 0
            ! Calculating delta and delta0
            do j = 1,ny
                do i = 1,nx
                    dp = dp + r(j,i)*s(i,j)
                end do
            end 
            
            delta_o = delta
            delta = dp


            ! Computing beta
            beta = delta/delta_o
            
            ! Updating search vector
            do j = 1,ny
                do i = 1,nx
                    d(i,j) = s(i,j) + beta*d(i,j)
                end do
            end do

            ! Update residual check	
		    rcurrent = SUM(SUM(ABS(r(1:nx,1,ny)),1),1) / ((nx-2)*(ny-2))	


            ! Updating timer and iteration counter
            time = time + dt
            iter = iter + 1

            WRITE(*,'(a20,i20.1)') '# Current Iteration = ', iter
            WRITE(*,'(a20,E20.6)') '# Current Residual = ', rcurrent
            WRITE(*,'(a20,f20.10)') '# Current Max Temp = ', MAXVAL(t)

        end do 
                
                





        end subroutine CGsolve



end module cjgradient
