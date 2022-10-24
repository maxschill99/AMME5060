module cgradient
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
    
        subroutine CGSolve(an,as,ae,aw,ap,b,T)

        ! Initialising variables
        Real(kind = 8), INTENT(IN) :: an(nx,ny), as(nx,ny), ae(nx,ny), aw(nx,ny), ap(nx,ny), b(nx,ny)
        Real(kind = 8), INTENT(INOUT) :: T(nx,ny)

        Real(kind = 8), dimension(nx,ny) :: Minv, d, res, q, s, Tin, Tout
        Real(kind = 8) :: delta, delta_o, dp, alpha, beta
        Real(kind = 8) :: rcurrent
        Integer :: i,j,niter_precon, il,jl,ih,jh
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

        ! Calling solution initialisations
        ! Initialising - ae,aw,as,an,ap,T,x,y
        ! call solninit(an,as,ae,aw,ap,b,T)

        ! initiliasing residuals
        call residcalc(aw,ae,an,as,ap,b,T,res)
          
        ! Calculate Domain averaged residual for stopping critterion
        ! rcurrent = SUM(SUM(ABS(r(il:ih,jl:jh)),1),1) / ((kx-2)*(ky-2))
        rcurrent = SUM(SUM(ABS(res(1:nx,1:ny)),1),1) / ((nx-2)*(ny-2))

        ! Calculating d matrix (search vector matrix)
        do j = 1,ny
        ! do j = jl,jh
            do i = 1,nx
            ! do i = il,ih
                Minv(i,j) = 1/ap(i,j)
            end do
        end do

        if (precon) then
            do j = 1,ny
            ! do j = jl,jh
                do i = 1,nx
                ! do i = il,ih
                    d(i,j) = Minv(i,j)*res(i,j)
                end do
            end do
        end if

        ! Pre-conditioning matrix
        if (precon) then
            call jacobiprecon(ae,aw,an,as,ap,b,Minv,d)
        else
            d = res
        end if

        ! Calculating delta and delta0
        do j = 1,ny
        ! do j = jl,jh
            do i = 1,nx
            ! do i = il,ih
                dp = dp + res(i,j)*d(i,j)
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
            do j = 2,(ny-1)
            ! do j = jl,jh
                do i = 2,(nx-1)
                ! do i = il,ih
                    q(i,j) = ae(i,j)*d(i-1,j) + aw(i,j)*d(i+1,j) + an(i,j)*d(i,j+1) + as(i,j)*d(i,j-1) &
                               + ap(i,j)*d(i,j)
                end do
            end do


            dp = 0
            ! compute alpha - note computing dot product in the numerator
            do j = 1,ny
            ! do j = jl,jh
                do i = 1,nx
                ! do i = il,ih
                    dp = dp + d(i,j)*q(i,j)
                end do
            end do
            alpha = delta/dp

            ! updating T
            do j = 1,ny
            ! do j = jl,jh
                do i = 1,nx
                ! do i = il,ih
                    T(i,j) = T(i,j) + alpha*d(i,j)
                end do
            end do


            ! Updating residual
            if ((MOD(iter,50).eq.0)) then
                call residcalc(aw,ae,an,as,ap,b,T,res)
            else
                do j = 1,ny
                ! do j = jl,jh
                    do i = 1,nx
                    ! do i = il,ih
                        res(i,j) = res(i,j) - alpha*q(i,j)
                    end do
                end do
            end if


            ! Pre-conditioning s matrix
            if (precon) then
                call jacobiprecon(ae,aw,an,as,ap,b,Minv,s)
            else
                s = res
            end if


            dp = 0
            ! Calculating delta and delta0
            do j = 1,ny
            ! do j = jl,jh
                do i = 1,nx
                ! do i = il,ih
                    dp = dp + res(j,i)*s(i,j)
                end do
            end do
            
            delta_o = delta
            delta = dp


            ! Computing beta
            beta = delta/delta_o
            
            ! Updating search vector
            do j = 1,ny
            ! do j = jl,jh
                do i = 1,nx
                ! do i = il,ih
                    d(i,j) = s(i,j) + beta*d(i,j)
                end do
            end do

            ! Update residual check	
		    rcurrent = SUM(SUM(ABS(res(1:nx,1:ny)),1),1) / ((nx-2)*(ny-2))	


            ! Updating timer and iteration counter
            time = time + dt
            iter = iter + 1

            WRITE(*,'(a20,i20.1)') '# Current Iteration = ', iter
            WRITE(*,'(a20,E20.6)') '# Current Residual = ', rcurrent
            WRITE(*,'(a20,f20.10)') '# Current Max Temp = ', MAXVAL(T)

        end do 
                
                

        end subroutine CGsolve



end module cgradient
