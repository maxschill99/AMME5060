module solinit

USE variablemodule

 CONTAINS

    !-------------------------------------------------------------------------------------------!
    !-------------------------------------------------------------------------------------------!
    !-------------------------------------------------------------------------------------------!
    !-------------------------------------------------------------------------------------------!
    !-------------------------------------------------------------------------------------------!
    !-------------------------------------------------------------------------------------------!
    ! subroutine to allocate variables
    subroutine allocatevars(an,as,ae,aw,ap,b,T,Told,Tn,res,il,ih,jl,jh)

        Real(kind = 8), allocatable, INTENT(INOUT) :: an(:,:), as(:,:), ae(:,:), aw(:,:), ap(:,:), &
                                                      b(:,:), T(:,:), Told(:,:), Tn(:,:), res(:,:)
        Integer, INTENT(IN) :: il,ih,jl,jh

        allocate(an(il:ih,jl:jh))
        allocate(as(il:ih,jl:jh))
        allocate(ae(il:ih,jl:jh))
        allocate(aw(il:ih,jl:jh))
        allocate(ap(il:ih,jl:jh))
        allocate(b(il:ih,jl:jh))
        allocate(T(il:ih,jl:jh))
        allocate(Told(il:ih,jl:jh))
        allocate(Tn(il:ih,jl:jh))
        allocate(res(il:ih,jl:jh))

    end subroutine allocatevars

    !-------------------------------------------------------------------------------------------!
    !-------------------------------------------------------------------------------------------!
    !-------------------------------------------------------------------------------------------!
    !-------------------------------------------------------------------------------------------!
    !-------------------------------------------------------------------------------------------!
    !-------------------------------------------------------------------------------------------!
    ! subroutine to initialise solution
    SUBROUTINE solutioninit(an,as,ae,aw,ap,b,T,il,ih,jl,jh)

        ! ! Declaring local variables
        ! Integer :: il,ih,jl,jh
        ! il = ind_low_x
        ! ih = ind_high_x
        ! jl = ind_low_y
        ! jh = ind_high_y

        Integer, INTENT(IN) :: il,ih,jl,jh
        ! Real(kind = 8), INTENT(OUT) :: an(nx,ny), as(nx,ny), ae(nx,ny), aw(nx,ny), ap(nx,ny), b(nx,ny)
        Real(kind = 8), INTENT(OUT) :: an(il:ih,jl:jh), as(il:ih,jl:jh), ae(il:ih,jl:jh), aw(il:ih,jl:jh), &
         ap(il:ih,jl:jh), b(il:ih,jl:jh), T(il:ih,jl:jh)

        
        ! Real(kind = 8) :: an(nx,ny), as(nx,ny), ae(nx,ny), aw(nx,ny), ap(nx,ny), b(nx,ny)
        ! Real(kind = 8) :: Tin(nx,ny), Tout(nx,ny)


        ! Computing A and B matrices - need matrices for conjugate gradient method
        ! do j = 1,ny
        do j = jl,jh
            ! do i = 1,nx
            do i = il,ih
                ! simplifying variable notation
                an(i,j) = (2*dt*alpha)/dy**2
                as(i,j) = (2*dt*alpha)/dy**2
                ae(i,j) = (2*dt*alpha)/dx**2
                aw(i,j) = (2*dt*alpha)/dx**2
                ap(i,j) = (-2/dx**2 - 2/dy**2)*2*dt*alpha
                b(i,j) = 0
            end do
        end do

        allocate(x(il:ih))
        allocate(y(jl:jh))

        ! x is in the j direction, y is in the i direction
        do i = il,ih
            y(i) = (i-1)*dy
        end do
        do j = jl,jh
            x(j) = (j-1)*dx
        end do

        ! Setting solver boundary conditions
        T(:,:) = 0

        do j = jl,jh
            do i = il,ih
                T(1,j) = sin(pi*x(j))/Lx
            end do
        end do
        T(:,1) = 0
        T(:,nx) = 0
        T(ny,:) = 0


    END SUBROUTINE solutioninit



end module solinit