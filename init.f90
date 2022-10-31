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
        Integer(kind = 8), INTENT(IN) :: il,ih,jl,jh

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

        Integer(kind = 8), INTENT(IN) :: il,ih,jl,jh
        Real(kind = 8), INTENT(OUT) :: an(il:ih,jl:jh), as(il:ih,jl:jh), ae(il:ih,jl:jh), aw(il:ih,jl:jh), &
         ap(il:ih,jl:jh), b(il:ih,jl:jh), T(il:ih,jl:jh)

        Real(kind = 8), allocatable :: x(:), y(:)


        ! Computing A and B matrices - need matrices for conjugate gradient method
        ! do j = 1,ny
        do j = jl,jh
            ! do i = 1,nx
            do i = il,ih
                ! simplifying variable notation
                an(i,j) = (dt*alpha)/dy**2
                as(i,j) = (dt*alpha)/dy**2
                ae(i,j) = (dt*alpha)/dx**2
                aw(i,j) = (dt*alpha)/dx**2
                ap(i,j) = 1 - (2*dt*alpha)/dx**2 - (2*dt*alpha)/dy**2
                b(i,j) = 0
            end do
        end do

        ! Allocating and initialising position arrays
        allocate(x(jl:jh))
        allocate(y(il:ih))
        ! x is in the j direction, y is in the i direction
        do i = il,ih
            y(i) = (i-1)*dy
        end do
        do j = jl,jh
            x(j) = (j-1)*dx
        end do

        ! Initialising temperature matrix
        T(:,:) = 0.0

        ! allocate(x(jl:jh))
        ! do j = jl,jh
        !     x(j) = (j-1)*dx
        ! end do

        ! Setting solver boundary conditions
        do j = jl,jh
            T(1,j) = sin((pi*x(j))/Lx)
        end do

        ! Edge boundary conditions
        if (ih.eq.ny) then
            T(ny,:) = 0
        elseif (jl.eq.1) then
            T(:,1) = 0
        elseif (jh.eq.nx) then
            T(:,nx) = 0
        end if


    END SUBROUTINE solutioninit



end module solinit