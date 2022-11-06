module solinit
! Initialises the solution

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

        ! Allocating sizes to variables
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
    subroutine solutioninit(an,as,ae,aw,ap,b,T,il,ih,jl,jh)

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
        ! Real(kind = 8), allocatable :: Tn(:,:)
        Real(kind = 8) :: hsq, k

        ! F = (dt*alpha)/dx**2
        ! Tn(il:ih,jl:jh) = 1
        hsq = dx**2
        k = dt


        ! Computing A and B matrices - need matrices for conjugate gradient method
        do j = jl,jh
            do i = il,ih
                ! IMPLICIT
                ap(i,j) = (1- (4*dt*alpha)/(dx*dx))
                an(i,j) = ((dt*alpha)/(dx*dx))
                as(i,j) = ((dt*alpha)/(dx*dx))
                ae(i,j) = ((dt*alpha)/(dx*dx))
                aw(i,j) = ((dt*alpha)/(dx*dx))
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

        PI=4.D0*DATAN(1.D0)
        ! Setting solver boundary conditions
        do j = jl,jh
            T(1,j) = sin(((pi*(j-1)*dx)/Lx))
        end do

        ! Edge boundary conditions
        if (ih.eq.ny) then
            T(ny,:) = 0
        elseif (jl.eq.1) then
            T(:,1) = 0
        elseif (jh.eq.nx) then
            T(:,nx) = 0
        end if


    end subroutine solutioninit



end module solinit