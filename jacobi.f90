module jacobi

    ! Calling modules
    USE variablemodule

    IMPLICIT NONE

    ! List of subroutines within module
    !-------------------------------------------------------------------------------------------!
    CONTAINS
    !-------------------------------------------------------------------------------------------!

    SUBROUTINE solninit()

        ! Real(kind = 8) INTENT(OUT) :: an,as,ae,aw,ap,b

        Real(kind = 8) :: an, as, ae, aw, ap, b

        ! Computing A and B matrices - need matrices for conjugate gradient method
        do j = 1,ny
            do i = 1,nx
                ! simplifying variable notation
                an(i,j) = (2*dt*alpha)/dy**2
                as(i,j) = (2*dt*alpha)/dy**2
                ae(i,j) = (2*dt*alpha)/dx**2
                aw(i,j) = (2*dt*alpha)/dx**2
                ap(i,j) = (-2/dx**2 - 2/dy**2)*2*dt*alpha
                b(i,j) = 0
            end do
        end do

    END SUBROUTINE solninit

    SUBROUTINE solver(Tin, Tout)

        IMPLICIT NONE
        ! Getting variables from varmod module
            ! nx - number of x points
            ! ny - number of y points
            ! dx - discretisation in x
            ! dy - discretisation in y
            ! dt - temporal discretisation
            ! dt - temporal discretisation
            ! alpha - thermal diffusivity
            ! t_final - total time

        REAL(kind=4), INTENT(IN) :: Tin(nx,ny)
        REAL(kind=4), INTENT(OUT) :: Tout(nx,ny)

        Real(kind = 8) :: time
        Real(kind = 8) :: T(nx,ny), Tn(nx,ny), Told(nx,ny)
        Integer :: i,j, iter


        !------------------------------------------------------------------------------------!
        !------------------------------------------------------------------------------------!

        ! Getting solution initialisation variables
        CALL solninit()

        ! Inititialising old temperature array
        Told(:,:) = 0
        T = Tin

        ! Initialising time counter
        time = 0
        iter = 0

        ! Iteratively solving the jacobi equation
        ! Solving unsteady 2D Heat diffusion - dT/dt = alpha*(d^2T/dx^2 + d^2T/dy^2)
        ! do while ((t<ttot).and.(res>rmax))
        do while (time<t_final)
            do j = jl,jh
                do i = il,ih
                    Tn(i,j) = T(i+1,j)*ae(i,j) + T(i-1,j)*aw(i,j) + T(i,j+1)*an(i,j) + T(i,j-1)*as(i,j) + T(i,j)*ap(i,j) - Told(i,j)
                end do
            end do

            if (mod(iter,100).eq.0) then
                write(*,*) iter
            end if

            ! updating old and new temperature values
            Told = T
            T = Tn

            ! updating counter
            time = time + dt
            iter = iter + 1

            ! computing residuals


        end do

        ! Outputting the temp array
        Tout = T

    END SUBROUTINE solver


end module jacobi