module jacobi

    ! Calling modules
    USE variablemodule

    IMPLICIT NONE

    ! List of subroutines within module
    !-------------------------------------------------------------------------------------------!
    CONTAINS
    !-------------------------------------------------------------------------------------------!

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

        REAL, INTENT(IN) :: Tin(nx,ny)
        REAL, INTENT(OUT) :: Tout(nx,ny)

        Real(kind = 8) :: an, as, ae, aw, ap
        Real(kind = 8) :: time
        Real(kind = 8) :: T(nx,ny), Tn(nx,ny), Told(nx,ny)
        Integer :: i,j


        !------------------------------------------------------------------------------------!
        !------------------------------------------------------------------------------------!
        

        ! simplifying variable notation
        an = (2*dt*alpha)/dy**2
        as = (2*dt*alpha)/dy**2
        ae = (2*dt*alpha)/dx**2
        aw = (2*dt*alpha)/dx**2
        ap = (-2/dx**2 - 2/dy**2)*2*dt*alpha

        ! Inititialising old temperature array
        Told(:,:) = 0
        T = Tin

        ! Initialising time counter
        time = 0

        ! Iteratively solving the jacobi equation
        ! Solving unsteady 2D Heat diffusion - dT/dt = alpha*(d^2T/dx^2 + d^2T/dy^2)
        ! do while ((t<ttot).and.(res>rmax))
        do while (time<t_final)
            do i = 2,(nx-1)
                do j = 2,(ny-1)
                    Tn(i,j) = T(i+1,j)*ae + T(i-1,j)*aw + T(i,j+1)*an + T(i,j-1)*as + T(i,j)*ap - Told(i,j)
                end do
            end do

            ! updating old and new temperature values
            Told = T
            T = Tn

            ! updating counter
            time = time + dt

            ! computing residuals


        end do

        ! Outputting the temp array
        Tout = T

    END SUBROUTINE


end module jacobi