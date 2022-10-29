! Just the tecplot code we are given (modified)

MODULE outputmodule

		! nx 	number of points in x direction
		! ny 	number of points in y direction
		! x		vector of coordinates in x direction
		! y		vector of coordinates in y direction
		! T 	2D matrix with temperatures in x and y directions
	
	! Ensuring all names in the program must be explicitly declared
	IMPLICIT NONE
	
	! --------------------------------------------
	CONTAINS
	! --------------------------------------------


	SUBROUTINE tecplot_2D ( iunit, nx, ny, x, y, T )

		IMPLICIT NONE

		Integer :: iunit,nx,ny,i,j,ierr
		Real ( kind = 8 ) T(ny,nx)
		Real ( kind = 8 ) x(nx)
		Real ( kind = 8 ) y(ny)

		Character(80), parameter ::  file_name = 'TecPlot2D.tec'

		open ( unit = iunit, file = file_name, form = 'formatted', access = 'sequential', status = 'replace', iostat = ierr )

		if ( ierr /= 0 ) then
			write ( *, '(a)' ) '  Error opening file : tecplot_2D '
			stop
		end if
		   
		write ( iunit, '(a)' ) 'Title="' // trim ( 'Temperature Data' ) // '"'
		write ( iunit, '(a)' ) 'Variables=' // trim ( '"X","Y","T"' )

		write ( iunit, '(a)' ) ' '
		write ( iunit, '(a,i6,a,i6,a)' ) 'Zone I=', ny, ', J=', nx, ', F=POINT'
		 
		do i = 1, nx ! loop through columns left to right
			do j = 1, ny ! loop through rows top to bottom
				write ( iunit, '(2f10.3,g15.6)' ) x(i), y(j), T(j,i)
			end do
		end do
		  
		close ( unit = iunit )

	END SUBROUTINE tecplot_2D
	
END MODULE outputmodule