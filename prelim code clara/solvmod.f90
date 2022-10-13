! Solving subroutines: residual obtaining and red-black function

MODULE solvingmodule

	USE variablemodule
	
	! Ensuring all names in the program must be explicitly declared
	IMPLICIT NONE
	
	! --------------------------------------------
	CONTAINS
	! --------------------------------------------
	
	! In the below subroutines:
	! Looping column then row so data of consecutive array locations are contiguous in memory
		
	! --------------------------------------------
	
	SUBROUTINE get_residual_finalval()
		IMPLICIT NONE
		
		! Constant that doesn't change in the loop
		c1 = -2/(dx**2) - 2/(dy**2)
			
		DO j = node_low_x+1, ind_high_x-1
			DO i = ind_low_y+1,ind_high_y-1
				c2 = (alpha*U(i,j))/(2*dx)  -  1/(dx**2)
				c5 = (alpha*V(i,j))/(2*dy)  -  1/(dy**2)
				c4 = (-1)*(alpha*U(i,j))/(2*dx)  -  1/(dx**2)
				c3 = (-1)*(alpha*V(i,j))/(2*dy)  -  1/(dy**2)
				! c3 and c5 are flipped because lower j-indices correspond to high y-values
			
				Res(i,j) = (-1)*( (-1)*c1*Tk(i,j) + c2*Tk(i,j+1) + c3*Tk(i+1,j) + c4*Tk(i,j-1) + c5*Tk(i-1,j) )
			END DO
		END DO
		
		! Take absolute value and sum over all residuals for final value for that processor
		! Divide by total number of grid points
		res_proc_val = (1/(DBLE(nx)*DBLE(ny))) * SUM(ABS(Res))
		
	END SUBROUTINE get_residual_finalval
	
	! --------------------------------------------
	
	SUBROUTINE get_Tk1_Red()
		! Red = chequerboard from 2,2
		
		IMPLICIT NONE
		
		! Constant that doesn't change in the loop
		c1 = -2/(dx**2) - 2/(dy**2)
		
		! Loop through all columns
		DO j = ind_low_x+1, ind_high_x-1
		
			! Loop through rows but only every second one
			DO i = ind_low_y + 1 + MOD(j,2), ind_high_y-1, 2
			
				c2 = (alpha*U(i,j))/(2*dx)  -  1/(dx**2)
				c5 = (alpha*V(i,j))/(2*dy)  -  1/(dy**2)
				c4 = (-1)*(alpha*U(i,j))/(2*dx)  -  1/(dx**2)
				c3 = (-1)*(alpha*V(i,j))/(2*dy)  -  1/(dy**2)
				! c3 and c5 are flipped because lower j-indices correspond to high y-values
				
				Tk_1(i,j) = (c2/c1)*Tk(i,j+1) + (c3/c1)*Tk(i+1,j) + (c4/c1)*Tk(i,j-1) + (c5/c1)*Tk(i-1,j)
				
			END DO

		END DO
	
	END SUBROUTINE get_Tk1_Red
	
	! - - - - - - - - - - - - - - - - - - - - -
	
	SUBROUTINE get_Tk1_Black()
		! Black = chequerboard from 3,2
		
		IMPLICIT NONE
		
		! Constant that doesn't change in the loop
		c1 = -2/(dx**2) - 2/(dy**2)
		
		DO j = ind_low_x+1, ind_high_x-1
		
			! then loop through by row (only middle section!!)
			DO i = ind_low_y +2 -MOD(j,2), ind_high_y-1, 2
			
				c2 = (alpha*U(i,j))/(2*dx)  -  1/(dx**2)
				c5 = (alpha*V(i,j))/(2*dy)  -  1/(dy**2)
				c4 = (-1)*(alpha*U(i,j))/(2*dx)  -  1/(dx**2)
				c3 = (-1)*(alpha*V(i,j))/(2*dy)  -  1/(dy**2)
				! c3 and c5 are flipped because lower j-indices correspond to high y-values
				
				! Use Tk_1 own self because the new additions of the reds have been added in now.
				Tk_1(i,j) = (c2/c1)*Tk_1(i,j+1) + (c3/c1)*Tk_1(i+1,j) + (c4/c1)*Tk_1(i,j-1) + (c5/c1)*Tk_1(i-1,j)
				
			END DO

		END DO
	
	END SUBROUTINE get_Tk1_Black
	
	! --------------------------------------------
	
	! Jacobi subroutine, if needed
	SUBROUTINE get_Tk1_Jacobi()
		IMPLICIT NONE
		
		! Constant that doesn't change in the loop
		c1 = -2/(dx**2) - 2/(dy**2)
		
		DO j = ind_low_x+1, ind_high_x-1
		
			! then loop through by row (only middle section!!)
			DO i = ind_low_y+1,ind_high_y-1
			
				c2 = (alpha*U(i,j))/(2*dx)  -  1/(dx**2)
				c5 = (alpha*V(i,j))/(2*dy)  -  1/(dy**2)
				c4 = (-1)*(alpha*U(i,j))/(2*dx)  -  1/(dx**2)
				c3 = (-1)*(alpha*V(i,j))/(2*dy)  -  1/(dy**2)
				! c3 and c5 are flipped because lower j-indices correspond to high y-values
				
				Tk_1(i,j) = (c2/c1)*Tk(i,j+1) + (c3/c1)*Tk(i+1,j) + (c4/c1)*Tk(i,j-1) + (c5/c1)*Tk(i-1,j)
				
			END DO
				
		END DO
	
	END SUBROUTINE get_Tk1_Jacobi
	

END MODULE solvingmodule