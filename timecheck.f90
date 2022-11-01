program timecheck
    real(kind = 8) :: dt,t_final, time

    dt		= 1.86e-5		! Time step size [s]
    t_final = 5.				! [s] &&& just a guess for now - CHANGE TO WHATEVER'S APPROPRIATE

    time = 0
    do while (time<t_final)

        if ((time-int(time))<dt) then
            write(*,*) time, int(time)
        end if 
        
    time = time + dt

    end do

    write(*,*) time, int(time), 0.1*1e-5

end program

