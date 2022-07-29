module TIME_STEP
use NAMELIST
use TRANSPORT
use WRITE_READ_DATA
implicit none
public :: time_stepper,temp_after_time_step
private

contains
!******************************************
! Steps the transport equation through time
!******************************************
function time_stepper(n_times) result(output)
    integer, intent(in) :: n_times
    real*8, dimension(N_LATS,N_LONS) :: T_init_surf, T_init_deep,output
    integer :: time
    T_init_surf = read_file(INITIAL_SURFACE_TEMP_DATA,N_LATS,N_LONS)
    T_init_deep = read_file(INITIAL_DEEP_TEMP_DATA,N_LATS,N_LONS)

    do time = 1,n_times
        output = no_transport_no_deep(T_init_surf)
    end do
end function time_stepper
!*******************************************************
! Writes temperature output to files 'surf_temp_out.dat'
!*******************************************************
subroutine temp_after_time_step(N_STEPS)
    integer, intent(in):: N_STEPS
    real*8, dimension(N_LATS,N_LONS) :: output
    output = time_stepper(N_STEPS)

    call write_file('surf_temp_out.dat',output,N_LATS,N_LONS)


end subroutine

end module TIME_STEP
