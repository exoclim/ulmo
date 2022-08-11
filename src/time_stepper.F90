module time_stepper
use NAMELIST
use fgsl
use MASS_FLUX
use MATRIX_CALC
use READ_DATA
use process_output_data
use, intrinsic :: iso_fortran_env
implicit none
public :: t_stepper
private

contains
!***********************************************************************************
! Subroutine that steps through selected time. First creates the horizontal grid,
! then read in temperatures, and steps through n_times time steps. Writes the final
! Temperature data to files, and plots a selected plotting script
!***********************************************************************************
subroutine t_stepper(n_times)!,version)
    !integer(int64), dimension(N_LATS,N_LONS), intent(in) :: land_mask
    !real(int64), dimension(N_LATS,N_LONS) :: M_THETA,M_PHI
    real(int64), dimension(2,N_LATS,N_LONS) :: T,T_new!,M
    real(real64), dimension(N_LATS,N_LONS) :: T_surf, T_deep
    integer(int64), intent(in) :: n_times !, version
    real(real64) :: time, days
    integer(int64) :: ti

    T_surf= read_file(INITIAL_SURFACE_TEMP_DATA,N_LATS,N_LONS)
    T_deep = read_file(INITIAL_DEEP_TEMP_DATA,N_LATS,N_LONS)

    T(1,:,:) = T_surf
    T(2,:,:) = T_deep


    !Grid_vals *grid = xmalloc(sizeof(Grid_vals)) I think thi is unique to c
    ! need construct grid statement
    !temperature data
!    time = 0
!    days = 0

    do ti = 1,n_times
        time = (ti+1.)*DELTA_T
        !M_THETA = calculate_mass_flux_THETA()
        !M_PHI = calculate_mass_flux_PHI()
        !M(1,:,:) = mass_flux_THETA
        !M(2,:,:) = mass_flux_PHI:
        call calculate_new_T(T,T_new)

        if (mod(time,TIME_OUTPUT_FREQ)<TIME_OUTPUT_TOL) then
            days = T_OFFSET+time/(HOURS_PER_DAY*MINUTES_PER_HOUR*SECONDS_PER_MINUTE)
            print*, "Days passed = %lg\n", days
        end if
        if (mod(time,DATA_OUTPUT_FREQ)<TIME_OUTPUT_TOL) then
            call process_output(T_new,time)!,M)
        end if
        !T = T_new need something like this to re initialise T as T_new
    end do


end subroutine

end module
