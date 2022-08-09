module time_stepper
use NAMELIST
use fgsl
use MATRIX_CALC

use, intrinsic :: iso_fortran_env
implicit none
public :: calculate_matrix
private

contains
!***********************************************************************************
! Subroutine that steps through selected time. First creates the horizontal grid,
! then read in temperatures, and steps through n_times time steps. Writes the final
! Temperature data to files, and plots a selected plotting script
!***********************************************************************************
!subroutine t_stepper(n_times,version,T,M,land_mask)
!    integer(int64), dimension(N_LATS,N_LONS), intent(in) :: land_mask
!    real(int64), dimension(2,N_LATS,N_LONS),intent(in) :: T,M
!    integer(int64), intent(in) :: n_times, version
!    real(real64) :: time, days
!    integer(int64) :: ti
!
!    !Grid_vals *grid = xmalloc(sizeof(Grid_vals)) I think thi is unique to c
!    ! need construct grid statement
!    !temperature data
!
!    do ti = 1,n_times
!        time = (ti+1.)*DELTA_T
!        if (mod(time,TIME_OUTPUT_FREQ)<TIME_OUTPUT_TOL) then
!            days = T_OFFSET+time/(HOURS_PER_DAY*MINUTES_PER_HOUR*SECONDS_PER_MINUTE)
!            print*, "Days passed = %lg\n", days
!        end if
!        if (mod(time,DATA_OUTPUT_FREQ)<TIME_OUTPUT_TOL) then
!            !process_output_data(T,M,grid,N_LATS,N_LONS,time)
!        end if
!    end do
!
!
!end subroutine

end module
