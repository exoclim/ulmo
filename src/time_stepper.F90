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
    !real(int64), dimension(2,N_LATS,N_LONS) :: T_init!,M
    real(real64),dimension(2,N_LATS,N_LONS) :: T
    real(real64), dimension(N_LATS,N_LONS) :: T_surf,T_deep
    integer(int64), intent(in) :: n_times !, version
    real(real64) :: time, days !,avg_T_surf_init,avg_T_surf_new
    integer(int64) :: ti

    T_surf= read_file(INITIAL_SURFACE_TEMP_DATA,N_LATS,N_LONS)
    T_deep= read_file(INITIAL_DEEP_TEMP_DATA,N_LATS,N_LONS)

    T(1,:,:) = T_surf
    T(2,:,:) = T_deep


    !Grid_vals *grid = xmalloc(sizeof(Grid_vals)) I think thi is unique to c
    ! need construct grid statement
    !temperature data
!    time = 0
!    days = 0
!    allocate(T(2,N_LATS,N_LONS))
!    do ti = 1,n_times
!        time = (ti+1.)*DELTA_T
!        print*, T(1,1,1)
!        call calculate_new_T(T)
!        if (mod(time,TIME_OUTPUT_FREQ)<TIME_OUTPUT_TOL) then
!            days = T_OFFSET+time/(HOURS_PER_DAY*MINUTES_PER_HOUR*SECONDS_PER_MINUTE)
!            print*, "Days passed = %lg\n", days
!        end if
!
!        if (mod(time,DATA_OUTPUT_FREQ)<TIME_OUTPUT_TOL) then
!            call process_output(T,time)!,M)
!        end if
!
!        !M_THETA = calculate_mass_flux_THETA()
!        !M_PHI = calculate_mass_flux_PHI()
!        !M(1,:,:) = mass_flux_THETA
!        !M(2,:,:) = mass_flux_PHI:
!
!        !allocate(T(2,N_LATS,N_LONS))
!
!!
!!        if (ti == 1) then
!!
!!
!!            print*,'avg_T_surf_init',time,'=', sum(T_surf_init)/12960
!!            allocate(T_new(2,N_LATS,N_LONS))
!!
!!            call calculate_new_T(T_init,T_new)
!!
!!            call process_output(T_new,time)
!!
!!            print*, 'avg_T_surf_0',time,'=', sum(T_new(1,:,:))/12960
!!
!!            T = T_new
!!            !print*,T(1,1,1)
!!            !print*,T_new(1,1,1)
!!
!!            deallocate(T_new)
!!
!!
!!        else
!!!            print*, T_new(1,1,1) ! gives me right error
!!            print*, T(1,1,1)
!!
!!            allocate(T_new(2,N_LATS,N_LONS))
!!            call calculate_new_T(T,T_new)
!!            print*,T_new(1,1,1)
!!
!!            !print*,'avg T surf',time,sum(T_new(1,:,:))/12960
!!
!!            if (mod(time,TIME_OUTPUT_FREQ)<TIME_OUTPUT_TOL) then
!!                days = T_OFFSET+time/(HOURS_PER_DAY*MINUTES_PER_HOUR*SECONDS_PER_MINUTE)
!!                print*, "Days passed = %lg\n", days
!!            end if
!!
!!            if (mod(time,DATA_OUTPUT_FREQ)<TIME_OUTPUT_TOL) then
!!                call process_output(T_new,time)!,M)
!!            end if
!!
!!            T = T_new
!!            deallocate(T_new)
!!
!!
!!
!!        end if
!
!    !deallocate(T)
!
!
!
!
!
!
!    end do

!    deallocate(T)


end subroutine

end module
