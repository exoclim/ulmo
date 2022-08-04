module process_output_data
 use NAMELIST
 use Other_FLUXES
 use MASS_FLUX
 use DIV_M
 use WRITE_READ_DATA
 implicit none
 public :: process_output
 private
 contains

 subroutine process_output(T,M,time)
    !
    ! need to add time string to write file
    !
    real,dimension(2,N_LATS,N_LONS),intent(in) :: T,M
    real,dimension(2,N_LATS,N_LONS) :: Sv_flow
    real,dimension(N_LATS,N_LONS) :: upward_Q_flux,Sv_flow_PHI,Sv_flow_THETA,div_M
    real, intent(in) :: time
    real :: days
    integer :: day_int, i , j
    character(len=10) :: time_str

    days = T_OFFSET+time/(HOURS_PER_DAY*MINUTES_PER_HOUR*SECONDS_PER_MINUTE)
    days = nint(days) ! rounding to nearest integer from real
    print*,"Data outputted at %lg\n",days
    day_int = days+TIME_OUTPUT_TOL
    print*,time_str,"%d",day_int

    upward_Q_flux = calc_Q_flux(T)
    call write_file(OUTPUT_UPWARD_Q_FLUX,upward_Q_flux,N_LATS,N_LONS)

    Sv_flow_PHI = calculate_flow_sv_PHI()
    Sv_flow_THETA = calculate_flow_sv_THETA()
    call write_file(OUTPUT_X_MASS_FLUX_DATA,Sv_flow_PHI,N_LATS,N_LONS)
    call write_file(OUTPUT_Y_MASS_FLUX_DATA,Sv_flow_THETA,N_LATS,N_LONS)

    do i = 1,N_LATS
        do j = 1,N_LONS
            div_M(i,j) = calculate_div_M(i,j)
        end do
    end do

    call write_file(OUTPUT_VERITCAL_FLUX_DATA,div_M,N_LATS,N_LONS)

    call write_file(OUTPUT_SURFACE_TEMP_DATA,T(1,N_LATS,N_LONS),N_LATS,N_LONS)
    call write_file(OUTPUT_DEEP_TEMP_DATA,T(2,N_LATS,N_LONS),N_LATS,N_LONS)





 end subroutine
end module
