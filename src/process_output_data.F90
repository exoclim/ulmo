!**********************************************************
! This module processes the output data and writes to files
!**********************************************************
module process_output_data
 use Constants
 use READ_DATA
 use heat_fluxes
 use, intrinsic :: iso_fortran_env
 implicit none
 public :: process_output
 private
 contains
!********************************************************************
! Function constructs output file name from base name and output time
!********************************************************************
function construct_fname(base,time_str) result(fname)
    character(len=100) :: fname
    character(len=50), intent(in) :: base
    character(len=6), intent(in) :: time_str
    fname = trim(base)//trim(time_str)//OUTPUT_DATA_EXT
    !print*,'fname',fname
end function
!****************************************************************************
! Subroutine that writes output data to unique file name including time stamp
!****************************************************************************
subroutine write_file_time(base,time_str,array,rows,cols)
    integer(int64), intent(in) :: rows, cols
    real(real64), intent(in),dimension(:,:) :: array
    integer(int64) :: i,j,iu
    character(len=50), intent(in) :: base
    character(len=6), intent(in) :: time_str
    character(len=100) :: file_name

    file_name = construct_fname(base,time_str)
    !print*, file_name

    open(newunit=iu, file = file_name, status = 'replace', action='write')
    do i = 1,rows
        write(iu,*) (array(i,j),j=1,cols)
    end do

    close(iu)

end subroutine write_file_time
!*********************************************
! Function that converts a integer to a string
!*********************************************
function int_to_str(i) result(res)
  character(len=6) :: res
  integer(int64),intent(in) :: i
  character(range(i)+2) :: tmp
  write(tmp,'(i0)') i
  res = trim(tmp)
end function

!**********************************************************
! Subroutine to Process the output data for each time stamp
!**********************************************************
 subroutine process_output(T,upward_Q_flux,time)!,time)!,M)

    real(real64),dimension(:,:,:),intent(in) :: T !,M
    !real(real64),dimension(2,N_LATS,N_LONS) :: Sv_flow
    real(real64),dimension(:,:),intent(in) :: upward_Q_flux !,Sv_flow_PHI,Sv_flow_THETA,div_M
    real(real64) :: days
    integer(int64) :: days_int
    character(len=6) :: time_str
    real(real64), intent(in) :: time


!    days = T_OFFSET+time/(HOURS_PER_DAY*MINUTES_PER_HOUR*SECONDS_PER_MINUTE)
!    print*,'days = ',days
!    days_int = nint(days) ! rounding to nearest integer from real
!    print*,'days_int',days_int
!    print*,"Data outputted at:",days_int
!    !day_int = days_int+TIME_OUTPUT_TOL ! not sure about this !
!    time_str = int_to_str(day_int)
!    print*, 'time_str=',time_str
    days = T_OFFSET+time/(HOURS_PER_DAY*MINUTES_PER_HOUR*SECONDS_PER_MINUTE)
    days = nint(days)
    print*, 'Data outputted at = ', days
    days_int = int(days)
    time_str = int_to_str(days_int)

    !**Temperature outputs**!
    call write_file_time(OUTPUT_DEEP_TEMP_DATA,time_str,T(2,:,:),N_LATS,N_LONS)
    call write_file_time(OUTPUT_SURFACE_TEMP_DATA,time_str,T(1,:,:),N_LATS,N_LONS)
    call write_file_time(OUTPUT_UPWARD_Q_FLUX,time_str,upward_Q_flux,N_LATS,N_LONS)

    !**Mass Flux outputs**!
    !call write_file_time(OUTPUT_X_MASS_FLUX_DATA,time_str,Sv_flow_PHI,N_LATS,N_LONS)
    !call write_file_time(OUTPUT_Y_MASS_FLUX_DATA,time_str,Sv_flow_THETA,N_LATS,N_LONS)

    !do i = 1,N_LATS
    !    do j = 1,N_LONS
    !        div_M(i,j) = calculate_div_M(i,j)
    !    end do
    !end do

    !call write_file_time(OUTPUT_VERITCAL_FLUX_DATA,time_str,div_M,N_LATS,N_LONS)

!    call write_file('output_data/ProCb/OUTPUT_UPWARD_Q_FLUX.dat',upward_Q_flux,N_LATS,N_LONS)
!    call write_file('output_data/ProCb/OUTPUT_SURFACE_TEMP_DATA.dat',T(1,:,:),N_LATS,N_LONS)
!    call write_file('output_data/ProCb/OUTPUT_DEEP_TEMP_DATA.dat',T(2,:,:),N_LATS,N_LONS)



 end subroutine
end module
