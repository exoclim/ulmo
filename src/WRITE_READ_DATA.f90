module WRITE_READ_DATA
use NAMELIST
implicit none

public :: write_file,read_file
private

contains
!***********************************************
! Subroutine to write output arrays to data file
!***********************************************
subroutine write_file(file_name,array,rows, cols)
    integer, intent(in) :: rows, cols
    real*8, intent(in),dimension(rows,cols) :: array
    integer :: i,j,iu
    character(len=*), intent(in) :: file_name

    open(newunit=iu, file = file_name, status = 'replace', action='write')
    print "(a,i5)", 'Writing output...'
    do i = 1,rows
        write(iu,*) (array(i,j),j=1,cols)
    end do

    close(iu)



end subroutine write_file
!**************************
! function to read in data
!**************************
function read_file(data_file,rows,cols) result(output)
    character(len=*), intent(in) :: data_file
    real*8, dimension(rows,cols) :: output
    real*8, dimension(rows,cols) :: array
    integer, intent(in) :: rows,cols
    integer :: lat,lon,iu



    open(newunit=iu,file=data_file,status='old',action='read')
    do lat = 1,rows
        read(iu,*) (array(lat,lon),lon=1,cols) ! reading in input data file

    end do
    close(iu)
    output = array
!    print*,data_file,'(1,1)', output(1,1) !!check!!

end function read_file

subroutine process_output_data()

end subroutine

end module WRITE_READ_DATA
