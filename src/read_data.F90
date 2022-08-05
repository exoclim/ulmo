module READ_DATA
use NAMELIST
implicit none
!********************************************************************
! This module also contains the original subroutine for writing files
!********************************************************************
public :: read_file
private

contains
!**************************
! Function to read in data
!**************************
function read_file(data_file,rows,cols) result(array)
    character(len=*), intent(in) :: data_file
    integer, intent(in) :: rows,cols

    real, dimension(rows,cols) :: array

    integer :: lat,lon,iu



    open(newunit=iu,file=data_file,status='old',action='read')
    do lat = 1,rows
        read(iu,*) (array(lat,lon),lon=1,cols) ! reading in input data file

    end do
    close(iu)

!    print*,data_file,'(1,1)', output(1,1) !!check!!

end function read_file
!**************************************************************
! Subroutine to write output arrays to data file (old/not used)
!**************************************************************
subroutine write_file(file_name,array,rows, cols)
    integer, intent(in) :: rows, cols
    real, intent(in),dimension(rows,cols) :: array
    integer :: i,j,iu
    character(len=*), intent(in) :: file_name

    open(newunit=iu, file = file_name, status = 'replace', action='write')
    print "(a,i5)", 'Writing output...'
    do i = 1,rows
        write(iu,*) (array(i,j),j=1,cols)
    end do

    close(iu)



end subroutine write_file



end module READ_DATA
