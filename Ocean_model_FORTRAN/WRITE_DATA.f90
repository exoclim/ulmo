module WRITE_DATA
use NAMELIST
implicit none

public :: write_file
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
    do i = 1,rows
        write(iu,*) (array(i,j),j=1,cols)
    end do

    close(iu)



end subroutine write_file

end module WRITE_DATA
