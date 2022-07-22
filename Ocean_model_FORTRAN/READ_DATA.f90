!************************************************
! This module reads in  data and outputs an array
!************************************************
module READ_DATA
implicit none
public :: read_file
private
!
contains
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

end module READ_DATA
