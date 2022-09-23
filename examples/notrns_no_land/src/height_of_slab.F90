!********************************************
! This module contains depth calaculations
!********************************************
module HEIGHT_OF_SLAB
use Constants
use, intrinsic :: iso_fortran_env
implicit none
public:: h_slab
private
contains
!****************************************************
!Function returns the height of a layer given index h
!****************************************************
subroutine h_slab(h,thickness)
    integer(real64), intent(in) :: h
    real(real64),intent(out) :: thickness
    if (h == 1) then !using index h = 1 for surface
        thickness = H_S
    elseif (h == 2) then ! using index h = 2 for deep
        thickness = H_D
    else
        print "(a)",'Index for depth is out of range'

    end if
end subroutine h_slab

end module HEIGHT_OF_SLAB
