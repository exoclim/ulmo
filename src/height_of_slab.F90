!********************************************
! Returns the height of a layer given index h
!********************************************
module HEIGHT_OF_SLAB
use Constants
use fgsl
use, intrinsic :: iso_fortran_env
implicit none
public:: h_slab
private
contains
function h_slab(h) result(ans)
    integer(real64), intent(in) :: h
    real(real64) :: ans
    if (h == 1) then !using index 1 for surface
        ans = H_S
    elseif (h == 2) then ! using index 2 for deep
        ans = H_D
    else
        print "(a)",'Index for depth is out of range'
        stop ! I think this stop statement will work


    end if
end function
end module HEIGHT_OF_SLAB
