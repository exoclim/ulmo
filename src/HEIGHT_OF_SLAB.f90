!********************************************
! Returns the height of a layer given index h
!********************************************
module HEIGHT_OF_SLAB
use NAMELIST
implicit none
public:: h_slab
private
contains
function h_slab(h) result(ans)
    integer :: h
    real*8 :: ans
    if (h == 1) then !using index 1 for surface
        ans = H_S
    elseif (h == 0) then !using index 0 for deep
        ans = H_D
    else
        print "(a)",'Index for depth is out of range'
        stop ! I think this stop statement will work


    end if
end function
end module HEIGHT_OF_SLAB