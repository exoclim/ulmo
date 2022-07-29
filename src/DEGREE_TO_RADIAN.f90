!*************************************************************
! This module takes a angle in degrees and converts to radians
!*************************************************************
module DEGREE_TO_RADIAN
use NAMELIST
implicit none
public :: deg_to_rad
private

contains
!*********************************************
! Function convert angle in degrees to radians
!*********************************************
function deg_to_rad(angle) result(ans)
    real, intent(in) :: angle
    real :: ans
    ans = angle*pi/180
end function deg_to_rad


end module DEGREE_TO_RADIAN
