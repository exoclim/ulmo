!*************************************************************
! This module takes a angle in degrees and converts to radians
!*************************************************************
module DEGREE_TO_RADIAN
use Constants
use, intrinsic :: iso_fortran_env
implicit none
public :: deg_to_rad
private

contains
!*********************************************
! Function convert angle in degrees to radians
!*********************************************
function deg_to_rad(angle) result(ans)
    real(real64), intent(in) :: angle
    real(real64) :: ans
    ans = angle*pi/180
end function deg_to_rad


end module DEGREE_TO_RADIAN
