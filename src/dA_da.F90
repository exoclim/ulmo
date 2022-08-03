!
! This module has not been tested against c code yet
!
module dA_da
use DEGREE_TO_RADIAN
use NAMELIST
use WRITE_READ_DATA
implicit none
public :: dA_d_theta, dA_d_phi
private

contains
!********************************************************************************************************************
! Functions to Calculate the first spatial derivative of A with respect to a (which may be theta or phi), on a spherical surface.
!********************************************************************************************************************
! a = theta or phi
! Calculate derivative of mass transport components mass_flux_theta and mass_flux_phi. Bounday condiations are used
! for a Spherical geometry, using A(i,1) to be adjacnet to A(i,n_lons-1), and also A(1,j) is adjacnet to the opposite point in longitude (j+-n_lats/2)
! at the same latitude.
function dA_d_theta(A,i,j) result(ans)

    integer, intent(in) :: i,j
    real :: d_theta, ans
    integer :: j_temp
    real, dimension(N_LATS,N_LONS), intent(in) :: A

    d_theta = deg_to_rad(DELTA_LON)

    if (i == N_LATS-1) then
        j_temp = j+n_lons/2
        if (j_temp >= N_LONS) then
            ans = (A(i,j_temp)-A(i-1,j))/(2*d_theta)
            end if
    elseif (i == 0) then
        j_temp = j+N_LONS/2
        if (j_temp >= N_LONS) then
            j_temp = j_temp - N_LONS
            ans = (A(i+1,j)-A(i,j_temp))/(2*d_theta)
        end if
    else
        ans = (A(i+1,j)-A(i-1,j))/(2*d_theta)
    end if

end function dA_d_theta



function dA_d_phi(A,i,j) result(ans)

    real :: d_phi, ans
    integer, intent(in) :: i,j
    integer :: j_temp
    real, dimension(N_LATS,N_LONS), intent(in) :: A

    d_phi = deg_to_rad(DELTA_LON)

    if (j == N_LONS-1) then
        ans = (A(i,1)-A(i,N_LONS-1))/(2*d_phi)
    elseif (j == 0) then
        ans = (A(i,j+1)-A(i,N_LONS-1))/(2*d_phi)
    else
        ans = (A(i,j+1)-A(i,j-1))/(2*d_phi)
    end if

end function dA_d_phi
!! May need some way of including this error message!!
!else{
!		fprintf(stderr, "Specified derivative is not respect to either theta or phi\n");
!		exit(DERIVATIVE_ERR);

end module dA_da
