!*********************************************************
! Module to calculate the derivative of A wrt a coordinate
!*********************************************************
module dA_da
use DEGREE_TO_RADIAN
use NAMELIST
use READ_DATA
use, intrinsic :: iso_fortran_env
implicit none
public :: dA_d_theta, dA_d_phi
private

contains
!********************************************************************************************************************
! Subroutines to Calculate the first spatial derivative of A with respect to a (which may be theta or phi), on a spherical surface.
!********************************************************************************************************************
! a = theta or phi
! Calculate derivative of mass transport components mass_flux_theta and mass_flux_phi. Bounday condiations are used
! for a Spherical geometry, using A(i,1) to be adjacent to A(i,n_lons-1), and also A(1,j) is adjacnet to the opposite point in longitude (j+-n_lats/2)
! at the same latitude.
subroutine dA_d_theta(A,i,j,dA_theta_dtheta)

    integer(int64), intent(in) :: i,j
    real(real64) :: d_theta
    integer(int64) :: j_temp
    real(real64), dimension(:,:), intent(in) :: A
    real(real64), intent(out) :: dA_theta_dtheta

    d_theta = deg_to_rad(DELTA_LON)

    if (i == N_LATS) then
        j_temp = j+n_lons/2
        if (j_temp >= N_LONS) then
            j_temp = j_temp-N_LONS
        end if
        dA_theta_dtheta = (A(i,j_temp)-A(i-1,j))/(2*d_theta)
    elseif (i == 1) then
        j_temp = j+N_LONS/2
        if (j_temp >= N_LONS) then
            j_temp = j_temp - N_LONS
        end if
        dA_theta_dtheta = (A(i+1,j)-A(i,j_temp))/(2*d_theta)
    else
        dA_theta_dtheta = (A(i+1,j)-A(i-1,j))/(2*d_theta)
    end if

end subroutine dA_d_theta

subroutine dA_d_phi(A,i,j,dA_phi_dphi)

    real(real64) :: d_phi
    real(real64), intent(out) :: dA_phi_dphi
    integer(int64), intent(in) :: i,j
    real(real64), dimension(:,:), intent(in) :: A

    d_phi = deg_to_rad(DELTA_LON)

    if (j == N_LONS) then
        dA_phi_dphi = (A(i,1)-A(i,N_LONS-1))/(2*d_phi)
    elseif (j == 1) then
        dA_phi_dphi = (A(i,j+1)-A(i,N_LONS-1))/(2*d_phi)
    else
        dA_phi_dphi = (A(i,j+1)-A(i,j-1))/(2*d_phi)
    end if

end subroutine dA_d_phi


end module dA_da
