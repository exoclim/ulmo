!**********************************************************************
! Module contains calculations of  the derivative of A wrt a coordinate
!**********************************************************************
module dA_da
use DEGREE_TO_RADIAN
use Constants
use, intrinsic :: iso_fortran_env
implicit none
public :: calculate_dA_d_theta, calculate_dA_d_phi
private

contains
!********************************************************************************************************************
! Subroutines to Calculate the first spatial derivative of A (T or M) with respect to a (theta or phi), on a spherical surface.
!********************************************************************************************************************
! Calculate derivative of mass transport components mass_flux_theta and mass_flux_phi. Boundary conditions are used
! for a Spherical geometry, using A(i,1) to be adjacent to A(i,n_lons-1), and also A(1,j) is adjacent to the opposite point in longitude (j+-n_lats/2)
! at the same latitude.
subroutine calculate_dA_d_theta(A,i,j,dA_theta_dtheta)

    integer(int64), intent(in) :: i,j
    real(real64) :: d_theta
    integer(int64) :: j_temp
    real(real64), dimension(:,:), intent(in) :: A
    real(real64), intent(out) :: dA_theta_dtheta

    d_theta = deg_to_rad(DELTA_LON)

    if (i == N_LATS) then
        j_temp = j+n_lons/2
        if (j_temp > N_LONS) then
            j_temp = j_temp-N_LONS
        end if
        dA_theta_dtheta = (A(i,j_temp)-A(i-1,j))/(2*d_theta)
    elseif (i == 1) then
        j_temp = j+N_LONS/2
        if (j_temp > N_LONS) then
            j_temp = j_temp - N_LONS
        end if
        dA_theta_dtheta = (A(i+1,j)-A(i,j_temp))/(2*d_theta)
    else
        dA_theta_dtheta = (A(i+1,j)-A(i-1,j))/(2*d_theta)
    end if

end subroutine calculate_dA_d_theta

subroutine calculate_dA_d_phi(A,i,j,dA_phi_dphi)

    real(real64) :: d_phi
    real(real64), intent(out) :: dA_phi_dphi
    integer(int64), intent(in) :: i,j
    real(real64), dimension(:,:), intent(in) :: A

    d_phi = deg_to_rad(DELTA_LON)

    if (j == N_LONS) then
        dA_phi_dphi = (A(i,1)-A(i,N_LONS))/(2*d_phi)
    elseif (j == 1) then
        dA_phi_dphi = (A(i,j+1)-A(i,N_LONS))/(2*d_phi)
    else
        dA_phi_dphi = (A(i,j+1)-A(i,j-1))/(2*d_phi)
    end if

end subroutine calculate_dA_d_phi


end module dA_da
