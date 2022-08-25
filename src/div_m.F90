!****************************************************
! Module to calculate the divergence of the mass flux
!****************************************************
module div_m
use Constants
use DEGREE_TO_RADIAN
use dA_da
use MASS_FLUX
use, intrinsic :: iso_fortran_env
implicit none
public :: calculate_div_M
private
contains
!***********************************************************************************
!Subroutine returns the horizontal divergence of quanitity M. determines whether there is an up
!or downward flow between the ocean layers
!************************************************************************************
subroutine calculate_div_M(lat,lon,lats,mass_flux_phi,mass_flux_theta,div_M)
    integer(int64), intent(in) :: lat,lon
    real(real64), dimension(:,:), intent(in) :: mass_flux_theta ,mass_flux_phi
    real(real64) :: theta, dM_phi_dphi ,dM_theta_dtheta,lat_value
    real(real64), dimension(:,:),intent(in) :: lats ! This can be converted to a grid input
    real(real64), intent(out) :: div_M

    lat_value = lats(lat,1) ! This can be converted to a grid input
    theta = deg_to_rad(lat_value)

    call dA_d_theta(mass_flux_theta,lat,lon,dM_theta_dtheta)
    call dA_d_phi(mass_flux_phi,lat,lon,dM_phi_dphi)

    div_M = 1./R_PLANET*(dM_theta_dtheta-tan(theta)*mass_flux_theta(lat,lon)+1./cos(theta)*dM_phi_dphi)

end subroutine calculate_div_M

end module div_m
