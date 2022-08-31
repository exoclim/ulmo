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
subroutine calculate_div_M(lat,lon,theta,mass_flux_phi,mass_flux_theta,div_M)
    integer(int64), intent(in) :: lat,lon
    real(real64), intent(in) :: theta
    real(real64), dimension(:,:), intent(in) :: mass_flux_theta ,mass_flux_phi
    real(real64) ::  dM_phi_dphi ,dM_theta_dtheta
    real(real64), intent(out) :: div_M



    call calculate_dA_d_theta(mass_flux_theta,lat,lon,dM_theta_dtheta)
    call calculate_dA_d_phi(mass_flux_phi,lat,lon,dM_phi_dphi)

    div_M = 1./R_PLANET*(dM_theta_dtheta-tan(theta)*mass_flux_theta(lat,lon)+1./cos(theta)*dM_phi_dphi)

end subroutine calculate_div_M

end module div_m
