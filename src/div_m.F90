!
!This module has not been tested against c code yet
!
module DIV_M
use NAMELIST

use DEGREE_TO_RADIAN
use WRITE_READ_DATA
use dA_da
use MASS_FLUX
implicit none
public :: calculate_div_M
private
contains
!***********************************************************************************
! This function outputs the horizontal divergence of quantity M. determines whether
! there is an up or downward flow between the ocean layers
!************************************************************************************
function calculate_div_M(i,j) result(div_M)
    real :: theta, div_M , dM_phi_dphi ,dM_theta_dtheta
    real, dimension(N_LATS,2) :: lats_data_file
    real ::  lat
    integer, intent(in) :: i,j
    real, dimension(N_LATS,N_LONS) :: M_Theta ,M_Phi

    lats_data_file = read_file(LATS_FILE,N_LATS,2)
    lat = lats_data_file(i,2) ! Second column in lats data file includes the latitude points

    theta = deg_to_rad(lat)

    M_Theta = calculate_mass_flux_THETA()
    M_Phi = calculate_mass_flux_PHI()


    dM_theta_dtheta = dA_d_theta(M_Theta,i,j)
    dM_phi_dphi = dA_d_phi(M_Phi,i,j)

    div_M = 1./R_PLANET*(dM_theta_dtheta-tan(theta)*M_Theta(i,j)+1./cos(theta)*dM_phi_dphi)

end function calculate_div_M

end module DIV_M
