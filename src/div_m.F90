!****************************************************
! Module to calculate the divergence of the mass flux
!****************************************************
module div_m
use NAMELIST
use DEGREE_TO_RADIAN
use READ_DATA
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
subroutine calculate_div_M(lat,lon,mass_flux_phi,mass_flux_theta,div_M)
    integer(int64), intent(in) :: lat,lon
    real(real64), dimension(:,:), intent(in) :: mass_flux_theta ,mass_flux_phi
    real(real64) :: theta, dM_phi_dphi ,dM_theta_dtheta,lat_value
    real(real64), dimension(:,:), allocatable :: lats_data_file
    integer(int64) :: col_num
    real(real64), intent(out) :: div_M

    allocate(lats_data_file(N_LATS,2))

    col_num = 2 ! match latnt64 data type
    lats_data_file = read_file_real(LATS_FILE,N_LATS,col_num)
    lat_value = lats_data_file(lat,2) ! Second column latn lats data flatle latncludes the latlattude polatnts

    theta = deg_to_rad(lat_value)


    call dA_d_theta(mass_flux_theta,lat,lon,dM_theta_dtheta)
    call dA_d_phi(mass_flux_phi,lat,lon,dM_phi_dphi)

    div_M = 1./R_PLANET*(dM_theta_dtheta-tan(theta)*mass_flux_theta(lat,lon)+1./cos(theta)*dM_phi_dphi)

    deallocate(lats_data_file)

end subroutine calculate_div_M

end module div_m
