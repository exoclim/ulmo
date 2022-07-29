!**********************************************************
! This module contains all the Mass Flux functions for UlMO
!**********************************************************
module MASS_FLUX
use TRANSPORT
use WRITE_READ_DATA
use DEGREE_TO_RADIAN
use NAMELIST
implicit none
public :: calculate_flow_sv_THETA  , calculate_mass_flux_THETA, calculate_mass_flux_PHI, calculate_flow_sv_PHI
private
!
contains
!**********************************************************
! Function to calculate the surface stress in PHI direction
!**********************************************************
function calculate_surface_stress_PHI() result(tau_PHI)
     real*8, dimension(n_lats,n_lons) :: u_wind
     real*8, dimension(n_lats,n_lons) :: tau_PHI
     integer :: i,j

     u_wind = read_file(U_WIND_DATA,N_LATS,N_LONS)

!     do coord !! not sure how this works / what output do i need? !!
        do i = 1,N_LATS ! over lats points
            do j = 1,N_LONS ! over lons points
                tau_PHI(i,j) = C_D*RHO_AIR*u_wind(i,j)*abs(u_wind(i,j))
                end do
            end do
 !       end do


end function calculate_surface_stress_PHI
!************************************************************
! Function to calculate the surface stress in THETA direction
!************************************************************
function calculate_surface_stress_THETA() result(tau_THETA)
     real*8, dimension(n_lats,n_lons) :: v_wind
     real*8, dimension(n_lats,n_lons) :: tau_THETA
     integer :: i,j

     v_wind = read_file(V_WIND_DATA,N_LATS,N_LONS)

        do i = 1,N_LATS ! over lats points
            do j = 1,N_LONS ! over lons points
                tau_THETA(i,j) = C_D*RHO_AIR*v_wind(i,j)*abs(v_wind(i,j))
                end do
            end do


end function calculate_surface_stress_THETA
!*********************************************
! Function to calculate the Coriolis parameter
!*********************************************
function calculate_coriolis_parameter(theta_deg) result(f_out)
    real*8 :: f_out
    real*8, intent(in)  :: theta_deg
    f_out = 2*OMEGA*sin(deg_to_rad(theta_deg))

end function calculate_coriolis_parameter
!*************************************************************************************
! Functions to calculate the mass flux in THETA and PHI directions from surface stress
!*************************************************************************************
function calculate_mass_flux_THETA() result(mass_flux_THETA)
    integer :: i,j
    real*8, dimension(N_LATS,N_LONS) :: tau_PHI
    real*8, dimension(N_LATS,N_LONS) :: tau_THETA
    real*8, dimension(N_LATS,N_lONS) :: mass_flux_THETA
    real*8, dimension(N_LATS,1)      :: f
    real*8, dimension(N_LATS,1)      :: lats
    real*8, dimension(N_LATS,2)      :: lats_data_file
    real*8, dimension(N_LATS,N_LONS) :: land_mask


    tau_PHI = calculate_surface_stress_PHI()
    tau_THETA = calculate_surface_stress_THETA()
    lats_data_file = read_file(LATS_FILE,N_LATS,2)
    lats(:,1)= lats_data_file(:,2) ! Second column in lats data file includes the latitude points
    land_mask = read_file(LAND_MASK_DATA,N_LATS,N_LONS)


    do i = 1,N_LATS

        f(i,1) = calculate_coriolis_parameter(lats(i,1))

        do j = 1,N_LONS

            !! land mask if statement !!
            if (land_mask(i,j) == 1) then
                mass_flux_THETA(i,j) = 0
            else
                mass_flux_THETA(i,j) = (EPSILON*tau_THETA(i,j)+f(i,1)*tau_PHI(i,j))/(EPSILON**2.0+f(i,1)**2)
            end if


        end do
    end do

end function calculate_mass_flux_THETA

function calculate_mass_flux_PHI() result(mass_flux_PHI)
    integer :: i,j
    real*8, dimension(N_LATS,N_LONS) :: tau_PHI
    real*8, dimension(N_LATS,N_LONS) :: tau_THETA
    real*8, dimension(N_LATS,N_lONS) :: mass_flux_PHI
    real*8, dimension(N_LATS,1)      :: f
    real*8, dimension(N_LATS,1)      :: lats
    real*8, dimension(N_LATS,2)      :: lats_data_file
    real*8, dimension(N_LATS,N_LONS) :: land_mask


    tau_PHI = calculate_surface_stress_PHI()
    tau_THETA = calculate_surface_stress_THETA()
    lats_data_file = read_file(LATS_FILE,N_LATS,2)
    lats(:,1)= lats_data_file(:,2) ! Second column in lats data file includes the latitude points
    land_mask = read_file(LAND_MASK_DATA,N_LATS,N_LONS)


    do i = 1,N_LATS

        f(i,1) = calculate_coriolis_parameter(lats(i,1))

        do j = 1,N_LONS

            !! land mask if statement !!
            if (land_mask(i,j) == 1) then
                mass_flux_PHI(i,j) = 0
            else
                mass_flux_PHI(i,j) = (EPSILON*tau_PHI(i,j)-f(i,1)*tau_THETA(i,j))/(EPSILON**2.0+f(i,1)**2.0)
            end if

        end do
    end do

end function calculate_mass_flux_PHI
!!****************************************************************************
!! Functions to convert mass flux into the unit of Sverdrups for Theta and Phi
!!****************************************************************************
function calculate_flow_sv_PHI() result(sv_flow_PHI)
    integer :: i,j
    real*8, dimension(N_LATS,N_lONS) :: mass_flux_PHI
    real*8, dimension(N_LATS,N_lONS) :: sv_flow_PHI

    mass_flux_PHI = calculate_mass_flux_PHI()

    do i = 1,N_LATS
        do j= 1,N_LONS
            sv_flow_PHI = mass_flux_PHI * R_PLANET * deg_to_rad(DELTA_LAT)  / RHO_WATER /1.0e6

        end do
    end do

end function calculate_flow_sv_PHI

function calculate_flow_sv_THETA() result(sv_flow_THETA)
    integer :: i,j
    real*8, dimension(N_LATS,2)      :: lats_data_file
    real*8, dimension(N_LATS,1)      ::  lats
    real*8, dimension(N_LATS,1)      :: theta
    real*8, dimension(N_LATS,N_lONS) :: mass_flux_THETA
    real*8, dimension(N_LATS,N_lONS) :: sv_flow_THETA


    lats_data_file = read_file(LATS_FILE,N_LATS,2)
    lats(:,1)= lats_data_file(:,2) ! Second column in lats data file includes the latitude points

    mass_flux_THETA = calculate_mass_flux_THETA()

    do i = 1,N_LATS
        theta(i,1) = deg_to_rad(lats(i,1))
        do j= 1,N_LONS
            sv_flow_THETA(i,j) = mass_flux_THETA(i,j) * R_PLANET * cos(theta(i,1)) * deg_to_rad(DELTA_LON) / RHO_WATER /1.0e6

        end do
    end do

end function calculate_flow_sv_THETA



end module MASS_FLUX
