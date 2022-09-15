!**********************************************************
! This module contains all the Mass Flux functions for UlMO
!**********************************************************
module MASS_FLUX
use DEGREE_TO_RADIAN, only: deg_to_rad
use Constants
use, intrinsic :: iso_fortran_env
implicit none
public :: calculate_surface_stress_PHI, calculate_surface_stress_THETA,calculate_flow_sv_THETA , calculate_mass_flux_THETA &
            , calculate_mass_flux_PHI, calculate_flow_sv_PHI
private

contains
!**********************************************************
! Subroutine to calculate the surface stress in PHI direction
!**********************************************************
subroutine calculate_surface_stress_PHI(u_wind,tau_PHI)
    real(real64), dimension(:,:),intent(in) :: u_wind
    real(real64), dimension(:,:),intent(out) :: tau_PHI
    integer(int64) :: i,j

    do i = 1,N_LATS
        do j = 1,N_LONS
            tau_PHI(i,j) = C_D*RHO_AIR*u_wind(i,j)*abs(u_wind(i,j))
        end do
    end do


end subroutine calculate_surface_stress_PHI
!************************************************************
! Subroutine  to calculate the surface stress in THETA direction
!************************************************************
subroutine calculate_surface_stress_THETA(v_wind,tau_THETA)
    real(real64), dimension(:,:),intent(in):: v_wind(:,:)
    real(real64), dimension(:,:),intent(out) :: tau_THETA(:,:)
    integer(int64) :: i,j
    do i = 1,N_LATS
        do j = 1,N_LONS
            tau_THETA(i,j) = C_D*RHO_AIR*v_wind(i,j)*abs(v_wind(i,j))
        end do
    end do

end subroutine calculate_surface_stress_THETA
!*********************************************
! Subroutine  to calculate the Coriolis parameter
!*********************************************
function calculate_coriolis_parameter(theta_deg) result(ans)
    real(real64) :: ans
    real(real64), intent(in)  :: theta_deg
    ans = 2*OMEGA*sin(deg_to_rad(theta_deg))

end function calculate_coriolis_parameter
!*************************************************************************************
! Subroutines to calculate the mass flux in THETA and PHI directions from surface stress
!*************************************************************************************
subroutine calculate_mass_flux_THETA(lats,tau_PHI,tau_THETA,land_mask,mass_flux_THETA)
    integer(int64) :: i,j
    real(real64) :: f
    integer(int64), dimension(:,:),intent(in) :: land_mask
    real(real64), dimension(:,:),intent(in) :: tau_PHI,tau_THETA,lats
    real(real64), dimension(:,:),intent(out) :: mass_flux_THETA

    do i = 1,N_LATS
        f = calculate_coriolis_parameter(lats(i,1))
        do j = 1,N_LONS
            if (land_mask(i,j) == 1) then
                mass_flux_THETA(i,j) = 0
            else
                mass_flux_THETA(i,j) = (EPSI*tau_THETA(i,j)+f*tau_PHI(i,j))/(EPSI**2.0+f**2)
            end if
        end do
    end do

end subroutine calculate_mass_flux_THETA

subroutine calculate_mass_flux_PHI(lats,tau_PHI,tau_THETA,land_mask,mass_flux_PHI)
    integer(int64) :: i,j
    real(real64) :: f
    integer(int64), dimension(:,:),intent(in) :: land_mask
    real(real64), dimension(:,:),intent(in) :: tau_PHI,tau_THETA,lats
    real(real64), dimension(:,:),intent(out) :: mass_flux_PHI

    do i = 1,N_LATS
        f = calculate_coriolis_parameter(lats(i,1))
        do j = 1,N_LONS
            !! land mask if statement !!
            if (land_mask(i,j) == 1) then
                mass_flux_PHI(i,j) = 0
            else
                mass_flux_PHI(i,j) = (EPSI*tau_PHI(i,j)-f*tau_THETA(i,j))/(EPSI**2.0+f**2.0)
            end if
        end do
    end do

end subroutine calculate_mass_flux_PHI
!****************************************************************************
! Subroutines to convert mass flux into the unit of Sverdrups for Theta and Phi
!****************************************************************************
subroutine calculate_flow_sv_PHI(mass_flux_PHI,sv_flow_PHI)
    integer(int64) :: i,j
    real(real64), dimension(:,:),intent(out) :: sv_flow_PHI
    real(real64),dimension(:,:),intent(in) :: mass_flux_PHI

    do i = 1,N_LATS
        do j= 1,N_LONS
            sv_flow_PHI = mass_flux_PHI * R_PLANET * deg_to_rad(DELTA_LAT)  / RHO_WATER /1.0e6
        end do
    end do

end subroutine calculate_flow_sv_PHI

subroutine calculate_flow_sv_THETA(lats,mass_flux_THETA,sv_flow_THETA)
    integer(int64) :: i,j
    real(real64) :: theta
    real(real64), dimension(:,:),intent(in)  :: mass_flux_THETA,lats
    real(real64), dimension(:,:),intent(out) :: sv_flow_THETA

    do i = 1,N_LATS
        theta = deg_to_rad(lats(i,1))
        do j= 1,N_LONS
            sv_flow_THETA(i,j) = mass_flux_THETA(i,j) * R_PLANET * cos(theta) * deg_to_rad(DELTA_LON) / RHO_WATER /1.0e6
        end do
    end do


end subroutine calculate_flow_sv_THETA



end module MASS_FLUX
