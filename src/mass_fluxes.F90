!**********************************************************
! This module contains all the Mass Flux functions for UlMO
!**********************************************************
module MASS_FLUX
use READ_DATA
use DEGREE_TO_RADIAN
use NAMELIST
use, intrinsic :: iso_fortran_env
implicit none
public :: calculate_flow_sv_THETA  , calculate_mass_flux_THETA, calculate_mass_flux_PHI, calculate_flow_sv_PHI
private

contains
!**********************************************************
! Function to calculate the surface stress in PHI direction
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
! Function to calculate the surface stress in THETA direction
!************************************************************
subroutine calculate_surface_stress_THETA(v_wind,tau_THETA)
    real(real64), dimension(:,:),intent(in):: v_wind
    real(real64), dimension(:,:),intent(out) :: tau_THETA
    integer(int64) :: i,j
    do i = 1,N_LATS
        do j = 1,N_LONS
            tau_THETA(i,j) = C_D*RHO_AIR*v_wind(i,j)*abs(v_wind(i,j))
        end do
    end do

end subroutine calculate_surface_stress_THETA
!*********************************************
! Function to calculate the Coriolis parameter
!*********************************************
subroutine calculate_coriolis_parameter(theta_deg,f)
    real(real64),intent(out) :: f
    real(real64), intent(in)  :: theta_deg
    f = 2*OMEGA*sin(deg_to_rad(theta_deg))

end subroutine calculate_coriolis_parameter
!*************************************************************************************
! Functions to calculate the mass flux in THETA and PHI directions from surface stress
!*************************************************************************************
subroutine calculate_mass_flux_THETA(u_wind,v_wind,land_mask,mass_flux_THETA)
    integer(int64) :: i,j,col_num
    real(real64) :: f
    real(real64), dimension(:,:),allocatable :: tau_PHI,tau_THETA,lats,lats_data
    integer(int64), dimension(:,:),intent(in) :: land_mask
    real(real64), dimension(:,:),intent(in) :: u_wind,v_wind
    real(real64), dimension(:,:), intent(out) :: mass_flux_THETA

    allocate(tau_PHI(N_LATS,N_LONS))
    allocate(tau_THETA(N_LATS,N_LONS))

    allocate(lats(N_LATS,1))
    allocate(lats_data(N_LATS,2))

    col_num = 2
    call calculate_surface_stress_PHI(u_wind,tau_PHI)
    call calculate_surface_stress_THETA(v_wind,tau_THETA)
    lats_data = read_file_real(LATS_FILE,N_LATS,col_num)
    lats(:,1)= lats_data(:,2) ! Second column in lats data file includes the latitude points

    do i = 1,N_LATS
        call calculate_coriolis_parameter(lats(i,1),f)
        do j = 1,N_LONS
            !! land mask if statement !!
            if (land_mask(i,j) == 1) then
                mass_flux_THETA(i,j) = 0
            else
                mass_flux_THETA(i,j) = (EPSILON*tau_THETA(i,j)+f*tau_PHI(i,j))/(EPSILON**2.0+f**2)
            end if
        end do
    end do

    deallocate(tau_PHI,tau_THETA,lats,lats_data)
!    deallocate(tau_THETA)
!    deallocate(f)
!    deallocate(lats)
!    deallocate(lats_data)

end subroutine calculate_mass_flux_THETA

subroutine calculate_mass_flux_PHI(u_wind,v_wind,land_mask,mass_flux_PHI)
    integer(int64) :: i,j,col_num
    real(real64) :: f
    real(real64), dimension(:,:),allocatable :: tau_PHI,tau_THETA,lats,lats_data
    integer(int64), dimension(:,:),intent(in) :: land_mask
    real(real64), dimension(:,:),intent(in) :: u_wind, v_wind
    real(real64), dimension(:,:), intent(out) :: mass_flux_PHI

    allocate(tau_PHI(N_LATS,N_LONS))
    allocate(tau_THETA(N_LATS,N_LONS))
    allocate(lats(N_LATS,1))
    allocate(lats_data(N_LATS,2))

    col_num = 2 ! match int64 data type
    call calculate_surface_stress_PHI(u_wind,tau_PHI)
    call calculate_surface_stress_THETA(v_wind,tau_THETA)
    lats_data = read_file_real(LATS_FILE,N_LATS,col_num)
    lats(:,1)= lats_data(:,2) ! Second column in lats data file includes the latitude points

    do i = 1,N_LATS
        call calculate_coriolis_parameter(lats(i,1),f)
        do j = 1,N_LONS
            !! land mask if statement !!
            if (land_mask(i,j) == 1) then
                mass_flux_PHI(i,j) = 0
            else
                mass_flux_PHI(i,j) = (EPSILON*tau_PHI(i,j)-f*tau_THETA(i,j))/(EPSILON**2.0+f**2.0)
            end if
        end do
    end do

    deallocate(tau_PHI,tau_THETA,lats,lats_data)

end subroutine calculate_mass_flux_PHI
!!****************************************************************************
!! Functions to convert mass flux into the unit of Sverdrups for Theta and Phi
!!****************************************************************************
subroutine calculate_flow_sv_PHI(u_wind,v_wind,land_mask,sv_flow_PHI)
    integer(int64) :: i,j
    integer(int64), dimension(:,:),intent(in) :: land_mask
    real(real64), dimension(:,:), intent(in) :: u_wind,v_wind
    real(real64), dimension(:,:), intent(out) :: sv_flow_PHI
    real(real64),dimension(:,:),allocatable :: mass_flux_PHI

    allocate(mass_flux_PHI(N_LATS,N_LONS))

    call calculate_mass_flux_PHI(u_wind,v_wind,land_mask,mass_flux_PHI)

    do i = 1,N_LATS
        do j= 1,N_LONS
            sv_flow_PHI = mass_flux_PHI * R_PLANET * deg_to_rad(DELTA_LAT)  / RHO_WATER /1.0e6

        end do
    end do

    deallocate(mass_flux_PHI)

end subroutine calculate_flow_sv_PHI

subroutine calculate_flow_sv_THETA(u_wind,v_wind,land_mask,sv_flow_THETA)
    integer(int64) :: i,j,col_num
    real(real64) :: theta
    integer(int64), dimension(:,:),intent(in) :: land_mask
    real(real64), dimension(:,:), intent(in) :: u_wind,v_wind
    real(real64), dimension(:,:), allocatable :: lats_data,lats,mass_flux_THETA
    real(real64), dimension(:,:),intent(out) :: sv_flow_THETA

    allocate(lats(N_LATS,1))
    allocate(lats_data(N_LATS,2))
    allocate(mass_flux_THETA(N_LATS,N_LONS))

    col_num = 2
    lats_data = read_file_real(LATS_FILE,N_LATS,col_num)
    lats(:,1)= lats_data(:,2) ! Second column in lats data file includes the latitude points

    call calculate_mass_flux_THETA(u_wind,v_wind,land_mask,mass_flux_THETA)

    do i = 1,N_LATS
        theta = deg_to_rad(lats(i,1))
        do j= 1,N_LONS
            sv_flow_THETA(i,j) = mass_flux_THETA(i,j) * R_PLANET * cos(theta) * deg_to_rad(DELTA_LON) / RHO_WATER /1.0e6

        end do
    end do

    deallocate(lats,lats_data,mass_flux_THETA)

end subroutine calculate_flow_sv_THETA



end module MASS_FLUX
