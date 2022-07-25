!! This module has not been tested yet !!
!***************************************************************
! This module contains all the radiation flux functions for ULMO
!***************************************************************
module RAD_FLUXES
use NAMELIST
use WRITE_READ_DATA
use DEGREE_TO_RADIAN
implicit none
public :: calc_Q_flux
private
contains
!
!
!
function calc_Q_flux(T_surf) result(upward_Q_flux)
    real*8, dimension(N_LATS,N_LONS) :: F_net_sw_down
    real*8, dimension(N_LATS,N_LONS) :: F_lw_down
    real*8, dimension(N_LATS,N_LONS) :: F_latent_up
    real*8, dimension(N_LATS,N_LONS) :: F_sensible_up
    real*8, dimension(N_LATS,N_LONS) :: upward_Q_flux
    real*8, dimension(N_LATS,N_LONS), intent(in) :: T_surf
    integer :: i,j

    F_net_sw_down = read_file(SW_FLUX_NET_DOWN_DATA,N_LATS,N_LONS)
    F_lw_down     = read_file(LW_FLUX_DOWN_DATA,N_LATS,N_LONS)
    F_latent_up   = read_file(LATENT_UP_FLUX_DATA,N_LATS,N_LONS)
    F_sensible_up = read_file(SENSIBLE_UP_FLUX_DATA,N_LATS,N_LONS)

    do i = 1,N_LATS
        do j = 1, N_LONS
            upward_Q_flux(i,j) = F_latent_up(i,j)+F_sensible_up(i,j) &
            + EMISSIVITY*SIGMA*T_surf(i,j)**4 - (F_net_sw_down(i,j)+F_lw_down(i,j))
        end do
    end do




end function calc_Q_flux
!************************************************************************************************
! Function to read in atmospheric fluxes, net stellar sw flux (subrtracts reflected sw radiation
! from the incident)and incoming lw flux from atmospheric emission
!************************************************************************************************
function calculate_F_a() result(F_a)
    real*8, dimension(N_LATS,N_LONS) :: F_a
    real*8, dimension(N_LATS,N_LONS) :: F_net_sw_down
    real*8, dimension(N_LATS,N_LONS) :: F_lw_down
    real*8, dimension(N_LATS,N_LONS) :: F_latent_up
    real*8, dimension(N_LATS,N_LONS) :: F_sensible_up
    integer :: i,j

    F_net_sw_down = read_file(SW_FLUX_NET_DOWN_DATA,N_LATS,N_LONS)
    F_lw_down = read_file(LW_FLUX_DOWN_DATA,N_LATS,N_LONS)
    F_latent_up = read_file(LATENT_UP_FLUX_DATA,N_LATS,N_LONS)
    F_sensible_up = read_file(SENSIBLE_UP_FLUX_DATA,N_LATS,N_LONS)

    do i = 1,N_LATS
        do j=1,N_LONS
            F_a(i,j) = (F_net_sw_down(i,j)+F_lw_down(i,j))-(F_latent_up(i,j)+F_sensible_up(i,j))
        end do
    end do
end function
!**************************************************************************
! Checks if convetion condition is met, then calculates the associated flux
!***************************************************************************
function calculate_F_c(T) result(F_c)
    integer :: i,j,DEEP,SURFACE
    real*8, dimension(N_LATS,N_LONS):: F_c
    real*8, dimension(2,N_LATS,N_LONS), intent(in) :: T
    DEEP = 0
    SURFACE = 1
    do i=1,N_LATS-1 ! -1 as j<N_LATS
        do j=1,N_LONS-1 ! -1 as j<N_LATS
            if (T(DEEP,i,j)>T(SURFACE,i,j)) then ! convection condition
                F_c(i,j) = HEAT_TRANSFER_COEFFICIENT*(T(DEEP,i,j)-T(SURFACE,i,j))
            else
                F_c(i,j) = 0.0
            end if
        end do
    end do

end function
end module RAD_FLUXES
