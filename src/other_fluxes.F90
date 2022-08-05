!! This module has not been tested yet !!
!***************************************************************
! This module contains all the radiation flux functions for ULMO
!***************************************************************
module Other_FLUXES
use NAMELIST
use READ_DATA
use DEGREE_TO_RADIAN
implicit none
public :: calc_Q_flux, calculate_F_c, calculate_F_a
private
contains
!
!
!
function calc_Q_flux(T) result(upward_Q_flux)
    real, dimension(N_LATS,N_LONS) :: F_net_sw_down, F_lw_down,F_latent_up,F_sensible_up,upward_Q_flux
    real, dimension(2,N_LATS,N_LONS), intent(in) :: T
    integer :: i,j


    F_net_sw_down = read_file(SW_FLUX_NET_DOWN_DATA,N_LATS,N_LONS)
    F_lw_down     = read_file(LW_FLUX_DOWN_DATA,N_LATS,N_LONS)
    F_latent_up   = read_file(LATENT_UP_FLUX_DATA,N_LATS,N_LONS)
    F_sensible_up = read_file(SENSIBLE_UP_FLUX_DATA,N_LATS,N_LONS)

    do i = 1,N_LATS
        do j = 1, N_LONS
            upward_Q_flux(i,j) = F_latent_up(i,j)+F_sensible_up(i,j) &
            + EMISSIVITY*SIGMA*T(1,i,j)**4 - (F_net_sw_down(i,j)+F_lw_down(i,j))
        end do
    end do




end function calc_Q_flux
!************************************************************************************************
! Function to read in atmospheric fluxes, net stellar sw flux (subrtracts reflected sw radiation
! from the incident)and incoming lw flux from atmospheric emission
!************************************************************************************************
function calculate_F_a() result(F_a)
    real, dimension(N_LATS,N_LONS) :: F_a, F_net_sw_down,F_lw_down,F_latent_up,F_sensible_up
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
! Checks if convection condition is met, then calculates the associated flux
!***************************************************************************
function calculate_F_c(T) result(F_c)
    integer :: i,j
    real, dimension(N_LATS,N_LONS):: F_c
    real, dimension(2,N_LATS,N_LONS), intent(in) :: T
    ! h = 1 surface
    ! h = 2 deep
    do i=1,N_LATS
        do j=1,N_LONS
            if (T(2,i,j)>T(1,i,j)) then ! convection condition
                F_c(i,j) = HEAT_TRANSFER_COEFFICIENT*(T(2,i,j)-T(1,i,j))
            else
                F_c(i,j) = 0.0
            end if
        end do
    end do

end function
end module Other_FLUXES
