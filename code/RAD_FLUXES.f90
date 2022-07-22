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

end module RAD_FLUXES
