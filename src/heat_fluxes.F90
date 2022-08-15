!! This module has not been tested yet !!
!***************************************************************
! This module contains all the radiation flux functions for ULMO
!***************************************************************
module heat_fluxes
use NAMELIST
use READ_DATA
use DEGREE_TO_RADIAN
use, intrinsic :: iso_fortran_env
implicit none
public :: calc_Q_flux, calculate_F_c, calculate_F_a
private
contains



subroutine calc_Q_flux(T,upward_Q_flux)
    real(real64), dimension(:,:), allocatable :: F_net_sw_down, F_lw_down,F_latent_up,F_sensible_up
    integer(int64) :: i,j
    real(real64), dimension(:,:,:), intent(in) :: T
    real(real64), dimension(:,:), intent(out) :: upward_Q_flux

    allocate(F_net_sw_down(N_LATS,N_LONS),F_lw_down(N_LATS,N_LONS),F_latent_up(N_LATS,N_LONS),F_sensible_up(N_LATS,N_LONS))


    F_net_sw_down = read_file_real(SW_FLUX_NET_DOWN_DATA,N_LATS,N_LONS)
    F_lw_down     = read_file_real(LW_FLUX_DOWN_DATA,N_LATS,N_LONS)
    F_latent_up   = read_file_real(LATENT_UP_FLUX_DATA,N_LATS,N_LONS)
    F_sensible_up = read_file_real(SENSIBLE_UP_FLUX_DATA,N_LATS,N_LONS)

    do i = 1,N_LATS
        do j = 1, N_LONS
            upward_Q_flux(i,j) = F_latent_up(i,j)+F_sensible_up(i,j) &
            + EMISSIVITY*SIGMA*T(1,i,j)**4 - (F_net_sw_down(i,j)+F_lw_down(i,j))
        end do
    end do

    deallocate(F_net_sw_down,F_lw_down,F_latent_up,F_sensible_up)


end subroutine calc_Q_flux
!************************************************************************************************
! Function to read in atmospheric fluxes, net stellar sw flux (subrtracts reflected sw radiation
! from the incident)and incoming lw flux from atmospheric emission
!************************************************************************************************
subroutine calculate_F_a(F_a)
    real(real64), dimension(:,:), allocatable :: F_net_sw_down,F_lw_down,F_latent_up,F_sensible_up
    integer(int64) :: i,j
    real(real64), dimension(N_LATS,N_LONS),intent(out) :: F_a

    allocate(F_net_sw_down(N_LATS,N_LONS),F_lw_down(N_LATS,N_LONS),F_latent_up(N_LATS,N_LONS),F_sensible_up(N_LATS,N_LONS))

    F_net_sw_down = read_file_real(SW_FLUX_NET_DOWN_DATA,N_LATS,N_LONS)
    F_lw_down = read_file_real(LW_FLUX_DOWN_DATA,N_LATS,N_LONS)
    F_latent_up = read_file_real(LATENT_UP_FLUX_DATA,N_LATS,N_LONS)
    F_sensible_up = read_file_real(SENSIBLE_UP_FLUX_DATA,N_LATS,N_LONS)

    do i = 1,N_LATS
        do j=1,N_LONS
            F_a(i,j) = (F_net_sw_down(i,j)+F_lw_down(i,j))-(F_latent_up(i,j)+F_sensible_up(i,j))
        end do
    end do

    deallocate(F_net_sw_down,F_lw_down,F_latent_up,F_sensible_up)

end subroutine calculate_F_a
!**************************************************************************
! Checks if convection condition is met, then calculates the associated flux
!***************************************************************************
subroutine calculate_F_c(T,F_c)
    integer(int64) :: i,j
    real(real64), dimension(:,:),intent(out):: F_c
    real(real64), dimension(:,:,:), intent(in) :: T
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

end subroutine calculate_F_c
end module heat_fluxes
