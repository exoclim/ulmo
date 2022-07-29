module TRANSPORT
use NAMELIST
use DEGREE_TO_RADIAN
use RAD_FLUXES
implicit none
public :: no_transport_no_deep
private

contains
!***************************************
! No heat transport, No deep ocean layer
!***************************************
function no_transport_no_deep(T_n) result(T_n1) ! result is the temperature at next time step T_n+1
    real*8, dimension(N_LATS,N_LONS) :: F_a,T_n1 ! start with ignoring the deep i.e F_c = 0
    real*8, dimension(N_LATS,N_LONS),intent(in) :: T_n
    integer :: i,j

    F_a = calculate_F_a()

    do i = 1,N_LATS
        do j = 1,N_LONS
            T_n1(i,j) = ( DELTA_T/(RHO_WATER*C_V*H_S) )*( F_a(i,j)-SIGMA*EPSILON*(T_n(i,j))**4.0 ) + T_n(i,j)

        end do
    end do
end function no_transport_no_deep
!***********************************************
! Diffusion no deep ocean layer
!**********************************************
function diffusion_no_deep(T_n) result(T_n1)
    real*8, dimension(N_LATS,N_LONS) :: F_a,T_nl
    real*8, dimension(N_LATS,N_LONS), intent(in) :: T_n
    integer :: i,j

    F_a = calculate_F_a()

    do i = 1,N_LATS
        do j = 1,N_LONS
            T_nl(i,j) = (DELTA_T/(RHO_WATER*C_V*H_S))*( F_a(i,j)+SIGMA*EPSILON*(T_n(i,j))**4.0 ) + D*(1/R_PLANET)*DELTA_T*
        end do
    end do








end function diffusion_no_deep

end module TRANSPORT
