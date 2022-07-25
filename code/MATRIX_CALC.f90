module MATRIX_CALC
use NAMELIST
use DEGREE_TO_RADIAN
use HEIGHT_OF_SLAB
use DIV_M
use dA_da
implicit none
public :: calculate_matrix_index,calculate_new_lon
private
contains
!********************************************************************************
! Calculates the index for the matrix used in the solver, before being compressed
!********************************************************************************
function calculate_matrix_index(lon,lat,height) result(ans)
    integer, intent(in) :: lon,lat,height
    integer :: ans

    ans = height*N_LATS*N_LONS+lat*N_LONS+lon
end function
!*****************************************************************************************************************************
! At maximum or minimum latitutde, the adjacent latitude point above or bellow, respectively, is located at the same latitude,
! but at the new longitudinal point (from lon to new_lon)
!**************************************************************************************
function calculate_new_lon(lon) result(new_lon)
    integer, intent(in) :: lon
    integer :: new_lon

    new_lon = lon+N_LONS/2
    if (new_lon>=N_LONS) then
        new_lon = new_lon - N_LONS
    end if
end function
!**************************************************************************************************************************************
! Calculates and sets the matrix elements for a 2 layer ocean with 144 longitude points and 90 latitude points.
! version determines whether the matrix is calculated for a no horizontal transport model, diffusion only model or the full Ekman model.
!**************************************************************************************************************************************
!function calculate_matrix() result()
!    real*8 ::d_theta, d_phi, s_dfsn, s_ekman,thickness,theta,Aij
!    real*8, dimension(N_LATS,N_LONS) :: dM_theta_dtheta,dM_phi_dphi,M_theta,M_phi,div_M
!    integer:: height
!
!
!    d_theta = deg_to_rad(DELTA_LAT)
!    d_phi   = deg_to_rad(DELTA_LON)
!    s_dfsn  = DELTA_T*D/(R_PLANET**2)
!    do height=SURFACE,height<N_DEPTHS
!        thickness = h_slab(height)
!        s_ekman = DELTA_T/(RHO_WATER*thickness*R_PLANET)
!        do lat =1, lat<N_LATS
!            theta = deg_to_rad()!need to read in latitude data
!        end do
!    end do
!
!
!end function
!**********************************************************************************************************************
! calculates and sets the vector b in matrix Ax=b. EMISSIVITY*SIGMA*T^4 is a blackbody emission from one of the layers,
! to ensure an energy exchange between the layers. 1.0/(RHO_WATER*C_V*H_S) prefactor converts a flux in W/m2 to K/s
!**********************************************************************************************************************

end module MATRIX_CALC
