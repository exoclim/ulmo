!*********************************************
! Module testing fgsl ready for implementation
!*********************************************
module mat_test
use NAMELIST
use DEGREE_TO_RADIAN
use, intrinsic :: iso_fortran_env
use fgsl
implicit none
private
contains

!********************************************************************************
! Calculates the index for the matrix used in the solver, before being compressed
!********************************************************************************
function calculate_matrix_index(lon,lat,height) result(ans)
    integer(int64), intent(in) :: lon,lat,height
    integer(int64) :: ans

    ans = height*N_LATS*N_LONS+lat*N_LONS+lon
end function
!*****************************************************************************************************************************
! At maximum or minimum latitutde, the adjacent latitude point above or bellow, respectively, is located at the same latitude,
! but at the new longitudinal point (from lon to new_lon)
!**************************************************************************************
function calculate_new_lon(lon) result(new_lon)
    integer(int64), intent(in) :: lon
    integer(int64) :: new_lon

    new_lon = lon+N_LONS/2
    if (new_lon>=N_LONS) then
        new_lon = new_lon - N_LONS
    end if
end function calculate_new_lon
!**************************************************************************************************************************************
! Calculates and sets the matrix elements for a 2 layer ocean with 144 longitude points and 90 latitude points.
! version determines whether the matrix is calculated for a no horizontal transport model, diffusion only model or the full Ekman model.
!**************************************************************************************************************************************
subroutine calculate_matrix()
    type(fgsl_spmatrix) :: A
    real(fgsl_double) ::Aij, Aii
    integer(int64):: lat,lon,height
    integer(fgsl_size_t) :: i
    integer(fgsl_int) :: status

    i = calculate_matrix_index(lat,lon,height)
    Aii = 1
    status = fgsl_spmatrix_set(A, i, i, Aii)


end subroutine


end module mat_test
