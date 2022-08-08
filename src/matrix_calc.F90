!******************************************
! Module for Matrix calculations using fgsl
!******************************************
module MATRIX_CALC
use NAMELIST
use DEGREE_TO_RADIAN
use READ_DATA
use :: fgsl
use, intrinsic :: iso_fortran_env
implicit none
public :: calculate_matrix
private

!! This works and has been tested !!

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
subroutine calculate_matrix(lat,lon,height)
  ! NO TRANSPORT !
  type(fgsl_spmatrix) :: A
  integer(fgsl_size_t) :: mat_index
  integer(fgsl_int) :: status
  real(fgsl_double) :: Aii
  integer(int64), intent(in):: lat,lon,height

  ! matrix size = 144x90x2+90x144+144 = 39024 (height*N_LATS*N_LONS+lat*N_LONS+lon)
  A = fgsl_spmatrix_alloc(39024_fgsl_size_t, 39024_fgsl_size_t)
  mat_index = calculate_matrix_index(lat,lon,height)
  Aii = 1 ! no transport makes all diagonal components 1

  ! build the sparse matrix
  status = fgsl_spmatrix_set(A, mat_index, mat_index, Aii)

  ! writing matrix output to check !
  write(output_unit, '(A)') 'printing all matrix elements:'
  write(*,*) 'A(',mat_index,',',mat_index,') = ',fgsl_spmatrix_get(A, mat_index, mat_index)

end subroutine calculate_matrix



end module MATRIX_CALC
