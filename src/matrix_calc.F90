!******************************************
! Module for Matrix calculations using fgsl
!******************************************
module MATRIX_CALC
use Constants
use DEGREE_TO_RADIAN, only: deg_to_rad
use heat_fluxes, only: calculate_F_c,calc_Q_flux
use HEIGHT_OF_SLAB, only: h_slab
use fgsl
use, intrinsic :: iso_fortran_env
implicit none
public :: calculate_matrix , calculate_vector_f_values,calculate_matrix_index !, calculate_new_T
private

contains

!********************************************************************************
! Calculates the index for the matrix used in the solver, before being compressed
!********************************************************************************
function calculate_matrix_index(lon,lat,height) result(ans)
    integer(int64), intent(in) :: lon,lat,height
    integer(int64) :: ans

    ans = (height)*N_LATS*N_LONS+(lat)*N_LONS+lon
end function
!*****************************************************************************************************************************
! At maximum or minimum latitude, the adjacent latitude point above or bellow, respectively, is located at the same latitude,
! but at the new longitudinal point (from lon to new_lon)
!******************************************************************************************************************************
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
subroutine calculate_matrix(land_mask,A)
  ! NO TRANSPORT !
  integer(int64),dimension(:,:),intent(in) :: land_mask
  type(fgsl_spmatrix),intent(out) :: A
  integer(fgsl_size_t) ::  i
  integer(int64):: lat,lon,height
  integer(fgsl_int) :: status
  real(fgsl_double) :: Aii
  !integer, dimension(25920) :: mat_indexes

  ! matrix size = 144x90x2+90x144+144 = 39024 (height*N_LATS*N_LONS+lat*N_LONS+lon)
  A = fgsl_spmatrix_alloc(25920_fgsl_size_t, 25920_fgsl_size_t)


  do height = 0,N_DEPTHS-1
    do lat = 0,N_LATS-1
        do lon = 0,N_LONS-1

            i = calculate_matrix_index(lon,lat,height)
            !print*,i
            if(land_mask(lat+1,lon+1) == 1) then
                Aii = 1
                status = fgsl_spmatrix_set(A,i,i,Aii)
                cycle
            end if
            !mat_indexes(i+1) = i

            Aii = 1 ! no transport makes all diagonal components 1
            status = fgsl_spmatrix_set(A, i, i, Aii) ! build the sparse matrix
            !print*,'A matrix','(',i,2,')','=' ,fgsl_spmatrix_get(A,i,2_fgsl_size_t)

        end do
    end do
  end do
!print*,'mat index = ',mat_indexes
!print*,'mat index = ',size(mat_indexes)
end subroutine calculate_matrix



!**********************************************************************************************************************
! calculates and sets the vector b in matrix Ax=b. EMISSIVITY*SIGMA*T^4 is a blackbody emission from one of the layers,
! to ensure an energy exchange between the layers. 1.0/(RHO_WATER*C_V*H_S) prefactor converts a flux in W/m2 to K/s
!**********************************************************************************************************************
subroutine calculate_vector_f_values(land_mask,T,f_f,F_a,F_c)!,B)
    real(fgsl_double),dimension(:,:,:),intent(in) :: T
    integer(int64), dimension(:,:), intent(in) :: land_mask
    real(real64), dimension(:,:),intent(in) :: F_a,F_c
    !real(real64), dimension(:),allocatable :: vec_data
    real(real64) :: depth,b_hij
    integer(int64) :: h,i ,j,SURFACE,DEEP
    ! fgsl !
    !integer(fgsl_size_t), parameter :: n = 25920
    real(fgsl_double),dimension(:),intent(out) :: f_f
    integer(fgsl_size_t) :: vec_index !,size_vec
    real(fgsl_double) :: vec_val

    do h = 1,N_DEPTHS
        do i = 1, N_LATS
            do j= 1,N_LONS

                if(land_mask(i,j) == 1) then
                    b_hij = T(h,i,j)
                else
                    depth = h_slab(h)
                    if (h==1) then
                        b_hij = DELTA_T/(RHO_WATER*C_V*depth)*(F_c(i,j)+F_a(i,j)-EMISSIVITY*SIGMA*((T(h,i,j))**4))+T(h,i,j)
                    else
                        b_hij = DELTA_T/(RHO_WATER*C_V*depth)*(-F_c(i,j))+T(h,i,j) !! double check this with Jake !!
                    end if

                end if
                vec_index = (h-1)*N_LATS*N_LONS+(i-1)*N_LONS+(j-1) ! fgsl index
                vec_val = b_hij ! making sure the type is fgsl_size_t
                f_f(vec_index+1) = vec_val

            end do
        end do
    end do



end subroutine calculate_vector_f_values

end module MATRIX_CALC




