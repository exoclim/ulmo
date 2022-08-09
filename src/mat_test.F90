!*********************************************
! Module testing fgsl ready for implementation
!*********************************************
module mat_test
use NAMELIST
use DEGREE_TO_RADIAN
use READ_DATA
use Other_FLUXES
use HEIGHT_OF_SLAB
use fgsl
use, intrinsic :: iso_fortran_env
implicit none
public :: calculate_matrix , calculate_vector_b , calculate_new_T
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
function calculate_matrix(lat,lon,height) result(A)
  ! NO TRANSPORT !
  type(fgsl_spmatrix) :: A
  integer(fgsl_size_t) :: mat_index
  integer(fgsl_int) :: status
  real(fgsl_double) :: Aij, Aii
  integer(int64), intent(in):: lat,lon,height

  ! matrix size = 144x90x2+90x144+144 = 39024 (height*N_LATS*N_LONS+lat*N_LONS+lon)
  A = fgsl_spmatrix_alloc(39024_fgsl_size_t, 39024_fgsl_size_t)
  mat_index = calculate_matrix_index(lat,lon,height)
  Aii = 1 ! no transport makes all diagonal components 1

  ! build the sparse matrix
  status = fgsl_spmatrix_set(A, mat_index, mat_index, Aii)

  ! writing matrix output to check !
  !write(output_unit, '(A)') 'printing all matrix elements:'
  !write(*,*) 'A(',1,',',1,') = ',fgsl_spmatrix_get(A, 1_fgsl_size_t, 1_fgsl_size_t)

end function calculate_matrix



!**********************************************************************************************************************
! calculates and sets the vector b in matrix Ax=b. EMISSIVITY*SIGMA*T^4 is a blackbody emission from one of the layers,
! to ensure an energy exchange between the layers. 1.0/(RHO_WATER*C_V*H_S) prefactor converts a flux in W/m2 to K/s
!**********************************************************************************************************************
function calculate_vector_b(T) result(B)
    real(real64),dimension(2,N_LATS,N_LONS),intent(in) :: T
    !integer(int64), dimension(N_LATS,N_LONS), intent(in) :: land_mask
    real(real64), dimension(N_LATS,N_LONS) :: F_a,F_c
    real(real64) :: depth
    integer(int64) :: h, N_DEPTHS,i ,j,SURFACE,DEEP,iu
    !real(real64), parameter:: tol = TOL
    !real(fgsl_size_t),parameter:: max_iter = MAX_ITER
    ! fgsl !
    real(real64) :: b_hij
    integer(fgsl_size_t), parameter :: ndim = 39024
    real(fgsl_double), target :: v(ndim)
    integer(fgsl_size_t) :: vec_index, size_vec
    type(fgsl_vector) :: B

    SURFACE = 1 ! position in 3D temperatures array (1,:,:)
    DEEP = 2    ! position in 3D temperatures array (2,:,:)

    !b_hij = 0.0

    F_a = calculate_F_a()
    F_c = calculate_F_c(T)

    !testing (works)!
!    h = 2
!    depth=h_slab(h)
!    b_hij = (DELTA_T/(RHO_WATER*C_V*depth))*(-F_c(90,144)+T(2,90,144))
!    print*,'deep = ',b_hij
!     h=1
!     depth= h_slab(h)
!     b_hij = (DELTA_T/(RHO_WATER*C_V*depth))*(F_c(90,144)+F_a(90,144)-EMISSIVITY*SIGMA*(T(SURFACE,90,144))**4)*T(SURFACE,90,144)
!     print*,'surf=',b_hij

    do h = 1,2 ! N_DEPTHS = 2
        do i = 1,N_LATS
            do j = 1,N_LONS
                depth=h_slab(h)
                if (h==SURFACE) then
                    b_hij = (DELTA_T/(RHO_WATER*C_V*depth))*(F_c(i,j)+F_a(i,j)-EMISSIVITY*SIGMA*(T(SURFACE,i,j))**4)*T(SURFACE,i,j)
                else
                    b_hij = (DELTA_T/(RHO_WATER*C_V*depth))*(-F_c(i,j)+T(h,i,j))
                end if
                vec_index = h*N_LATS*N_LONS+(i*N_LONS+j)
                v(vec_index) =  b_hij

            end do
        end do
    end do
    B = fgsl_vector_init(v)


    !testing!
!    size_vec =  fgsl_vector_get_size(B)
!    print*, 'size=', size_vec
!    print*, 'B(1)', v(:)
end function
!**************************************************************
! Subroutine that uses FGSL sparse linear algebra to solve Ax=b
!**************************************************************
subroutine calculate_new_T(T,lat,lon,height) !result(T_new)
    type(fgsl_spmatrix) :: A,C
    integer(int64), intent(in):: lat,lon,height
    real(real64), dimension(2,N_LATS,N_LONS), intent(in) :: T
    real(real64), dimension(2,N_LATS,N_LONS) :: T_new
    !integer(int64), dimension(N_LATS,N_LONS), intent(in) :: land_mask
    real(real64):: residual
    integer(fgsl_size_t),parameter:: n = N_DEPTHS*N_LATS*N_LONS
    integer(int64) :: h,i,j,status
    real(fgsl_double), target :: v(n),q(n),u_x(n)
    integer(fgsl_int) :: p
    real(fgsl_size_t) :: iter
    type(fgsl_vector) :: b,x
    type(fgsl_splinalg_itersolve) :: work

    A = fgsl_spmatrix_alloc(n,n)


!    ! initializing the values of b to 1 !
!    v(1:p) = (/ (1._fgsl_double, p=1,n) /)
!    b = fgsl_vector_init(v)
!
!    ! intializing the values of x to 0
    q(1:p) = (/(0._fgsl_double, p=1,n)/)
    x = fgsl_vector_init(q)

    A = calculate_matrix(lat,lon,height)
    b = calculate_vector_b(T)
    C = fgsl_spmatrix_compcol(A) ! compressed column format
    work =  fgsl_splinalg_itersolve_alloc(fgsl_splinalg_itersolve_gmres,n,0_fgsl_size_t)
    iter = 0
    !initial guess x = 0
    u_x = 0

    do
        status = fgsl_splinalg_itersolve_iterate(C, b, TOL, x, work)

        ! print out residual norm ||A*x-b||
        residual = fgsl_splinalg_itersolve_normr(work)
        write(output_unit, '(A,I2,A,G15.6)') 'iter ', iter, ' residual = ', residual

        if (status == FGSL_SUCCESS) then
            write(output_unit, '(A)') 'Converged'
        end if
        iter = iter + 1
        if (status /= FGSL_CONTINUE .or. iter >= MAX_ITER) exit
    end do

    ! output solution !
! issue here !
    do h = 1,N_DEPTHS ! sets the new temperature from the fgsl vector
        do i = 1, N_LATS
            do j = 1,N_LONS
                T_new(h,i,j) = q(h*N_LATS*N_LONS+(i*N_LONS+j))
            end do
        end do
    end do


end subroutine



end module mat_test
