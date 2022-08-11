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
public :: calculate_matrix , calculate_vector_b !, calculate_new_T
private

contains

!********************************************************************************
! Calculates the index for the matrix used in the solver, before being compressed
!********************************************************************************
function calculate_matrix_index(lon,lat,height) result(ans)
    integer(fgsl_size_t), intent(in) :: lon,lat,height
    integer(fgsl_size_t) :: ans

    ans = (height)*N_LATS*N_LONS+(lat)*N_LONS+lon
end function
!*****************************************************************************************************************************
! At maximum or minimum latitutde, the adjacent latitude point above or bellow, respectively, is located at the same latitude,
! but at the new longitudinal point (from lon to new_lon)
!**************************************************************************************
function calculate_new_lon(lon) result(new_lon)
    integer(fgsl_size_t), intent(in) :: lon
    integer(fgsl_size_t) :: new_lon

    new_lon = lon+N_LONS/2
    if (new_lon>=N_LONS) then
        new_lon = new_lon - N_LONS
    end if
end function calculate_new_lon
!**************************************************************************************************************************************
! Calculates and sets the matrix elements for a 2 layer ocean with 144 longitude points and 90 latitude points.
! version determines whether the matrix is calculated for a no horizontal transport model, diffusion only model or the full Ekman model.
!**************************************************************************************************************************************
function calculate_matrix() result(A)
  ! NO TRANSPORT !
  type(fgsl_spmatrix) :: A
  integer(fgsl_size_t) ::  i,lat,lon,height
  integer(fgsl_int) :: status
  real(fgsl_double) :: Aii
  integer, dimension(25920) :: mat_indexes

  ! matrix size = 144x90x2+90x144+144 = 39024 (height*N_LATS*N_LONS+lat*N_LONS+lon)
  A = fgsl_spmatrix_alloc(25920_fgsl_size_t, 25920_fgsl_size_t)


  do height = 0,N_DEPTHS-1
    do lat = 0,N_LATS-1
        do lon = 0,N_LONS-1


            i = calculate_matrix_index(lon,lat,height)
            !print*,i

            !mat_indexes(i+1) = i

            Aii = 1 ! no transport makes all diagonal components 1
            status = fgsl_spmatrix_set(A, i, i, Aii) ! build the sparse matrix
            !print*,'A matrix','(',i,2,')','=' ,fgsl_spmatrix_get(A,i,2_fgsl_size_t)

        end do
    end do
  end do
!print*,'mat index = ',mat_indexes
!print*,'mat index = ',size(mat_indexes)
end function calculate_matrix



!**********************************************************************************************************************
! calculates and sets the vector b in matrix Ax=b. EMISSIVITY*SIGMA*T^4 is a blackbody emission from one of the layers,
! to ensure an energy exchange between the layers. 1.0/(RHO_WATER*C_V*H_S) prefactor converts a flux in W/m2 to K/s
!**********************************************************************************************************************
subroutine calculate_vector_b(T,v,B)
    real(fgsl_double),dimension(2,N_LATS,N_LONS),intent(in) :: T
    !integer(int64), dimension(N_LATS,N_LONS), intent(in) :: land_mask
    real(real64), dimension(N_LATS,N_LONS) :: F_a,F_c
    real(real64), dimension(25920) :: vec_data
    real(real64) :: depth
    integer(int64) :: h,i ,j,SURFACE,DEEP,iu
    !real(real64), parameter:: tol = TOL
    !real(fgsl_size_t),parameter:: max_iter = MAX_ITER
    ! fgsl !
    real(real64) :: b_hij
    integer(fgsl_size_t), parameter :: n = 25920
    real(fgsl_double), target,intent(out) :: v(n)
    integer(fgsl_size_t) :: vec_index,size_vec
    real(fgsl_double) :: vec_val
    type(fgsl_vector),intent(out) :: B

    SURFACE = 1 ! position in 3D temperatures array (1,:,:)
    DEEP = 2    ! position in 3D temperatures array (2,:,:)

    !b_hij = 0.0

    F_a = calculate_F_a()
    F_c = calculate_F_c(T)

    !testing (works)!
!    h = 2
!    depth=h_slab(h)
!    print*,'depth deep',depth
!    b_hij = (DELTA_T/(RHO_WATER*C_V*depth))*(-F_c(1,1)+T(2,1,1))
!    print*,'deep = ',b_hij
!    h=1
!    depth= h_slab(h)
!    print*,'depth deep',depth
!    b_hij = (DELTA_T/(RHO_WATER*C_V*depth))*(F_c(1,1)+F_a(1,1)-EMISSIVITY*SIGMA*(T(SURFACE,1,1))**4)*T(SURFACE,1,1)
!    print*,'surf=',b_hij
    do h = 1,N_DEPTHS
        do i = 1, N_LATS
            do j= 1,N_LONS

                depth = h_slab(h)
                if (h==1) then
                    b_hij = (DELTA_T/(RHO_WATER*C_V*depth))*(F_c(i,j)+F_a(i,j)-EMISSIVITY*SIGMA*(T(h,i,j))**4)*T(h,i,j)
                else
                    b_hij = (DELTA_T/(RHO_WATER*C_V*depth))*(-F_c(i,j)+T(h,i,j))
                end if
                !print*,b_hij

                vec_index = (h-1)*N_LATS*N_LONS+(i-1)*N_LONS+(j-1) ! fgsl index
                vec_data(vec_index+1) = b_hij ! needs fortran index
                !print*, vec_index
                vec_val = b_hij ! making sure the type is fgsl_size_t
                !print*, vec_val
                v(vec_index+1) = vec_val



            end do
        end do
    end do


print*, 'vec = ',v
B = fgsl_vector_init(v)
size_vec =  fgsl_vector_get_size(B)
!print*, 'size=', size_vec


end subroutine
!!**************************************************************
!! Subroutine that uses FGSL sparse linear algebra to solve Ax=b
!!**************************************************************
!subroutine calculate_new_T(T,T_new) !result(T_new)
!    type(fgsl_spmatrix) :: A,C
!    real(fgsl_double), dimension(2,N_LATS,N_LONS), intent(in) :: T
!    real(fgsl_double), dimension(2,N_LATS,N_LONS),intent(out) :: T_new
!    !integer(int64), dimension(N_LATS,N_LONS), intent(in) :: land_mask
!    real(fgsl_double):: residual
!    integer(fgsl_size_t),parameter:: n = 25920
!    integer(int64) :: h,i,j,iu
!
!    integer(fgsl_size_t) :: size_vec
!    integer(fgsl_int) :: p,status
!    integer(fgsl_size_t) :: iter = 0
!    type(fgsl_vector) :: u,f
!    type(fgsl_splinalg_itersolve) :: work
!    real(fgsl_double), target :: u_f(n),f_f(n)
!
! !   equation: Au = f
!
!
!    A = fgsl_spmatrix_alloc(n,n)
!
!
!    work =  fgsl_splinalg_itersolve_alloc(fgsl_splinalg_itersolve_gmres,n,0_fgsl_size_t)
!
!
!!    ! initializing the values of b to 1 !
!!    v(1:p) = (/ (1._fgsl_double, p=1,n) /)
!!    b = fgsl_vector_init(v)
!!
!!    ! intializing the values of u to 0
!    u_f(0:n-1) = (/(0._fgsl_double, p=0,n-1)/)
!    u = fgsl_vector_init(u_f)
!    size_vec =  fgsl_vector_get_size(u)
!    !print*, 'size of u =', size_vec
!    !print*,'vector u= ', u_f
!
!
!
!
!    call calculate_vector_b(T,f_f,f)
!    !print*,'f vector=',f_f
!    size_vec =  fgsl_vector_get_size(f)
!    !print*, 'size of f=', size_vec
!
!
!    A = calculate_matrix()
!    C = fgsl_spmatrix_compcol(A) ! compressed column format
!!
!!
!!    !initial guess x = 0
!!    u_f = 0
!!
!    do
!        status = fgsl_splinalg_itersolve_iterate(C, f, TOL, u, work)
!
!         !print out residual norm ||A*x-b||
!        residual = fgsl_splinalg_itersolve_normr(work)
!        write(output_unit, '(A,I2,A,G15.6)') 'iter ', iter, ' residual = ', residual
!
!        if (status == FGSL_SUCCESS) then
!            write(output_unit, '(A)') 'Converged'
!        end if
!        iter = iter + 1
!        if (status /= FGSL_CONTINUE .or. iter >= MAX_ITER) exit
!    end do
!
!!     output solution !
!
!    do h = 0_fgsl_size_t,N_DEPTHS-1_fgsl_size_t ! sets the new temperature from the fgsl vector
!        do i = 0_fgsl_size_t, N_LATS-1_fgsl_size_t
!            do j = 0_fgsl_size_t,N_LONS-1_fgsl_size_t
!                T_new(h+1,i+1,j+1) =  u_f(h*N_LATS*N_LONS+(i*N_LONS)+j)
!            end do
!        end do
!    end do
!
!    print*,'T_new(1,1,1) = ',T_new(1,1,1)
!    print*,'T_old(1,1,1)=', T(1,1,1)
!
!    call fgsl_splinalg_itersolve_free(work)
!
!    call fgsl_spmatrix_free(A)
!    call fgsl_spmatrix_free(C)
!    call fgsl_vector_free(f)
!    call fgsl_vector_free(u)
!
!
!end subroutine



end module mat_test
