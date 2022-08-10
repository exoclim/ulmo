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
function calculate_matrix() result(A)
  ! NO TRANSPORT !
  type(fgsl_spmatrix) :: A
  integer(fgsl_size_t) :: mat_index, i
  integer(fgsl_int) :: status
  real(fgsl_double) :: Aii
  integer(int64):: lat,lon,height

  !real(real64), dimension(N_LATS,1)      :: lats
  !real(real64), dimension(N_LATS,2)      :: lats_data_file

  !lats_data_file = read_file(LATS_FILE,N_LATS,col_num)
  !lats(:,1)= lats_data_file(:,2)

  ! matrix size = 144x90x2+90x144+144 = 39024 (height*N_LATS*N_LONS+lat*N_LONS+lon)
  A = fgsl_spmatrix_alloc(39024_fgsl_size_t, 39024_fgsl_size_t)


  do height = 0,N_DEPTHS-1
    do lat = 0,N_LATS-1
        do lon = 0,N_LONS-1


            mat_index = calculate_matrix_index(lat,lon,height)

            Aii = 1 ! no transport makes all diagonal components 1
            status = fgsl_spmatrix_set(A, mat_index, mat_index, Aii) ! build the sparse matrix

        end do
    end do
  end do





  ! writing matrix output to check !
!  write(output_unit, '(A)') 'printing all matrix elements:'
!  do i = 0,39024-1
!    write(*,*) 'A(',i,',',i,') = ',fgsl_spmatrix_get(A, i, i)
!  end do



end function calculate_matrix



!**********************************************************************************************************************
! calculates and sets the vector b in matrix Ax=b. EMISSIVITY*SIGMA*T^4 is a blackbody emission from one of the layers,
! to ensure an energy exchange between the layers. 1.0/(RHO_WATER*C_V*H_S) prefactor converts a flux in W/m2 to K/s
!**********************************************************************************************************************
function calculate_vector_b(T) result(v)
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

    do h = 1,N_DEPTHS
        do i = 1,N_LATS
            do j = 1,N_LONS
                depth=h_slab(h)
                if (h==SURFACE) then
                    b_hij = (DELTA_T/(RHO_WATER*C_V*depth))*(F_c(i,j)+F_a(i,j)-EMISSIVITY*SIGMA*(T(SURFACE,i,j))**4)*T(SURFACE,i,j)
                else
                    b_hij = (DELTA_T/(RHO_WATER*C_V*depth))*(-F_c(i,j)+T(h,i,j))
                end if
                vec_index = (h-1)*N_LATS*N_LONS+((i-1)*N_LONS+(j-1))
                v(vec_index) =  b_hij
                !print*, 'within function val =',v(vec_index)

            end do
        end do
    end do
    B = fgsl_vector_init(v)



    !testing!
!    size_vec =  fgsl_vector_get_size(B)
!    print*, 'size=', size_vec
!    print*, 'B(1)', v(:)
end function
!!**************************************************************
!! Subroutine that uses FGSL sparse linear algebra to solve Ax=b
!!**************************************************************
subroutine calculate_new_T(T) !result(T_new)
    type(fgsl_spmatrix) :: A,C
    real(real64), dimension(2,N_LATS,N_LONS), intent(in) :: T
    real(real64), dimension(2,N_LATS,N_LONS) :: T_new
    !integer(int64), dimension(N_LATS,N_LONS), intent(in) :: land_mask
    real(real64):: residual
    integer(fgsl_size_t),parameter:: n = 39024
    integer(int64) :: h,i,j,status,iu

    integer(fgsl_size_t) :: size_vec
    integer(fgsl_int) :: p

    !integer(fgsl_int) :: p
    integer(fgsl_size_t) :: iter = 0
    type(fgsl_vector) :: u,f
    type(fgsl_splinalg_itersolve) :: work
    real(fgsl_double), target :: u_f(n),f_f

    !equation: Au = f


    A = fgsl_spmatrix_alloc(n,n)

    work =  fgsl_splinalg_itersolve_alloc(fgsl_splinalg_itersolve_gmres,n,0_fgsl_size_t)


!!    ! initializing the values of b to 1 !
!!    v(1:p) = (/ (1._fgsl_double, p=1,n) /)
!!    b = fgsl_vector_init(v)
!!
!!    ! intializing the values of u to 0
    u_f(0:n-1) = (/(0._fgsl_double, p=0,n-1)/)
    u = fgsl_vector_init(u_f)
    size_vec =  fgsl_vector_get_size(u)
    print*, 'size of u =', size_vec
    print*,'u(0)= ', u_f(0)



    !f = calculate_vector_b(T)
    f_f = calculate_vector_b(T)
    print*,'pre init=',f_f(0)
    f = fgsl_vector_init(f_f)
    size_vec =  fgsl_vector_get_size(f)
    print*, 'size of f=', size_vec
    print*, 'f(0)= ', f_f(0)


    A = calculate_matrix()
    C = fgsl_spmatrix_compcol(A) ! compressed column format


!    !initial guess x = 0
!    u_f = 0
!
    do
        status = fgsl_splinalg_itersolve_iterate(C, f, TOL, u, work)

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
    do h = 0,N_DEPTHS-1 ! sets the new temperature from the fgsl vector
        do i = 0, N_LATS-1
            do j = 0,N_LONS-1
                T_new(h+1,i+1,j+1) =  u_f(h*N_LATS*N_LONS+(i*N_LONS)+j)
            end do
        end do
    end do

    print*,'T_new(1,1,1) = ',T_new(1,1,1)
    print*,'T_old(1,1,1)=', T(1,1,1)

    call fgsl_splinalg_itersolve_free(work)

    call fgsl_spmatrix_free(A)
    call fgsl_spmatrix_free(C)
    call fgsl_vector_free(f)
    call fgsl_vector_free(u)


end subroutine



end module mat_test
