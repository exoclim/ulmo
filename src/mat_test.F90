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
public :: calculate_matrix, calculate_vector_b
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
subroutine calculate_matrix(lat,lon,height)
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
  write(output_unit, '(A)') 'printing all matrix elements:'
  write(*,*) 'A(',mat_index,',',mat_index,') = ',fgsl_spmatrix_get(A, mat_index, mat_index)

end subroutine calculate_matrix



!**********************************************************************************************************************
! calculates and sets the vector b in matrix Ax=b. EMISSIVITY*SIGMA*T^4 is a blackbody emission from one of the layers,
! to ensure an energy exchange between the layers. 1.0/(RHO_WATER*C_V*H_S) prefactor converts a flux in W/m2 to K/s
!**********************************************************************************************************************
subroutine calculate_vector_b(T)!,land_mask)
    real(real64),dimension(2,N_LATS,N_LONS),intent(in) :: T
    !integer(int64), dimension(N_LATS,N_LONS), intent(in) :: land_mask
    real(real64), dimension(N_LATS,N_LONS) :: F_a,F_c
    real(real64) :: depth
    integer(int64) :: h, N_DEPTHS,i ,j,SURFACE,DEEP
    ! fgsl !
    real(fgsl_double) :: b_hij
    integer(fgsl_size_t), parameter :: ndim = 39024
    real(fgsl_double), target :: v(ndim)
    integer(fgsl_size_t) :: vec_index, size_vec
    type(fgsl_vector) :: B

    SURFACE = 1 ! position in 3D temperatures array (1,:,:)
    DEEP = 2    ! position in 3D temperatures array (2,:,:)

    b_hij = 0.0

    F_a = calculate_F_a()
    F_c = calculate_F_c(T)



    do h = 1, N_DEPTHS
        do i = 1,N_LATS
            do j = 1,N_LONS
!                if (land_mask(i,j) == 1) then
!                    b_hij = T(h,i,j)
!                else
                    depth=h_slab(h)
                    if (h==SURFACE) then
                        b_hij = (DELTA_T/(RHO_WATER*C_V*depth))*(F_c(i,j)+F_a(i,j)-EMISSIVITY*SIGMA*(T(SURFACE,i,j))**4)*T(h,i,j)
                    else
                        b_hij = (DELTA_T/(RHO_WATER*C_V*depth))*(-F_c(i,j)+T(h,i,j))

                    end if

                !end if
                vec_index = calculate_matrix_index(i,j,h)
                v(vec_index) =  b_hij
            end do
        end do
    end do
B = fgsl_vector_init(v)
size_vec =  fgsl_vector_get_size(B)
write(*,*) 'element 1 =',v(1)


! This works !
!integer(fgsl_int) :: p
!integer(fgsl_size_t), parameter :: ndim = 12
!v(1:ndim) =  (/ (dble(p)+0.1_fgsl_double, p=1,ndim) /)
!B = fgsl_vector_init(v)
!write(*,*) v(3::3)
!call fgsl_vector_free(B)
!
!i.e v(1:ndim) = (/(p+0.1,p=1,ndim))

end subroutine
!
!subroutine calculate_new_T(T,land_mask)
!    integer(int64) :: n, N_DEPTHS
!    real(real64), dimension(N_LATS,N_LONS), intent(in) :: T, land_mask
!    real(real64):: tol, max_iter, residual
!    integer(int64) :: h,i,j
!!
!!    n = N_DEPTHS*n_lats*N_LONS
!!    fgsl_spmatrix *A = fgsl_spmatrix_alloc(n,n)
!!    fgsl_spmatrix *C
!!    fgsl_vector *b = fgsl_vector_alloc(n)
!!    fgsl_vector *x = fgsl_vector_alloc(n)
!!    call calculate_matrix()
!!    call calculate_vector_b(T,land_mask)
!!    C = fgsl_spmatrix_ccs(A)
!!
!!    ! now solve the system with the GMRES iterative solver
!!    tol = TOL ! solution relative tolerance
!!    max_iter = MAX_ITER ! maximum iterations
!!    gsl_splinalg_itersolve_type *S = gsl_splinalg_itersolve_gmres ! not sure
!!    gsl_splinalg_itersolve *work = gsl_splinalg_itersolve_alloc(S, n, 0)
!!    iter = 0
!!    ! initial guess x = 0
!!    gsl_vector_set_zero(x)
!!
!!    !solve the system A x = b
!!     do while (status == GSL_CONTINUE && ++iter < max_iter)
!!        status = gsl_splinalg_itersolve_iterate(C, b, tol, x, work)
!!
!!        if (status /= GSL_SUCCESS && iter==max_iter-1 ) then! solver exits if solution diverges
!!            residual = gsl_splinalg_itersolve_normr(work)
!!            ! print statements
!!            exit(SOLVER_ERR)
!!        end if
!!
!!     end do
!!    ! output solution
!!    do (h = 1,N_DEPTHS) ! sets the new temperature from gsl vector
!!        do (i = 1,N_LATS)
!!            do (j = 1,N_LONS)
!!                T(h,i,j) = gsl_vector_get(x,h*N_LATS*N_LONS+(i*N_LONS+j))
!!            end do
!!        end do
!!    end do
!!
!!    fgsl_splinalg_itersolve_free(work);
!!
!!	fgsl_spmatrix_free(A);
!!	fgsl_spmatrix_free(C);
!!	fgsl_vector_free(b);
!!	fgsl_vector_free(x);
!
!end subroutine
!
!subroutine time_stepper(n_times,version,T,M,land_mask)
!    integer(int64), dimension(N_LATS,N_LONS), intent(in) :: land_mask
!    real(int64), dimension(2,N_LATS,N_LONS),intent(in) :: T,M
!    integer(int64), intent(in) :: n_times, version
!    real(real64) :: time, days
!    integer(int64) :: ti
!
!    !Grid_vals *grid = xmalloc(sizeof(Grid_vals)) I think thi is unique to c
!    ! need construct grid statement
!    !temperature data
!
!    do ti = 1,n_times
!        time = (ti+1.)*DELTA_T
!        if (mod(time,TIME_OUTPUT_FREQ)<TIME_OUTPUT_TOL) then
!            days = T_OFFSET+time/(HOURS_PER_DAY*MINUTES_PER_HOUR*SECONDS_PER_MINUTE)
!            print*, "Days passed = %lg\n", days
!        end if
!        if (mod(time,DATA_OUTPUT_FREQ)<TIME_OUTPUT_TOL) then
!            !process_output_data(T,M,grid,N_LATS,N_LONS,time)
!        end if
!    end do
!
!
!end subroutine
!
!subroutine test
!
!  integer(fgsl_size_t), parameter :: ndim = 12
!  integer(fgsl_int), target :: v(ndim)
!  integer(fgsl_int), pointer :: p_slice(:), p_vec(:)
!  integer(fgsl_int) :: i, status
!  type(fgsl_vector_int) :: vec, slice
!
!  v(1:ndim) = [ (i, i=1,ndim) ]
!  slice = fgsl_vector_init(v(3:), stride=3_fgsl_size_t)
!  vec = fgsl_vector_init(v)
!  p_slice => fgsl_vector_to_fptr(slice)
!  p_vec => fgsl_vector_to_fptr(vec)
!  write(*, '(''Size of slice pointer is: '',i3)') size(p_slice)
!  write(*, '(''Size of complete pointer is: '',i3)') size(p_vec)
!  write(*, '(''Components: '',4i4)') p_slice(1:size(p_slice))
!  write(*, '(''Should be : '',4i4)') v(3::3)
!  v(6) = v(6) + 1
!  write(*, '(''Increase value of 2nd element of slice: '',2i4)') &
!       p_slice(2), p_vec(6)
!  call fgsl_vector_free(slice)
!  call fgsl_vector_free(vec)
!
!end subroutine test

end module mat_test
