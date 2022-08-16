!******************************************
! Module for Matrix calculations using fgsl
!******************************************
module MATRIX_CALC
use NAMELIST
use DEGREE_TO_RADIAN
use READ_DATA
use heat_fluxes
use HEIGHT_OF_SLAB
use fgsl
use dA_da
use div_m
use, intrinsic :: iso_fortran_env
implicit none
public :: calculate_matrix , calculate_vector_f_values , calculate_new_T
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
! At maximum or minimum latitutde, the adjacent latitude point above or bellow, respectively, is located at the same latitude,
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
subroutine calculate_matrix(land_mask,A)!,M,version)
  ! NO TRANSPORT !
  ! version: 0 = no transport, 1= Diffusion, 2 = Diffusion ekman
  !integer(int64), intent(in) :: version
  integer(int64), dimension(:,:), intent(in) :: land_mask
  !real(real64),dimension(:,:,:), intent(in) :: M
  !real(real64) :: thickness,s_ekman,theta,dM_theta_dtheta,dM_phi_dphi,M_theta,M_phi,div_M,d_theta,d_phi,s_dfsn
  !real(int64), dimension(:,:),allocatable :: lats_data,lats
  integer(int64):: lat,lon,height!,lon_j,lat_j
  integer(fgsl_size_t) ::  i !,j
  integer(fgsl_int) :: status
  real(fgsl_double) :: Aii,Aij
  type(fgsl_spmatrix),intent(out) :: A

  ! matrix size = 144x90x2+90x144+144 = 39024 (height*N_LATS*N_LONS+lat*N_LONS+lon)
  A = fgsl_spmatrix_alloc(25920_fgsl_size_t, 25920_fgsl_size_t)

  !allocate(lats_data(N_LATS,2),lats(N_LATS,1))

  !d_theta = deg_to_rad(DELTA_LAT)
  !d_phi = deg_to_rad(DELTA_LON)
  !s_dfsn = DELTA_T*D/(R_PLANET**2)

  do height = 0,N_DEPTHS-1
    !thickness = h_slab(height+1)
    !s_ekman = DELTA_T/(RHO_WATER*thickness*R_PLANET)
    do lat = 0,N_LATS-1
        !lats_data = read_file_real(LATS_FILE,N_LATS,2_int64) ! This can be converted to a grid input
        !lats(:,1)= lats_data(:,2) ! This can be converted to a grid input
        !theta = deg_to_rad(lats(lat+1,1))
        do lon = 0,N_LONS-1
            i = calculate_matrix_index(lon,lat,height)
            ! check if land !
            if(land_mask(lat+1,lon+1)==1) then
                Aii = 1
                status = fgsl_spmatrix_set(A,i,i,Aii)
                cycle ! ask about this
            end if

            !M_theta = M(1,lat+1,lon+1)
            !M_phi = M(2,lat+1,lon+1)
            !call dA_d_theta(M(1,:,:),lat+1,lon+1,dM_theta_dtheta) ! +1 as these are Fortran arrays
            !call dA_d_phi(M(2,:,:),lat+1,lon+1,dM_phi_dphi)
            !call calculate_div_M(lat+1,lon+1,M(2,:,:),M(1,:,:),div_M)
            Aij = 0.0
            ! central point !
            Aii = 1 ! no transport makes all diagonal components 1

!            if(version /= 0) then
!                Aii = Aii+2.*s_dfsn*(1./(d_theta**2)+1./((cos(theta)*d_phi)**2))
!                if(version == 2) then
!                    if(height == 0) then
!                        Aii = Aii + s_ekman*(dM_theta_dtheta-M_theta*tan(theta)+1./cos(theta)*dM_phi_dphi)
!                        if(div_M < 0.0) then
!                            Aii = Aii +(-DELTA_T/(RHO_WATER*thickness)*div_M)
!                        end if
!                    elseif(height == 1) then
!                        Aii = Aii + (-s_ekman)*(dM_theta_dtheta-M_theta*tan(theta)+1./cos(theta)*dM_phi_dphi)
!                        if(div_M > 0.0) then
!                            Aii = Aii+(+DELTA_T/(RHO_WATER*thickness)*div_M)
!                        end if
!                    end if
!                end if
!            end if

            status = fgsl_spmatrix_set(A, i, i, Aii) ! build the sparse matrix

            ! Neighbouring points !
!            if(version /= 0) then
!                ! neighbouring latitude - below !
!                if(lat == 0) then
!                    lon_j = calculate_new_lon(lon)
!                    lat_j = lat
!                    j = calculate_matrix_index(lon_j,lat_j,height)
!                else
!                    lon_j = lon
!                    lat_j = lat
!                    j = calculate_matrix_index(lon_j,lat_j,height)
!                end if
!
!                Aij = -s_dfsn*(1./(d_theta)**2+tan(theta)/(2*d_theta))
!
!                if(version == 2) then
!                    if(height == 0) then
!                        Aij = Aij + s_ekman*(-M_theta/(2.*d_theta))
!                    else
!                        Aij = Aij + (-s_ekman)*(-M_theta/(2.*d_theta))
!                    end if
!                end if
!
!                if(land_mask(lat_j+1,lon_j+1) == 1) then
!                    Aii = Aii + Aij
!                    Aij = 0.
!                end if
!
!                status = fgsl_spmatrix_set(A,i,j,Aij)
!
!                ! neighbouring longitude - left !
!
!                if(lon == 0) then
!                    lon_j = n_lons-1
!                    lat_j = lat
!                    j = calculate_matrix_index(lon_j,lat_j,height)
!                else
!                    lon_j = lon-1
!                    lat_j = lat
!                    j = calculate_matrix_index(lon_j,lat_j,height)
!                end if
!
!                Aij = -s_dfsn*1./(cos(theta)*d_phi)**2
!
!                if(version == 2) then
!                    if(height == 0) then
!                        Aij = Aij + s_ekman*(-M_phi/(2.*d_phi*cos(theta)))
!                    else
!                        Aij = Aij - s_ekman*(-M_phi/(2.*d_phi*cos(theta)))
!                    end if
!                end if
!
!                if(land_mask(lat_j+1,lon_j+1) == 1) then
!                    Aii = Aii + Aij
!                    Aij = 0.
!                end if
!
!                status = fgsl_spmatrix_set(A,i,j,Aij)
!
!                ! neighbouring longitude - right !
!
!                if(lon == N_LONS-1) then
!                    lon_j = 0
!                    lat_j = lat
!                    j = calculate_matrix_index(lon_j,lat_j,height)
!                else
!                    lon_j = lon+1
!                    lat_j = lat
!                    j = calculate_matrix_index(lon_j,lat_j,height)
!                end if
!
!                Aij = -s_dfsn*1./(cos(theta)*d_phi)**2
!
!                if(version == 2) then
!                    if(height == 0) then
!                        Aij = Aij + s_ekman*(M_phi/(2.*d_phi*cos(theta)))
!                    else
!                        Aij = Aij + (-s_ekman)*(M_phi/(2.*d_phi*cos(theta)))
!                    end if
!                end if
!
!                if(land_mask(lat_j+1,lon_j+1) == 1) then
!                    Aii = Aii + Aij
!                    Aij = 0.
!                end if
!
!                status = fgsl_spmatrix_set(A,i,j,Aij)
!                status = fgsl_spmatrix_set(A,i,i,Aii)
!
!                ! Interaction between deep and surface layer !
!
!                if(version == 2) then
!                    if(height == 0 .and. div_M > 0.0) then
!                        Aij = (-DELTA_T/(RHO_WATER*thickness)*div_M)
!                        j = calculate_matrix_index(lon,lat,1_int64)
!                        status = fgsl_spmatrix_set(A,i,j,Aij)
!                    elseif(height == 1 .and. div_M < 0.0) then
!                        Aij = (+DELTA_T/(RHO_WATER*thickness)*div_M)
!                        j = calculate_matrix_index(lon,lat,0_int64)
!                        status = fgsl_spmatrix_set(A,i,j,Aij)
!                    end if
!                end if

            !end if
        end do
    end do
  end do

  !deallocate(lats_data,lats)
  !call fgsl_spmatrix_free(A)

end subroutine calculate_matrix
!**********************************************************************************************************************
! calculates and sets the vector b in matrix Ax=b. EMISSIVITY*SIGMA*T^4 is a blackbody emission from one of the layers,
! to ensure an energy exchange between the layers. 1.0/(RHO_WATER*C_V*H_S) prefactor converts a flux in W/m2 to K/s
!**********************************************************************************************************************
subroutine calculate_vector_f_values(T,f_f,land_mask)
    real(real64),dimension(:,:,:),intent(in) :: T
    integer(int64), dimension(:,:), intent(in) :: land_mask
    real(real64), dimension(:,:),allocatable :: F_a,F_c
    real(real64), dimension(:),allocatable :: vec_data
    real(real64) :: depth,b_hij
    integer(int64) :: h,i ,j,SURFACE,DEEP
    ! fgsl !
    integer(fgsl_size_t), parameter :: n = 25920
    real(fgsl_double),dimension(:),intent(out) :: f_f
    integer(fgsl_size_t) :: vec_index
    real(fgsl_double) :: vec_val


    allocate(F_c(N_LATS,N_LONS),F_a(N_LATS,N_LONS),vec_data(n))

    SURFACE = 1 ! position in 3D temperatures array (1,:,:)
    DEEP = 2    ! position in 3D temperatures array (2,:,:)

    call calculate_F_a(F_a)
    call calculate_F_c(T,F_c)

    do h = 1,N_DEPTHS
        do i = 1, N_LATS
            do j= 1,N_LONS
                if(land_mask(i,j)==1) then
                    b_hij = T(h,i,j)
                else
                    depth = h_slab(h)
                    if (h==1) then
                        b_hij = DELTA_T/(RHO_WATER*C_V*depth)*(F_c(i,j)+F_a(i,j)-EMISSIVITY*SIGMA*((T(h,i,j))**4))+T(h,i,j)
                    else
                        b_hij = DELTA_T/(RHO_WATER*C_V*depth)*(-F_c(i,j))+T(h,i,j) !! double check this with Jake !!
                end if
                vec_index = (h-1)*N_LATS*N_LONS+(i-1)*N_LONS+(j-1) ! fgsl index
                vec_data(vec_index+1) = b_hij ! needs fortran index
                vec_val = b_hij ! making sure the type is fgsl_size_t
                f_f(vec_index+1) = vec_val

                end if
            end do
        end do
    end do

    deallocate(F_a,F_c,vec_data)

end subroutine calculate_vector_f_values
!!**************************************************************
!! Subroutine that uses FGSL sparse linear algebra to solve Ax=b
!!**************************************************************
subroutine calculate_new_T(T,land_mask)!,M,version)

    !integer(int64),intent(in) :: version
    integer(int64), dimension(:,:), intent(in) :: land_mask
    real(real64), dimension(:,:,:), intent(inout) :: T!,M
    integer(int64) :: h,i,j
    ! fgsl !
    type(fgsl_spmatrix):: A,C
    real(fgsl_double):: residual
    integer(fgsl_size_t),parameter:: n = 25920
    integer(fgsl_int) :: status
    integer(fgsl_size_t) :: iter = 0
    real(fgsl_double), parameter :: tol = 1.0e-6
    type(fgsl_vector) :: u,f
    type(fgsl_splinalg_itersolve_type), parameter :: S = fgsl_splinalg_itersolve_gmres
    type(fgsl_splinalg_itersolve) :: work
    real(fgsl_double),dimension(:),target,allocatable :: u_f,f_f

 !   solving equation: Au = f


    !A = fgsl_spmatrix_alloc(n,n)
    allocate(u_f(n),f_f(n))




    work =  fgsl_splinalg_itersolve_alloc(S,n,0_fgsl_size_t)

    ! initializing vector u
    u = fgsl_vector_init(u_f)

    ! vector check !
    !size_vec =  fgsl_vector_get_size(u)
    !print*, 'size of u =', size_vec
    !print*,'vector u= ', u_f



    f = fgsl_vector_init(f_f)

    call calculate_vector_f_values(T,f_f,land_mask)



    ! vector check !
    !print*,'f vector=',f_f
    !size_vec =  fgsl_vector_get_size(f)
    !print*, 'size of f=', size_vec


    call calculate_matrix(land_mask,A)
    ! Matrix check !
!    do i = 0_fgsl_size_t,25919_fgsl_size_t
!        print*, fgsl_spmatrix_get(A,i,i) ! These should all be 1
!        print*, fgsl_spmatrix_get(A,i,2_fgsl_size_t) ! These should all be 0 except 2,2
!    end do

    C = fgsl_spmatrix_compcol(A) ! compressed column format
    ! Matrix check !
!    do i = 0_fgsl_size_t,25919_fgsl_size_t
!        print*, fgsl_spmatrix_get(C,i,i) ! These should all be 1
!        print*, fgsl_spmatrix_get(C,i,2_fgsl_size_t) ! These should all be 0 except 2,2
!    end do


!    !initial guess x = 0
    u_f = 0

    do
        status = fgsl_splinalg_itersolve_iterate(C, f, tol, u, work)

         !print out residual norm ||A*x-b||
        residual = fgsl_splinalg_itersolve_normr(work)
        !write(output_unit, '(A,I2,A,G15.6)') 'iter ', iter, ' residual = ', residual

        if (status == FGSL_SUCCESS) then
        !    write(output_unit, '(A)') 'Converged'
        endif
        iter = iter + 1
        if (status /= FGSL_CONTINUE .or. iter >= MAX_ITER) exit
    end do

!   output solution !

    do h = 1,N_DEPTHS
        do i = 1, N_LATS
            do j = 1,N_LONS
                T(h,i,j) =  u_f(h*N_LATS*N_LONS+(i*N_LONS)+j)
            end do
        end do
    end do

    !print*,'T_new(1,1,1) = ',T_new(1,1,1)
    !print*,'T_old(1,1,1)=', T(1,1,1)

    deallocate(u_f,f_f)
    call fgsl_splinalg_itersolve_free(work)
    call fgsl_spmatrix_free(A)
    call fgsl_spmatrix_free(C)
    call fgsl_vector_free(f)
    call fgsl_vector_free(u)





end subroutine



end module MATRIX_CALC
