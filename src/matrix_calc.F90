module MATRIX_CALC
use NAMELIST
use DEGREE_TO_RADIAN
use HEIGHT_OF_SLAB
use DIV_M
use dA_da
use WRITE_READ_DATA
use Other_FLUXES
use fgsl
implicit none
public :: calculate_matrix_index,calculate_new_lon,calculate_matrix,calculate_new_T,calculate_vector_b
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
end function calculate_new_lon
!**************************************************************************************************************************************
! Calculates and sets the matrix elements for a 2 layer ocean with 144 longitude points and 90 latitude points.
! version determines whether the matrix is calculated for a no horizontal transport model, diffusion only model or the full Ekman model.
!**************************************************************************************************************************************
subroutine calculate_matrix()
    real ::d_theta, d_phi, s_dfsn, s_ekman,thickness,Aij, Aii
    real, dimension(N_LATS,N_LONS) :: dM_theta_dtheta,dM_phi_dphi,M_theta,M_phi,div_M,land_mask
    real, dimension(N_LATS,2)      :: lats_data_file
    real, dimension(N_LATS,1)      ::  lats,theta
    integer:: height,SURFACE,DEEP,lat,lon,N_DEPTHS,i,j
    SURFACE = 1
    DEEP = 2

    lats_data_file = read_file(LATS_FILE,N_LATS,2)
    lats(:,1)= lats_data_file(:,2)

    d_theta = deg_to_rad(DELTA_LAT)
    d_phi   = deg_to_rad(DELTA_LON)
    s_dfsn  = DELTA_T*D/(R_PLANET**2)
    do height=SURFACE,N_DEPTHS
        thickness = h_slab(height)
        s_ekman = DELTA_T/(RHO_WATER*thickness*R_PLANET)
        do lat =1, N_LATS
            theta(i,1) = deg_to_rad(lats(i,1))!need to read in latitude data
            do lon = 1,N_LONS
                Aii = 0.0
                i = calculate_matrix_index(lat,lon,height)
                ! central point - check if land
                if (land_mask(i,j) == 1) then
                    Aii=1
                    !fgsl_spmatrix_set(A,i,i,Aii)
                    cycle
                end if
                ! define mass flux for ekman terms
                Aij = 0.0
                ! central point !
                Aii = 1
                ! if statements for with transport versions
                !fgsl_spmatrix_set(A,i,i,Aii)
                ! if statements for with transport versions !

            end do
        end do
    end do
end subroutine
!**********************************************************************************************************************
! calculates and sets the vector b in matrix Ax=b. EMISSIVITY*SIGMA*T^4 is a blackbody emission from one of the layers,
! to ensure an energy exchange between the layers. 1.0/(RHO_WATER*C_V*H_S) prefactor converts a flux in W/m2 to K/s
!**********************************************************************************************************************
subroutine calculate_vector_b(T,land_mask)
    real,dimension(2,N_LATS,N_LONS),intent(in) :: T
    real, dimension(N_LATS,N_LONS), intent(in) :: land_mask
    real, dimension(N_LATS,N_LONS) :: F_a,F_c
    real :: b_hij,depth
    integer :: h, N_DEPTHS,i ,j,SURFACE,DEEP
    real, dimension(N_LATS,1) :: b

    SURFACE = 1 ! position in 3D temperatures array (1,:,:)
    DEEP = 2    ! position in 3D temperatures array (2,:,:)

    b_hij = 0.0

    F_a = calculate_F_a()
    F_c = calculate_F_c(T)
    N_DEPTHS = 2


    do h = 1, N_DEPTHS
        do i = 1,N_LATS
            do j = 1,N_LONS
                if (land_mask(i,j) == 1) then
                    b_hij = T(h,i,j)
                else
                    depth=h_slab(h)
                    if (h==SURFACE) then
                        b_hij = (DELTA_T/(RHO_WATER*C_V*depth))*(F_c(i,j)+F_a(i,j)-EMISSIVITY*SIGMA*(T(SURFACE,i,j))**4)*T(h,i,j)
                    else
                        b_hij = (DELTA_T/(RHO_WATER*C_V*depth))*(-F_c(i,j)+T(h,i,j))

                    end if

                end if
                b =  fgsl_vector_init(b,(h*N_LATS*N_LONS+(i*N_LONS+j)),b_hij)
            end do
        end do
    end do

end subroutine

subroutine calculate_new_T(T,land_mask)
    integer :: n, N_DEPTHS
    real, dimension(N_LATS,N_LONS), intent(in) :: T, land_mask
    real:: tol

    n = N_DEPTHS*n_lats*N_LONS
!    fgsl_spmatrix *A = fgsl_spmatrix_alloc(n,n)
!    fgsl_spmatrix *C
!    fgsl_vector *b = fgsl_vector_alloc(n)
!    fgsl_vector *x = fgsl_vector_alloc(n)
!    call calculate_matrix()
!    call calculate_vector_b(T,land_mask)
!    C = fgsl_spmatrix_ccs(A)

    ! now solve the system with the GMRES iterative solver
    tol = TOL




end subroutine
end module MATRIX_CALC
