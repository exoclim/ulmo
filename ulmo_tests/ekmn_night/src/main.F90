program main

    !**This version of Ulmo uses the fgsl linear algebra solver to solve differential equations**!
    use READ_DATA, only: read_file_real,read_file_int,write_file
    use Constants
    use heat_fluxes, only: calc_Q_flux,calculate_F_a,calculate_F_c
    use MATRIX_CALC, only: calculate_matrix_index,calculate_matrix,calculate_vector_f_values
    use fgsl
    use, intrinsic :: iso_fortran_env
    use process_output_data
implicit none

!*********************
!! ULMO MAIN SCRIPT !!
!*********************

!**Heat flux variables**!
real(real64), dimension(:,:), allocatable :: F_net_sw_down, F_lw_down,F_latent_up,F_sensible_up
real(real64), dimension(:,:,:), allocatable :: T
real(real64), dimension(:,:), allocatable :: upward_Q_flux, F_a, F_c
real(real64), dimension(:,:), allocatable :: T_surf_init,T_deep_init

!**Land mask and coordinate variables**!
integer(int64), dimension(:,:), allocatable :: land_mask
real(real64), dimension(:,:),allocatable :: lats,lats_data

!**Matrix/Vector fgsl variables**!
type(fgsl_spmatrix) :: A,C
integer(fgsl_size_t), parameter :: n = 25920
real(fgsl_double),dimension(1:n),target :: f_f, u_f ! target is used in git hub example
type(fgsl_vector) :: f,u
integer(int64) :: h,i,j
real(fgsl_double):: residual
integer(fgsl_int) :: status
integer(fgsl_size_t) :: iter = 0
real(fgsl_double), parameter :: tol = 1.0e-6
type(fgsl_splinalg_itersolve_type), parameter :: S = fgsl_splinalg_itersolve_gmres
type(fgsl_splinalg_itersolve) :: work

!**Time stepping varaibles**!
integer(int64) :: n_step,n_times
real(real64) :: days,time



!**Land Mask and coords**!
allocate(land_mask(N_LATS,N_LONS))
land_mask = read_file_int(LAND_MASK_DATA,N_LATS,N_LONS)

allocate(lats_data(N_LATS,2),lats(N_LATS,1))

lats_data = read_file_real(LATS_FILE,N_LATS,2_int64)
lats(:,1)= lats_data(:,2) ! Second column in lats data file includes the latitude points

deallocate(lats_data)


!**Heat flux calculations**!
allocate(F_net_sw_down(N_LATS,N_LONS),F_lw_down(N_LATS,N_LONS),F_latent_up(N_LATS,N_LONS),F_sensible_up(N_LATS,N_LONS))
allocate(T(2,N_LATS,N_LONS))
allocate(upward_Q_flux(N_LATS,N_LONS),F_a(N_LATS,N_LONS),F_c(N_LATS,N_LONS))
allocate(T_surf_init(N_LATS,N_LONS),T_deep_init(N_LATS,N_LONS))


F_net_sw_down = read_file_real(SW_FLUX_NET_DOWN_DATA,N_LATS,N_LONS)
F_lw_down     = read_file_real(LW_FLUX_DOWN_DATA,N_LATS,N_LONS)
F_latent_up   = read_file_real(LATENT_UP_FLUX_DATA,N_LATS,N_LONS)
F_sensible_up = read_file_real(SENSIBLE_UP_FLUX_DATA,N_LATS,N_LONS)
T_surf_init   = read_file_real(INITIAL_SURFACE_TEMP_DATA,N_LATS,N_LONS)
T_deep_init   = read_file_real(INITIAL_DEEP_TEMP_DATA,N_LATS,N_LONS)
T(1,:,:) = T_surf_init(:,:)
T(2,:,:) = T_deep_init(:,:)

deallocate(T_surf_init,T_deep_init)


!**Calculating F_a**!
call calculate_F_a(F_net_sw_down,F_lw_down,F_latent_up,F_sensible_up,F_a) ! Only depends on um inputs

!**Initialsing vector f**!
f = fgsl_vector_init(f_f)


!**Initialising vector u**!
u = fgsl_vector_init(u_f)
u_f = 0. !initial guess x = 0

!**Calculating Matrix**!
call calculate_matrix(land_mask,A) ! allocated using subroutine
C = fgsl_spmatrix_compcol(A) ! compressed column format



!**CALCULATING NEW TEMP USING FGSL**!
n_times = 2000 !number of time steps 20 steps = 10 days
work =  fgsl_splinalg_itersolve_alloc(S,n,0_fgsl_size_t)
do n_step = 1,n_times

    call calculate_F_c(T,F_c)
    !call calc_Q_flux(T,upward_Q_flux,F_net_sw_down,F_lw_down,F_latent_up,F_sensible_up)
    call calculate_vector_f_values(land_mask,T,f_f,F_a,F_c)


    !**Solving eqaution Au = f using fgsl**!
    do
        status = fgsl_splinalg_itersolve_iterate(C, f, tol, u, work)

        ! print out residual norm ||A*x-b||
        residual = fgsl_splinalg_itersolve_normr(work)
        write(output_unit, '(A,I2,A,G15.6)') 'iter ', iter, ' residual = ', residual

        if (status == FGSL_SUCCESS) then
            !write(output_unit, '(A)') 'Converged'
        endif
        iter = iter + 1
        if (status /= FGSL_CONTINUE .or. iter >= 10) exit
    end do

    !**output solution**!

    do h = 0,N_DEPTHS-1
        do i = 0, N_LATS-1
            do j = 0,N_LONS-1
                T(h+1,i+1,j+1) = u_f(calculate_matrix_index(j,i,h)+1)
            end do
        end do
    end do



    !**Outputting data**!
    time = n_step*DELTA_T
    if (mod(time,TIME_OUTPUT_FREQ)<TIME_OUTPUT_TOL) then
        days = T_OFFSET+time/(HOURS_PER_DAY*MINUTES_PER_HOUR*SECONDS_PER_MINUTE)
        !print*, "Days passed =", days,'Avg Surf Temp = ', sum(T(1,:,:))/(144*90)
        print*, 'Days passed = ', days
    end if

    if (mod(time,DATA_OUTPUT_FREQ)<TIME_OUTPUT_TOL) then
        call process_output(T,upward_Q_flux,time)
    end if


end do


!**Deallocation of Matrices and vectors**!
call fgsl_splinalg_itersolve_free(work)
call fgsl_vector_free(f) ! deallocate vector
call fgsl_spmatrix_free(A) ! deallocated matrix
call fgsl_spmatrix_free(C)
deallocate(F_net_sw_down,F_lw_down,F_latent_up,F_sensible_up)
deallocate(T)
deallocate(land_mask)
deallocate(upward_Q_flux,F_a,F_c)


end program main
