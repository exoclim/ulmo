program ULMO

    !**This version of Ulmo uses the fgsl linear algebra solver to solve differential equations**!

    use MASS_FLUX, only: calculate_surface_stress_PHI,calculate_surface_stress_THETA,calculate_flow_sv_THETA &
    ,calculate_mass_flux_THETA,calculate_mass_flux_PHI,calculate_flow_sv_PHI
    use READ_DATA, only: read_file_real,read_file_int,write_file
    use Constants
    use heat_fluxes, only: calc_Q_flux,calculate_F_a,calculate_F_c
    use MATRIX_CALC, only: calculate_matrix_index,calculate_matrix,calculate_vector_f_values
    use fgsl
    use memory_usage
    use, intrinsic :: iso_fortran_env
    !**TESTING**!
    use process_output_data
implicit none

!*********************
!! ULMO MAIN SCRIPT !!
!*********************

!**Mass flux variables**!
real(real64), dimension(:,:),allocatable :: u_wind,v_wind,tau_PHI,tau_THETA,mass_flux_THETA &
                                            ,mass_flux_PHI,sv_flow_PHI,sv_flow_THETA
real(real64), dimension(:,:,:), allocatable :: M

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
!real(fgsl_double),dimension(:),allocatable :: f_f,u_f ! could make these a target
real(fgsl_double),dimension(1:n),target :: f_f, u_f ! target is used in git hub example
type(fgsl_vector) :: f,u
integer(int64) :: h,i,j
real(fgsl_double):: residual
integer(fgsl_int) :: status
integer(fgsl_size_t) :: iter = 0
!real(fgsl_double), parameter :: tol = 1.0e-6
type(fgsl_splinalg_itersolve_type), parameter :: S = fgsl_splinalg_itersolve_gmres
type(fgsl_splinalg_itersolve) :: work

!**Time stepping varaibles**!
integer(int64) :: n_step,n_times

!**Memory usage varaiables**!
integer :: valueRSS
real(real64),dimension(:), allocatable :: Rss_data

!**CPU Time varaiables**!
real(real64) :: start, finish
real(real64),dimension(:,:),allocatable:: cpu_data

!**SYSTEM ANALYSIS**!
call cpu_time(start)
allocate(Rss_data(2005))
call system_mem_usage(valueRSS)
Rss_data(1) = valueRSS
write (*,"(a18,i5)") 'valueRSS start =',valueRSS

!**Land Mask allocation and load**!
allocate(land_mask(N_LATS,N_LONS))
land_mask = read_file_int(LAND_MASK_DATA,N_LATS,N_LONS)

!**Calculating surface stresses from u and v wind**!
allocate(u_wind(N_LATS,N_LONS),v_wind(N_LATS,N_LONS),tau_PHI(N_LATS,N_LONS),tau_THETA(N_LATS,N_LONS))

u_wind = read_file_real(U_WIND_DATA,N_LATS,N_LONS)
v_wind = read_file_real(V_WIND_DATA,N_LATS,N_LONS)

call calculate_surface_stress_PHI(u_wind,tau_PHI)
call calculate_surface_stress_THETA(v_wind,tau_THETA)

deallocate(u_wind,v_wind)

!**Calculating mass fluxes from surface stresses**!
allocate(lats_data(N_LATS,2),lats(N_LATS,1))

allocate(mass_flux_THETA(N_LATS,N_LONS),mass_flux_PHI(N_LATS,N_LONS))

lats_data = read_file_real(LATS_FILE,N_LATS,2_int64)
lats(:,1)= lats_data(:,2) ! Second column in lats data file includes the latitude points

deallocate(lats_data)

call calculate_mass_flux_THETA(lats,tau_PHI,tau_THETA,land_mask,mass_flux_THETA)
call calculate_mass_flux_PHI(lats,tau_PHI,tau_THETA,land_mask,mass_flux_PHI)

deallocate(tau_PHI,tau_THETA)

!**Calculating mass flux in Sverdrups**!!
allocate(sv_flow_PHI(N_LATS,N_LONS),sv_flow_THETA(N_LATS,N_LONS))

call calculate_flow_sv_PHI(mass_flux_PHI,sv_flow_PHI)
call calculate_flow_sv_THETA(lats,mass_flux_THETA,sv_flow_THETA)

!**Writing mass flux initial output to files**!
!call write_file('output_data/sv_flow_PHI_init.dat',sv_flow_PHI,N_LATS,N_LONS)
!call write_file('output_data/sv_flow_Theta_init.dat',sv_flow_THETA,N_LATS,N_LONS)

allocate(M(2,N_LATS,N_LONS))

M(1,:,:) = mass_flux_PHI(:,:) !PHI = 1
M(2,:,:) = mass_flux_THETA(:,:) !THETA = 2

deallocate(mass_flux_PHI,mass_flux_THETA,lats,sv_flow_PHI,sv_flow_THETA)

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

!**SYSTEM ANALYSIS**!
call system_mem_usage(valueRSS)
write (*,"(a49,i5)") 'valueRSS before allocating vectors and matrix =',valueRSS
Rss_data(2) = valueRSS

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


!**SYSTEM ANALYSIS**!
call system_mem_usage(valueRSS)
write (*,"(a49,i5)") 'valueRSS before time stepping =',valueRSS
Rss_data(3) = valueRSS
call cpu_time(finish)
print '("CPU time before time stepping = ",f6.3," seconds.")',finish-start
allocate(cpu_data(n_times,1_int64))


!**CALCULATING NEW TEMP**!
n_times = 2000 !number of time steps 20 steps = 10 days
do n_step = 1,n_times

    !**SYSTEM ANALYSIS**!
    call cpu_time(start)
    call system_mem_usage(valueRSS)
    write(*,"(a40,i5,a1,i6)") 'valueRSS before solving F_c etc',n_step,'=',valueRSS
    Rss_data(n_step+2) = valueRSS


    call calculate_F_c(T,F_c)
    !call calc_Q_flux(T,upward_Q_flux,F_net_sw_down,F_lw_down,F_latent_up,F_sensible_up)
    call calculate_vector_f_values(land_mask,T,f_f,F_a,F_c)

    !**SYSTEM ANALYSIS**!
    call system_mem_usage(valueRSS)
    write(*,"(a40,i5,a1,i6)") 'valueRSS before solving AU=f at n_step',n_step,'=',valueRSS
    Rss_data(n_step+3) = valueRSS

    !**Solving eqaution Au = f using fgsl**!

    work =  fgsl_splinalg_itersolve_alloc(S,n,0_fgsl_size_t)

    do
        status = fgsl_splinalg_itersolve_iterate(C, f, tol, u, work)

        ! print out residual norm ||A*x-b||
        residual = fgsl_splinalg_itersolve_normr(work)
        !write(output_unit, '(A,I2,A,G15.6)') 'iter ', iter, ' residual = ', residual

    !    if (status == FGSL_SUCCESS) then
    !        !write(output_unit, '(A)') 'Converged'
    !    endif
        iter = iter + 1
        if (status /= FGSL_CONTINUE .or. iter >= MAX_ITER) exit
    end do

    !**output solution**!

    do h = 0,N_DEPTHS-1
        do i = 0, N_LATS-1
            do j = 0,N_LONS-1
                T(h+1,i+1,j+1) = u_f(calculate_matrix_index(j,i,h)+1)
            end do
        end do
    end do


    call fgsl_splinalg_itersolve_free(work) !**CRUCIAL THAT THIS STATEMENT IS INSIDE THE LOOP**!

    !**SYSTEM ANALYSIS**!
    call cpu_time(finish)
    print '("CPU time at step = ",I6.3,f6.3," seconds.")',n_step,finish-start
    cpu_data(n_step,1) = finish-start

end do

!**SYSTEM ANALYSIS**!
call system_mem_usage(valueRSS)
write (*,"(a31,i5)")  'valueRSS before deallocating =',valueRSS
Rss_data(n_times+4) = valueRSS

!**TESTING**!
call process_output(T,upward_Q_flux)

!**Deallocation of Matrices and vectors**!
call fgsl_vector_free(f) ! deallocate vector
call fgsl_spmatrix_free(A) ! deallocated matrix
call fgsl_spmatrix_free(C)
deallocate(F_net_sw_down,F_lw_down,F_latent_up,F_sensible_up)
deallocate(T,M)
deallocate(land_mask)
deallocate(upward_Q_flux,F_a,F_c)


!**SYSTEM ANALYSIS**!
call system_mem_usage(valueRSS)
write (*,"(a16,i5)")  'valueRSS end =',valueRSS
Rss_data(n_times+5) = valueRSS
call write_file('output_data/rss.dat',Rss_data,n_times+5,1_int64)
print*, 'rss_size =', n_times+5
call write_file('output_data/cpu_data.dat',cpu_data,n_times,1_int64)
deallocate(Rss_data)
deallocate(cpu_data)

!integer(int64) :: n_times ! time stepper testing
!!real(fgsl_double),dimension(:),target,allocatable :: f_f ! vector testing
!!type(fgsl_vector) :: f ! vector testing
!!integer(fgsl_size_t),parameter:: n = 25920 ! vector testing
!real(real64) :: time, days !,avg_T_surf_init,avg_T_surf_new ! timestepper testing
!integer(int64) :: ti ! time stepper testing
!
!n_times = 20
!
!allocate(T(2,N_LATS,N_LONS),T_surf(N_LATS,N_LONS),T_deep(N_LATS,N_LONS),T_surf_out(n_times,1))
!
!!allocate(f_f(n)) ! vector testing
!
!T_surf= read_file_real(INITIAL_SURFACE_TEMP_DATA,N_LATS,N_LONS)
!T_deep= read_file_real(INITIAL_DEEP_TEMP_DATA,N_LATS,N_LONS)
!
!T(1,:,:) = T_surf
!T(2,:,:) = T_deep

!! Calculate_matrix !! ****WORKS****
!lat = 1
!lon = 1
!height =1
!A = calculate_matrix()
!print*,'ulmo main print A(1,1):', fgsl_spmatrix_get(A, 1_fgsl_size_t, 1_fgsl_size_t)


!! calculate_vector_f !! ****WORKS****
!f = fgsl_vector_init(f_f)
!call calculate_vector_f_values(T,f_f)
!print*, f_f
!deallocate(f_f)
!b = fgsl_vector_init(v)
!print*,'1st value', v(1)

!! calculate_new_T !! ****WORKS****
!print*, 'T_1(1,1,1) = ', T(1,2,7)
!call calculate_new_T(T)
!print*, 'T_2(1,1,1) = ', T(1,2,7)
!call calculate_new_T(T)
!print*, 'T_3(1,1,1)=',T(1,2,7)
!ti = 3
!time = (ti+1)*DELTA_T

!do ti = 1,n_times
!    call calculate_new_T(T)
!    T_surf_out(ti,1) = sum(T(1,:,:))/(144*90)
!end do
!print*, T_surf_out
!call write_file('output_data/ProCb/T_surf.dat',T_surf_out,n_times,1_int64)


!! t_stepper !!




    !Grid_vals *grid = xmalloc(sizeof(Grid_vals)) I think thi is unique to c
    ! need construct grid statement
    !temperature data

!do ti = 1,n_times
!    time = (ti+1.)*DELTA_T
!    call calculate_new_T(T)
!    if (mod(time,TIME_OUTPUT_FREQ)<TIME_OUTPUT_TOL) then
!        days = T_OFFSET+time/(HOURS_PER_DAY*MINUTES_PER_HOUR*SECONDS_PER_MINUTE)
!        print*, "Days passed =", days,'Avg Surf Temp = ', sum(T(1,:,:))/(144*90)
!    end if
!
!    if (mod(time,DATA_OUTPUT_FREQ)<TIME_OUTPUT_TOL) then
!        call process_output(T,time)!,M)
!    end if
!
!end do

!!! MEMORY CHECKING !!!
!call system_mem_usage(1)


!deallocate(T,T_surf,T_deep,T_surf_out)




end program ULMO

