
program ULMO
    use MASS_FLUX
    use READ_DATA
    use NAMELIST
    use MATRIX_CALC
    use fgsl
    !use time_stepper
    use process_output_data
    use, intrinsic :: iso_fortran_env
implicit none
!*********************
!! ULMO MAIN SCRIPT !!
!*********************
!!! Matrix calculation testing !!!
!! Mass flux testing !!
!real(real64), dimension(:,:),allocatable :: u_wind,v_wind,sv_flow_phi,sv_flow_theta
!integer(int64), dimension(:,:), allocatable :: land_mask
!
!allocate(u_wind(N_LATS,N_LONS),v_wind(N_LATS,N_LONS),land_mask(N_LATS,N_LONS),sv_flow_phi(N_LATS,N_LONS)&
!,sv_flow_theta(N_LATS,N_LONS))
!
!u_wind = read_file_real(U_WIND_DATA,N_LATS,N_LONS)
!v_wind = read_file_real(V_WIND_DATA,N_LATS,N_LONS)
!land_mask = read_file_int(LAND_MASK_DATA,N_LATS,N_LONS)
!
!call calculate_flow_sv_PHI(u_wind,v_wind,land_mask,sv_flow_phi)
!call calculate_flow_sv_THETA(u_wind,v_wind,land_mask,sv_flow_theta)
!
!call write_file('output_data/sv_flow_Phi.dat',sv_flow_phi,N_LATS,N_LONS)
!call write_file('output_data/sv_flow_Theta.dat',sv_flow_theta,N_LATS,N_LONS)
!
!deallocate(u_wind,v_wind,land_mask,sv_flow_phi,sv_flow_theta)

real(real64),dimension(:,:,:),allocatable :: T
real(real64), dimension(:,:),allocatable :: T_surf,T_deep
!real(fgsl_double),dimension(:),target,allocatable :: f_f ! vector testing
!type(fgsl_vector) :: f ! vector testing
!integer(fgsl_size_t),parameter:: n = 25920 ! vector testing

!integer(int64) :: n_times !, version
!real(real64) :: time, days !,avg_T_surf_init,avg_T_surf_new
!integer(int64) :: ti

allocate(T(2,N_LATS,N_LONS),T_surf(N_LATS,N_LONS),T_deep(N_LATS,N_LONS))

!allocate(f_f(n)) ! vector testing

T_surf= read_file_real(INITIAL_SURFACE_TEMP_DATA,N_LATS,N_LONS)
T_deep= read_file_real(INITIAL_DEEP_TEMP_DATA,N_LATS,N_LONS)

T(1,:,:) = T_surf
T(2,:,:) = T_deep

!! Calculate_matrix !! ****WORKS****
!lat = 1
!lon = 1
!height =1
!A = calculate_matrix()
!print*,'ulmo main print A(1,1):', fgsl_spmatrix_get(A, 1_fgsl_size_t, 1_fgsl_size_t)


!! calculate_vector_f !! ****WORKS****
!f = fgsl_vector_init(f_f)
!call calculate_vector_f(T,f_f)
!print*, f_f
!deallocate(f_f)
!b = fgsl_vector_init(v)
!print*,'1st value', v(1)

!! calculate_new_T !!
print*, 'T_1(1,1,1) = ', T(1,1,1)
call calculate_new_T(T)
print*, 'T_2(1,1,1) = ', T(1,1,1)
call calculate_new_T(T)
print*, 'T_3(1,1,1)=',T(1,1,1)

!print*,DELTA_T
!write(400, '(''Should be : '',4F12.5)') v(3::3)
!print*, T_new()

!! t_stepper !!

!    n_times = 1200000
!
!
!    !Grid_vals *grid = xmalloc(sizeof(Grid_vals)) I think thi is unique to c
!    ! need construct grid statement
!    !temperature data
!!    time = 0
!!    days = 0
!!    allocate(T(2,N_LATS,N_LONS))
!    do ti = 1,n_times
!        time = (ti+1.)*DELTA_T
!        print*,ti, T(1,1,1)
!        call calculate_new_T(T)
!        if (mod(time,TIME_OUTPUT_FREQ)<TIME_OUTPUT_TOL) then
!            days = T_OFFSET+time/(HOURS_PER_DAY*MINUTES_PER_HOUR*SECONDS_PER_MINUTE)
!            print*, "Days passed = %lg\n", days
!        end if
!
!        if (mod(time,DATA_OUTPUT_FREQ)<TIME_OUTPUT_TOL) then
!            call process_output(T,time)!,M)
!        end if
!
!!time = 1200000.
!!call process_output(T,time)
!        T = T
!    end do


deallocate(T,T_surf,T_deep)




end program ULMO




