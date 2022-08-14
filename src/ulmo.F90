
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
!real(real64), dimension(N_LATS,N_LONS) :: T_surf, T_deep,mass_flux_THETA,mass_flux_PHI,sv_mass_flux_PHI,sv_mass_flux_THETA
!real(real64),dimension(2,N_LATS,N_LONS) :: T, M, T_new
!integer(int64), dimension(N_LATS,N_LONS) :: land_mask
!real(real64) :: time
!
!integer(fgsl_size_t), parameter :: n = 25920
!real(fgsl_double), target :: v(n)
!
!
!integer(int64):: lon,lat,height,n_times
!type(fgsl_spmatrix) :: A
!type(fgsl_vector) :: B


!character(len=100) :: a


! COMBINING T_SURF AND T_DEEP INTO ONE 3D ARRAY FOR EASIER ANALYSIS !


! COMBINING M_THETA and M_PHI INTO ONE 3D ARRAY FOR EASIER ANALYSIS !
!mass_flux_THETA = calculate_mass_flux_THETA()
!mass_flux_PHI = calculate_mass_flux_PHI()
!sv_mass_flux_THETA =calculate_flow_sv_THETA()
!sv_mass_flux_PHI =calculate_flow_sv_PHI()
!
!M(1,:,:) = mass_flux_THETA
!M(2,:,:) = mass_flux_PHI

!! Writing mass flux outputs to files !!
!call write_file('output_data/sv_mass_flux_Phi.dat',sv_mass_flux_Phi,N_LATS,N_LONS)
!call write_file('output_data/mass_flux_PHI.dat',mass_flux_PHI,N_LATS,N_LONS)
!call write_file('output_data/sv_mass_flux_Theta.dat',sv_mass_flux_THETA,N_LATS,N_LONS)
!call write_file('output_data/mass_flux_Theta.dat',mass_flux_THETA,N_LATS,N_LONS)

!!! TESTING THE TRASNPORT AND TIME STEPPER WITH NO DEEP OCEAN !!!
!print"(a,i5)", 'Running Program...'
!call temp_after_time_step(N_STEPS=TIME_STEPS) ! Time_steps value is specified in NAMELIST.f90 file
!print"(a,i5)", 'Finished!'

!!! Matrix calculation testing !!!

!! Calculate_matrix !! ****WORKS****
!lat = 1
!lon = 1
!height =1
!A = calculate_matrix()
!print*,'ulmo main print A(1,1):', fgsl_spmatrix_get(A, 1_fgsl_size_t, 1_fgsl_size_t)


!! calculate_vector_b !! ****WORKS****
!call calculate_vector_b(T,v,B)
!print*, v
!b = fgsl_vector_init(v)
!print*,'1st value', v(1)

!! calculate_new_T !!
!call calculate_new_T(T,T_new)
!print*,DELTA_T
!write(400, '(''Should be : '',4F12.5)') v(3::3)
!print*, T_new()

!! t_stepper !!
!    real(real64),dimension(2,N_LATS,N_LONS) :: T
!    real(real64), dimension(N_LATS,N_LONS) :: T_surf,T_deep
!    integer(int64) :: n_times !, version
!    real(real64) :: time, days !,avg_T_surf_init,avg_T_surf_new
!    integer(int64) :: ti
!
!    T_surf= read_file(INITIAL_SURFACE_TEMP_DATA,N_LATS,N_LONS)
!    T_deep= read_file(INITIAL_DEEP_TEMP_DATA,N_LATS,N_LONS)
!
!    T(1,:,:) = T_surf
!    T(2,:,:) = T_deep
!
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

! Mass flux testing !
real(real64), dimension(:,:),allocatable :: u_wind,v_wind,sv_flow_phi,sv_flow_theta
integer(int64), dimension(:,:), allocatable :: land_mask

allocate(u_wind(N_LATS,N_LONS),v_wind(N_LATS,N_LONS),land_mask(N_LATS,N_LONS),sv_flow_phi(N_LATS,N_LONS)&
,sv_flow_theta(N_LATS,N_LONS))

u_wind = read_file_real(U_WIND_DATA,N_LATS,N_LONS)
v_wind = read_file_real(V_WIND_DATA,N_LATS,N_LONS)
land_mask = read_file_int(LAND_MASK_DATA,N_LATS,N_LONS)

call calculate_flow_sv_PHI(u_wind,v_wind,land_mask,sv_flow_phi)
call calculate_flow_sv_THETA(u_wind,v_wind,land_mask,sv_flow_theta)

call write_file('output_data/sv_flow_Phi.dat',sv_flow_phi,N_LATS,N_LONS)
call write_file('output_data/sv_flow_Theta.dat',sv_flow_theta,N_LATS,N_LONS)

deallocate(u_wind,v_wind,land_mask,sv_flow_phi,sv_flow_theta)





end program ULMO




