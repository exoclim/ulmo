
program ULMO
    use MASS_FLUX
    use READ_DATA
    use NAMELIST
    use MATRIX_CALC
    use fgsl
    use process_output_data
    use, intrinsic :: iso_fortran_env
implicit none
!*********************
!! ULMO MAIN SCRIPT !!
!*********************

!! Version: 0 = no transport, 1 = diffusion only, 2 = Ekman and diffusion !!

real(real64), dimension(:,:),allocatable :: u_wind,v_wind,mass_flux_phi,mass_flux_theta,sv_flow_phi,sv_flow_theta
integer(int64), dimension(:,:), allocatable :: land_mask
real(real64),dimension(:,:,:),allocatable :: T, M
real(real64), dimension(:,:),allocatable :: T_surf,T_deep,T_surf_out
!type(fgsl_spmatrix) :: A !matrix testing

integer(int64) :: n_times ! time stepper testing
!real(fgsl_double),dimension(:),target,allocatable :: f_f ! vector testing
!type(fgsl_vector) :: f ! vector testing
!integer(fgsl_size_t),parameter:: n = 25920 ! vector testing
!real(real64) :: time, days !,avg_T_surf_init,avg_T_surf_new ! timestepper testing
integer(int64) :: ti ! time stepper testing

n_times = 3

allocate(T(2,N_LATS,N_LONS),M(2,N_LATS,N_LONS),T_surf(N_LATS,N_LONS),T_deep(N_LATS,N_LONS),T_surf_out(n_times,1) &
        ,u_wind(N_LATS,N_LONS),v_wind(N_LATS,N_LONS),land_mask(N_LATS,N_LONS),mass_flux_phi(N_LATS,N_LONS)       &
        , mass_flux_theta(N_LATS,N_LONS),sv_flow_phi(N_LATS,N_LONS),sv_flow_theta(N_LATS,N_LONS))

!allocate(f_f(n)) ! vector testing

land_mask = read_file_int(LAND_MASK_DATA,N_LATS,N_LONS)

T_surf= read_file_real(INITIAL_SURFACE_TEMP_DATA,N_LATS,N_LONS)
T_deep= read_file_real(INITIAL_DEEP_TEMP_DATA,N_LATS,N_LONS)

! Surface = 1 Deep = 2
T(1,:,:) = T_surf
T(2,:,:) = T_deep

u_wind = read_file_real(U_WIND_DATA,N_LATS,N_LONS)
v_wind = read_file_real(V_WIND_DATA,N_LATS,N_LONS)
call calculate_mass_flux_PHI(u_wind,v_wind,land_mask,mass_flux_phi)
call calculate_mass_flux_THETA(u_wind,v_wind,land_mask,mass_flux_theta)

! Theta = 1 Phi = 2
M(1,:,:) = mass_flux_theta
M(2,:,:) = mass_flux_phi

! Initial mass flux tests !
call calculate_flow_sv_PHI(u_wind,v_wind,land_mask,sv_flow_phi)
call calculate_flow_sv_THETA(u_wind,v_wind,land_mask,sv_flow_theta)
call write_file('output_data/sv_flow_Phi.dat',sv_flow_phi,N_LATS,N_LONS)
call write_file('output_data/sv_flow_Theta.dat',sv_flow_theta,N_LATS,N_LONS)




!! Calculate_matrix !! ****WORKS**** For version = 0
!lat = 1
!lon = 1
!height =1
!call calculate_matrix(A)
!print*, fgsl_spmatrix_get(A,91_fgsl_size_t,91_fgsl_size_t)
!print*,'ulmo main print A(1,1):', fgsl_spmatrix_get(A, 1_fgsl_size_t, 1_fgsl_size_t)


!! calculate_vector_f !! ****WORKS****
!f = fgsl_vector_init(f_f)
!call calculate_vector_f_values(T,f_f)
!print*, f_f(1)
!deallocate(f_f)



!! calculate_new_T !! ****WORKS****
!print*, 'T_1(1,1,1) = ', T(1,2,7)
!call calculate_new_T(T)
!print*, 'T_2(1,1,1) = ', T(1,2,7)
!call calculate_new_T(T)
!print*, 'T_3(1,1,1)=',T(1,2,7)
!ti = 3
!time = (ti+1)*DELTA_T

do ti = 1,n_times
    !print*, 'T_surf_avg =' ,sum(T(1,:,:))/(144*90)
    print*, 'T(1,1,1) = ', T(1,1,1)
    call calculate_new_T(T)!,land_mask)
end do
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


deallocate(T,M,T_surf,T_deep,T_surf_out,u_wind,v_wind,mass_flux_phi,mass_flux_theta,land_mask &
            , sv_flow_phi,sv_flow_theta)




end program ULMO




