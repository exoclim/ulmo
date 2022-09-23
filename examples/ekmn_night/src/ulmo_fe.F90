program ULMO

    !****************************************************************************************!
    !**This version of Ulmo uses the Forward Euler method to solving differential equations**!
    !****************************************************************************************!

    use MASS_FLUX, only: calculate_surface_stress_PHI,calculate_surface_stress_THETA,calculate_flow_sv_THETA &
    ,calculate_mass_flux_THETA,calculate_mass_flux_PHI,calculate_flow_sv_PHI
    use READ_DATA, only: read_file_real,read_file_int,write_file
    use Constants
    use heat_fluxes, only: calc_Q_flux,calculate_F_a,calculate_F_c
    use HEIGHT_OF_SLAB
    use DEGREE_TO_RADIAN
    use dA_da,only: calculate_dA_d_theta, calculate_dA_d_phi
    use div_m,only: calculate_div_M
    use calc_new_T_fe,only: calc_new_T_surf_diff, calc_new_T_deep_diff,calc_new_T_surf_notrns, &
                            calc_new_T_deep_notrns,calc_new_T_surf_ekman,calc_new_T_deep_ekman
    use, intrinsic :: iso_fortran_env
    use process_output_data

implicit none

!*********************
!! ULMO MAIN SCRIPT !!
!*********************

!**Mass flux variables**!
real(real64), dimension(:,:),allocatable :: u_wind,v_wind,tau_PHI,tau_THETA,mass_flux_THETA &
                                            ,mass_flux_PHI,sv_flow_PHI,sv_flow_THETA,Vert_mass_flux
real(real64), dimension(:,:,:), allocatable :: M
real(real64) :: div_m


!**Land mask and coordinate variables**!
integer(int64), dimension(:,:), allocatable :: land_mask
real(real64), dimension(:,:),allocatable :: lats,lats_data
real(real64) :: d_phi,d_theta,theta,thickness

!**Heat flux variables**!
real(real64), dimension(:,:), allocatable :: F_net_sw_down, F_lw_down,F_latent_up,F_sensible_up
real(real64), dimension(:,:,:), allocatable :: T,T_new
real(real64), dimension(:,:), allocatable :: upward_Q_flux, F_a, F_c
real(real64), dimension(:,:), allocatable :: T_surf_init,T_deep_init


!**Forward Euler/time stepping variables**!
integer(int64) :: n_times,n_step,h,i,j,version,diff_coef
real(real64) :: days,time,s_dfsn,s_ekmn,s_notrns



!**Land Mask and coordinates allocation and load **!
allocate(land_mask(N_LATS,N_LONS))
land_mask = read_file_int(LAND_MASK_DATA,N_LATS,N_LONS)

allocate(lats_data(N_LATS,2),lats(N_LATS,1))
lats_data = read_file_real(LATS_FILE,N_LATS,2_int64)
lats(:,1)= lats_data(:,2) ! Second column in lats data file includes the latitude points
deallocate(lats_data)

!**Calculating surface stresses from u and v wind**!
allocate(u_wind(N_LATS,N_LONS),v_wind(N_LATS,N_LONS),tau_PHI(N_LATS,N_LONS),tau_THETA(N_LATS,N_LONS))

u_wind = read_file_real(U_WIND_DATA,N_LATS,N_LONS)
v_wind = read_file_real(V_WIND_DATA,N_LATS,N_LONS)

call calculate_surface_stress_PHI(u_wind,tau_PHI)
call calculate_surface_stress_THETA(v_wind,tau_THETA)

deallocate(u_wind,v_wind)

!**Calculating mass fluxes from surface stresses**!
allocate(mass_flux_THETA(N_LATS,N_LONS),mass_flux_PHI(N_LATS,N_LONS))

call calculate_mass_flux_THETA(lats,tau_PHI,tau_THETA,land_mask,mass_flux_THETA)
call calculate_mass_flux_PHI(lats,tau_PHI,tau_THETA,land_mask,mass_flux_PHI)

deallocate(tau_PHI,tau_THETA)

!**Calculating mass flux in Sverdrups**!!
allocate(sv_flow_PHI(N_LATS,N_LONS),sv_flow_THETA(N_LATS,N_LONS))

call calculate_flow_sv_PHI(mass_flux_PHI,sv_flow_PHI)
call calculate_flow_sv_THETA(lats,mass_flux_THETA,sv_flow_THETA)


allocate(M(2,N_LATS,N_LONS))

M(1,:,:) = mass_flux_THETA(:,:) !THETA = 1
M(2,:,:) = mass_flux_PHI(:,:) !PHI = 2

deallocate(mass_flux_PHI,mass_flux_THETA)

!**Calculating vertical mass flux**!
allocate(Vert_mass_flux(N_LATS,N_LONS))

do i = 1,N_LATS
    do j =1,N_LONS
        call calculate_div_M(i,j,theta,M(2,:,:),M(1,:,:),div_m)
        Vert_mass_flux(i,j) = div_m
    end do
end do

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
T(1,:,:) = T_surf_init(:,:) !SURFACE = 1
T(2,:,:) = T_deep_init(:,:) !DEEP = 2

deallocate(T_surf_init,T_deep_init)

!**Solving using the Forward Euler method**!

n_times = TIME_STEPS  !Number of time steps, 20 steps = 10 days.

!CONSTANTS-that don't depend on coordinates!
call calculate_F_a(F_net_sw_down,F_lw_down,F_latent_up,F_sensible_up,F_a)
d_theta = deg_to_rad(DELTA_LAT)
d_phi = deg_to_rad(DELTA_LON)
diff_coef = D
s_dfsn = DELTA_T*diff_coef/(R_PLANET)**2


!!****By this point the only arrays are: T, M, VERT_mass_flux, sv_flow_PHI ,sv_flow_Theta, F_a,
!!F_c ,F_net_sw_down,F_lw_down,F_latent_up and F_sensible_up****!!

!**Version Options hard code or ask via read input**!
!print*,'Enter version of ulmo (0,1 or 2):'
!read *,version
version = 2

allocate(T_new(2,N_LATS,N_LONS))

!** TIME STEPPING **!
do n_step = 1,n_times


    if(n_step==1) then
        print*,'Starting time stepping. This will take a while!'
    endif

    time = n_step*DELTA_T
    !print*,n_step


    do h = 1,N_DEPTHS
        do i = 1, N_LATS
            do j = 1,N_LONS



                call h_slab(h,thickness)
                s_notrns = DELTA_T/(RHO_WATER*C_V*thickness)
                s_ekmn = DELTA_T/(R_PLANET*RHO_WATER*thickness)
                theta = deg_to_rad(lats(i,1))

                call calculate_F_c(T,F_c)
                call calc_Q_flux(T,upward_Q_flux,F_net_sw_down,F_lw_down,F_latent_up,F_sensible_up)

                if(land_mask(i,j) == 1) then

                    T_new(h,i,j) = T(h,i,j)

                else
                    !**No transport**!
                    if(version == 0) then

                        if(h==1) then

                            call calc_new_T_surf_notrns(T,i,j,h,s_notrns,F_c,F_a,T_new)

                        elseif(h==2) then

                            call calc_new_T_deep_notrns(T,i,j,h,s_notrns,F_c,T_new)

                        else
                            print*, 'Index h is out of range'

                        end if

                    !**Diffusion**!
                    elseif(version==1) then

                        if(h==1) then

                            call calc_new_T_surf_diff(T,i,j,h,d_theta,d_phi,s_notrns,s_dfsn,F_c,F_a,theta,T_new)


                        elseif(h==2) then

                            call calc_new_T_deep_diff(T,i,j,h,d_theta,d_phi,s_notrns,s_dfsn,F_c,theta,T_new)

                        else

                            print*, 'Index h is out of range'

                        end if

                    !**Diffusion and Ekman transport**!
                    elseif(version==2) then

                        if(h==1) then
                            call calc_new_T_surf_ekman(T,M,i,j,h,d_theta,d_phi,s_notrns,s_dfsn,s_ekmn,F_c,F_a,&
                                                        theta,T_new)

                        elseif(h==2) then
                            call calc_new_T_deep_ekman(T,M,i,j,h,d_theta,d_phi,s_notrns,s_dfsn,s_ekmn,F_c,&
                                                        theta,T_new)

                        else
                            print*, 'Index h is out of range'

                        end if

                    end if

                end if

            end do
        end do
    end do

    T(:,:,:) = T_new(:,:,:)

    !**Writing outputs to files,(outputs every 100 days)**!
    if (mod(time,TIME_OUTPUT_FREQ)<TIME_OUTPUT_TOL) then
        days = T_OFFSET+time/(HOURS_PER_DAY*MINUTES_PER_HOUR*SECONDS_PER_MINUTE)
        print*, 'Days passed = ', days
    end if

    if (mod(time,DATA_OUTPUT_FREQ)<TIME_OUTPUT_TOL) then
        call process_output(T,upward_Q_flux,time,sv_flow_PHI,sv_flow_THETA,Vert_mass_flux)
    end if

end do

deallocate(Vert_mass_flux)
deallocate(sv_flow_PHI,sv_flow_THETA)
deallocate(T_new)
deallocate(land_mask,lats)
deallocate(upward_Q_flux,F_a,F_c)
deallocate(F_net_sw_down,F_lw_down,F_latent_up,F_sensible_up)
deallocate(T,M)

end program ULMO
