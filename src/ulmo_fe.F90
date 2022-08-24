program ULMO

    !**This version of Ulmo uses the Forward Euler method to solving differential equations**!

    use MASS_FLUX, only: calculate_surface_stress_PHI,calculate_surface_stress_THETA,calculate_flow_sv_THETA &
    ,calculate_mass_flux_THETA,calculate_mass_flux_PHI,calculate_flow_sv_PHI
    use READ_DATA, only: read_file_real,read_file_int,write_file
    use Constants
    use heat_fluxes, only: calc_Q_flux,calculate_F_a,calculate_F_c
    use HEIGHT_OF_SLAB
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

!**Land mask and coordinate variables**!
integer(int64), dimension(:,:), allocatable :: land_mask
real(real64), dimension(:,:),allocatable :: lats,lats_data

!**Heat flux variables**!
real(real64), dimension(:,:), allocatable :: F_net_sw_down, F_lw_down,F_latent_up,F_sensible_up
real(real64), dimension(:,:,:), allocatable :: T
real(real64), dimension(:,:), allocatable :: upward_Q_flux, F_a, F_c
real(real64), dimension(:,:), allocatable :: T_surf_init,T_deep_init

!**Forward Euler variables**!
integer(int64) :: n_times,n_step,h,i,j
real(real64) :: days,time


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
T(1,:,:) = T_surf_init(:,:) !SURFACE = 1
T(2,:,:) = T_deep_init(:,:) !DEEP = 2

deallocate(T_surf_init,T_deep_init)

!**Solving using the Forward Euler method**!

n_times = 40000  !Number of time steps, 20 steps = 10 days.

call calculate_F_a(F_net_sw_down,F_lw_down,F_latent_up,F_sensible_up,F_a) ! Only depends on um inputs

do n_step = 1,n_times

!    if(n_step==1) then
!        write(*,"(a12,i5)") 'Time Step = ', n_step
!    else
!        write(*,"(i17)") n_step
!    end if

    do h = 1,N_DEPTHS
        do i = 1, N_LATS
            do j = 1,N_LONS

                call calculate_F_c(T,F_c)
                call calc_Q_flux(T,upward_Q_flux,F_net_sw_down,F_lw_down,F_latent_up,F_sensible_up)

                if(land_mask(i,j) == 1) then

                    T(h,i,j) = T(h,i,j)

                else

                    if(h==1) then
                        T(h,i,j) = (DELTA_T/(RHO_WATER*C_V*h_slab(h)))*(F_a(i,j)+F_c(i,j)+SIGMA*EPSILON*(T(h,i,j))**4) + T(h,i,j)

                    elseif(h==2) then
                        T(h,i,j) = (DELTA_T/(RHO_WATER*C_V*h_slab(h)))*(-F_c(i,j))+T(h,i,j)

                    else
                        print*, 'Index h is out of range'

                    end if

                end if

            end do
        end do
    end do

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



deallocate(land_mask)
deallocate(upward_Q_flux,F_a,F_c)
deallocate(F_net_sw_down,F_lw_down,F_latent_up,F_sensible_up)
deallocate(T,M)

end program ULMO
