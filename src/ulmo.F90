
program ULMO
    use MASS_FLUX
    use WRITE_READ_DATA
    use NAMELIST
    use MATRIX_CALC
implicit none
!*********************
!! ULMO MAIN SCRIPT !!
!*********************
real, dimension(N_LATS,N_LONS) :: T_surf, T_deep,mass_flux_THETA,mass_flux_PHI,sv_mass_flux_PHI,sv_mass_flux_THETA
real,dimension(2,N_LATS,N_LONS) :: T, M
real,dimension(N_LATS,1) :: b

! COMBINING T_SURF AND T_DEEP INTO ONE 3D ARRAY FOR EASIER ANALYSIS !
T_surf= read_file(INITIAL_SURFACE_TEMP_DATA,N_LATS,N_LONS)
T_deep = read_file(INITIAL_DEEP_TEMP_DATA,N_LATS,N_LONS)

T(1,:,:) = T_surf
T(2,:,:) = T_deep

! COMBINING M_THETA and M_PHI INTO ONE 3D ARRAY FOR EASIER ANALYSIS !
mass_flux_THETA = calculate_mass_flux_THETA()
mass_flux_PHI = calculate_mass_flux_PHI()
sv_mass_flux_THETA =calculate_flow_sv_THETA()
sv_mass_flux_PHI =calculate_flow_sv_PHI()

M(1,:,:) = mass_flux_THETA
M(2,:,:) = mass_flux_PHI

!!! Writing mass flux outputs to files !!
!call write_file('output_data/sv_mass_flux_Phi.dat',sv_mass_flux_Phi,N_LATS,N_LONS)
!call write_file('output_data/mass_flux_PHI.dat',mass_flux_PHI,N_LATS,N_LONS)
!call write_file('output_data/sv_mass_flux_Theta.dat',sv_mass_flux_THETA,N_LATS,N_LONS)
!call write_file('output_data/mass_flux_Theta.dat',mass_flux_THETA,N_LATS,N_LONS)

!!! TESTING THE TRASNPORT AND TIME STEPPER WITH NO DEEP OCEAN !!!
!print"(a,i5)", 'Running Program...'
!call temp_after_time_step(N_STEPS=TIME_STEPS) ! Time_steps value is specified in NAMELIST.f90 file
!print"(a,i5)", 'Finished!'

!!! Matrix calculation testing !!!
!b = calculate_vector_b(T,LAND_MASK_DATA)
!call write_file('output_data/b.txt',b,N_LATS,1)


end program ULMO




