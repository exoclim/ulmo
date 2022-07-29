
program ULMO
    use MASS_FLUX
    use WRITE_READ_DATA
    use NAMELIST
    use HEIGHT_OF_SLAB
    use TIME_STEP
    use TRANSPORT
implicit none
!*******************
!!! CODE TESTING !!!
!*******************

! MASS FLUX TESTING !

!real*8,  dimension(N_LATS,N_LONS) :: mass_flux_THETA
!real*8,  dimension(N_LATS,N_LONS) :: mass_flux_PHI
!real*8,  dimension(N_LATS,N_LONS) :: sv_mass_flux_THETA
!real*8,  dimension(N_LATS,N_LONS) :: sv_mass_flux_PHI
!
!mass_flux_THETA = calculate_mass_flux_THETA()
!mass_flux_PHI = calculate_mass_flux_PHI()
!sv_mass_flux_THETA =calculate_flow_sv_THETA()
!sv_mass_flux_PHI =calculate_flow_sv_PHI()
!
!
!!! Writing mass flux outputs to files !!
!call write_file('sv_mass_flux_Phi.dat',sv_mass_flux_Phi,N_LATS,N_LONS)
!call write_file('mass_flux_PHI.dat',mass_flux_PHI,N_LATS,N_LONS)
!call write_file('sv_mass_flux_Theta.dat',sv_mass_flux_THETA,N_LATS,N_LONS)
!call write_file('mass_flux_Theta.dat',mass_flux_THETA,N_LATS,N_LONS)

! TESTING THE TRASNPORT AND TIME STEPPER WITH NO DEEP OCEAN !
print"(a,i5)", 'Running Program...'
call temp_after_time_step(N_STEPS=TIME_STEPS) ! Time_steps value is specified in NAMELIST.f90 file
print"(a,i5)", 'Finished!'
end program ULMO




