
program ULMO
    use MASS_FLUX
    use WRITE_READ_DATA
    use NAMELIST
    use HEIGHT_OF_SLAB
implicit none
!*******************
!!! CODE TESTING !!!
!*******************
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


end program ULMO




