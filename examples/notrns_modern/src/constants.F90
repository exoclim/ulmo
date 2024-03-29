!*********************************************************************************
! This Module contains all of the planetary constants and file names for the model
!*********************************************************************************
module Constants
use, intrinsic :: iso_fortran_env
implicit none
private

!********************
!Planetary constants
!********************
  real(real64), parameter,public :: C_D = 0.0013 ! drag coefficient, dimensionless
  real(real64), parameter,public :: C_V = 4000.0 ! specific heat capacity of water J/kg/K
  real(real64), parameter,public :: EMISSIVITY = 0.985 ! effective blackbody emissivity of ocean water
  real(real64), parameter,public  :: HEAT_TRANSFER_COEFFICIENT = 3000. ! heat transfer coefficient for water W/m2/K 50-3000
  real(real64), parameter,public  :: RHO_AIR =  1.22 ! air density kg/m3
  real(real64), parameter,public  :: RHO_WATER =  1027 ! km/m3
  real(real64), parameter,public  :: R_PLANET = 1.12384*6.371e6 ! radius of planet in m (as fraction of Earth radius)
  real(real64), parameter,public  :: EPSI = 0.00001 ! pretty good need to find appropriate value
  real(real64), parameter,public  :: OMEGA = 6.501e-6 ! angular frequency planet, about its central axis
!  real(real64), parameter,public  :: OMEGA = 7.272e-5 ! Earth value for omega
  real(real64), parameter,public  :: SIGMA = 5.67e-8 ! stefan boltzman constant
  real(real64), parameter,public  :: H_S = 50.0 ! thickness of surface layer in m
  real(real64), parameter,public  :: H_D = 150.0 ! thickness of deep layer below
  !real(real64), parameter,public  :: D = 25000.0 ! horizontal diffusion coefficient m2/s
  real(real64),parameter,public   :: D = 250.0 !needs to be this for Diffusion or Ekmn (stability)
  integer(int64), parameter,public  :: N_LATS = 90 ! number of latitude points
  integer(int64), parameter,public  :: N_LONS = 144 ! number of longitude points
  real(real64), parameter,public  :: LAT_MIN = -89.0
  real(real64), parameter,public  :: DELTA_LAT = 2.0
  real(real64), parameter,public  :: LON_MIN = 1.25
  real(real64), parameter,public  :: DELTA_LON = 2.5
  real(real64), parameter,public  :: HOURS_PER_DAY = 24.
  real(real64), parameter,public  :: MINUTES_PER_HOUR = 60.
  real(real64), parameter,public  :: SECONDS_PER_MINUTE = 60.
  real(real64), parameter,public  :: DELTA_T = 12.0*MINUTES_PER_HOUR*SECONDS_PER_MINUTE ! time step of 1 hours
  integer(int64), parameter,public  :: TIME_STEPS = 20000*2 ! for time step of 12 hours, 20,000 time steps needed for 10,000 day run

  !integer(int64), parameter,public  :: TIME_STEPS = 50
  real(real64), parameter,public  :: TIME_OUTPUT_FREQ = 25.0*HOURS_PER_DAY*MINUTES_PER_HOUR*SECONDS_PER_MINUTE ! print update frequency in seconds, but first term is number of days
  real(real64), parameter,public  :: DATA_OUTPUT_FREQ = 100.0*HOURS_PER_DAY*MINUTES_PER_HOUR*SECONDS_PER_MINUTE ! output frequency in seconds, but first term is number of days
  real(real64), parameter,public  :: TIME_OUTPUT_TOL = 1e-6 ! days
  real(real64), parameter,public  :: T_OFFSET = 0.0 ! IF you want to restart run, specify start time here, so it doesn't save over stuff, and change input temp files to output ones
  real(real64), parameter,public  :: TOL = 1e-6
  integer(int64), parameter,public  :: MAX_ITER = 100
  real(real64), parameter,public  :: pi = 4*atan(1.)
  integer(int64),parameter,public :: N_DEPTHS = 2 ! 1 deep and 1 surface layer
!***********************
! input data file names
!***********************
  character(len=*),parameter,public  :: LAND_MASK_DATA = "input_data/Earth_modern.dat"
  real(real64)            ,parameter,public  :: MAX_FILE_LINE_SIZE = 4000
  character(len=*),parameter,public  :: SW_FLUX_NET_DOWN_DATA = "input_data/ProCb/surface_net_downward_shortwave_flux.dat"
  character(len=*),parameter,public  :: LW_FLUX_DOWN_DATA = "input_data/ProCb/surface_downwelling_longwave_flux_in_air.dat"
!  character(len=*),parameter,public  :: LW_FLUX_DOWN_DATA = "input_data/ProCb/surface_net_downward_longwave_flux.dat"
  character(len=*),parameter,public  :: LATENT_UP_FLUX_DATA = "input_data/ProCb/surface_upward_latent_heat_flux.dat"
  character(len=*),parameter,public  :: SENSIBLE_UP_FLUX_DATA = "input_data/ProCb/surface_upward_sensible_heat_flux.dat"
!  character(len=*),parameter,public  :: INITIAL_SURFACE_TEMP_DATA = "output_data/ProCb/T_surf_1000_days.dat"
!  character(len=*),parameter,public  :: INITIAL_DEEP_TEMP_DATA "output_data/ProCb/T_deep_1000_days.dat"
  character(len=*),parameter,public  :: INITIAL_SURFACE_TEMP_DATA = "input_data/ProCb/surface_temperature.dat"
  character(len=*),parameter,public  :: INITIAL_DEEP_TEMP_DATA = "input_data/ProCb/surface_temperature.dat"
  character(len=*),parameter,public  :: U_WIND_DATA = "input_data/ProCb/x_wind.dat"
  character(len=*),parameter,public  :: V_WIND_DATA = "input_data/ProCb/y_wind.dat"
  character(len=*),parameter,public  :: X_STRESS_DATA = "input_data/ProCb/surface_downward_eastward_stress.dat"
  character(len=*),parameter,public  :: Y_STRESS_DATA = "input_data/ProCb/surface_downward_northward_stress.dat"
  character(len=*),parameter,public  :: LATS_FILE = "input_data/lats.dat"
  character(len=*),parameter,public  :: LONS_FILE = "input_data/lons.dat"
!******************
! output data files
!******************
  real(real64)            ,parameter,public  :: MAX_FNAME_CHAR = 100
  character(len=50),parameter,public  :: OUTPUT_UPWARD_Q_FLUX = "output_data/ProCb/upward_surface_Q_flux_"
  character(len=50),parameter,public  :: OUTPUT_SURFACE_TEMP_DATA = "output_data/ProCb/T_surf_"
  character(len=50),parameter,public  :: OUTPUT_DEEP_TEMP_DATA = "output_data/ProCb/T_deep_"
  character(len=50),parameter,public  :: OUTPUT_X_MASS_FLUX_DATA = "output_data/ProCb/eastward_flow_"
  character(len=50),parameter,public  :: OUTPUT_Y_MASS_FLUX_DATA = "output_data/ProCb/northward_flow_"
  character(len=50),parameter,public  :: OUTPUT_VERITCAL_FLUX_DATA = "output_data/ProCb/vertical_mass_flux_"
  character(len=50),parameter,public  :: OUTPUT_DATA_EXT = "_days.dat"
end module Constants
