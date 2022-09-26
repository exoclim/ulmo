#!/usr/bin/env python
# coding: utf-8

# In[1]:


import iris
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import warnings
warnings.filterwarnings("ignore")
from statistics import mean


# In[2]:


input_data_path = r"C:\Users\oakle\Ocean\src\input_data\ProCb\hadGEM3"
os.chdir(input_data_path)


# In[3]:


#Loading in .nc files
data_1 = iris.load_cube('ua_Amon_HadGEM3-GC31-LL_piControl_r1i1p1f1_gn_225001-234912.nc')
data_2 = iris.load_cube('va_Amon_HadGEM3-GC31-LL_piControl_r1i1p1f1_gn_225001-234912.nc')
data_3 = iris.load('oakley.nc')


# In[4]:


#Surface Temperature
air_temp = data_3[0]
surface_temperature = air_temp[1199]
surface_temperature_data = surface_temperature.data
np.savetxt('surface_temperature.dat',surface_temperature_data,delimiter="\t",fmt='%f8')


# In[5]:


#Shortwave
surface_net_downward_shortwave_flux = data_3[6]- data_3[2]
surface_net_downward_shortwave_flux = surface_net_downward_shortwave_flux[1199]
surface_net_downward_shortwave_flux_data = surface_net_downward_shortwave_flux.data
np.savetxt('surface_net_downward_shortwave_flux.dat',surface_net_downward_shortwave_flux_data,delimiter="\t",fmt='%f8')


# In[6]:


#Longwave
surface_downwelling_longwave_flux_in_air = data_3[1]
surface_downwelling_longwave_flux_in_air = surface_downwelling_longwave_flux_in_air[1199].data
np.savetxt('surface_downwelling_longwave_flux_in_air.dat',surface_downwelling_longwave_flux_in_air ,delimiter="\t",fmt='%f8')


# In[7]:


#latent heat
surface_upward_latent_heat_flux = data_3[3]
surface_upward_latent_heat_flux = surface_upward_latent_heat_flux[1199].data
np.savetxt('surface_upward_latent_heat_flux.dat',surface_upward_latent_heat_flux,delimiter="\t",fmt='%f8')


# In[8]:


#sensible heat flux
surface_upward_sensible_heat_flux = data_3[4]
surface_upward_sensible_heat_flux = surface_upward_sensible_heat_flux[1199].data
np.savetxt('surface_upward_sensible_heat_flux.dat',surface_upward_sensible_heat_flux,delimiter="\t",fmt='%f8')


# In[9]:


#x and y wind
x_wind = data_1[0][0]
y_wind = data_2[0][0]
lats = data_3[0].coord('latitude').points
lons = data_3[0].coord('longitude').points


# In[10]:


# Interpolating x and y wind from (144,192) grid to (145,192) grid

STASH_wind = iris.fileformats.pp.STASH(model=1, section=0, item=30)
def make_cube_of_zeros(nlat, nlon):
    """Make an iris cube of zeros with coordinates compatible with the UM grid."""

    # Create a numpy array of zeros.
    land_mask_data = np.zeros((nlat, nlon), dtype=np.int32)

    # Convert coordinate points to dimensional coordinates.
    # Only need 'coord_system' kwarg when using .mask files not .nc files (can remove for .nc files)
    lonc = iris.coords.DimCoord(lons, units="degrees", standard_name="longitude",
                                coord_system=iris.coord_systems.GeogCS(iris.fileformats.pp.EARTH_RADIUS))
    latc = iris.coords.DimCoord(lats, units="degrees", standard_name="latitude",
        coord_system=iris.coord_systems.GeogCS(iris.fileformats.pp.EARTH_RADIUS))

    # Assemble a cube.
    cube = iris.cube.Cube(
        land_mask_data,
        dim_coords_and_dims=((latc, 0), (lonc, 1)),
        standard_name="eastward_wind",
        units="1",
        attributes=dict(STASH=STASH_wind),
    )
    return cube

sample_points = make_cube_of_zeros(144,192)

x_wind.coord(axis='x').coord_system = sample_points.coord(axis='x').coord_system # removes coordinate mismatch
x_wind.coord(axis='y').coord_system = sample_points.coord(axis='y').coord_system # removes coordinate mismatch
y_wind.coord(axis='x').coord_system = sample_points.coord(axis='x').coord_system # removes coordinate mismatch
y_wind.coord(axis='y').coord_system = sample_points.coord(axis='y').coord_system # removes coordinate mismatch

#regridding land mask onto um grid
x_wind_new= x_wind.regrid(sample_points,iris.analysis.Linear())
iris.save(x_wind_new, "x_wind_new.nc")
y_wind_new= y_wind.regrid(sample_points,iris.analysis.Linear())
iris.save(y_wind_new, "y_wind_new.nc")


# In[11]:


# Since the interpolated data will have large numbers at continent boundaries due to values being filled with 1e20 at points 
# including Orography for original data now need to remove these and replace them with thier nearest neighbour values

#replacing 1e20 values with zero so only large numbers at continent boundaries remain.
data_x_wind = np.where(x_wind_new.data.data==1e20,0,x_wind_new.data.data)


# In[12]:


def calculate_new_lon(lon) :
    ''' At maximum or minimum latitude, the adjacent latitude point above or bellow, respectively, is located at the same latitude,
        but at the new longitudinal point (from lon to new_lon) '''
    N_LONS = 144
    new_lon = int(lon+N_LONS/2)
    if (new_lon >= N_LONS):
        new_lon = int(new_lon - N_LONS)
    return new_lon
    


# In[13]:


#replacing large numbers at continent boundaries with nearest neighbour values
for i in range(0,143):
    for j in range(0,191):
        
        if i == 0:
            if j == 0:
                if data_x_wind[i,j] > 1e8:
                    #looking for neighbouring longitude points
                    if data_x_wind[i,j+1] !=0 and data_x_wind[i,j+1] < 1e8:
                        data_x_wind[i,j] = data_x_wind[i,j+1]
                    elif data_x_wind[i,191] !=0 and data_x_wind[i,191] <1e8:
                        data_x_wind[i,j] = data_x_wind[i,191]
                    elif data_x_wind[i,j+2] !=0 and data_x_wind[i,j+2] <1e8:
                        data_x_wind[i,j] = data_x_wind[i,j+2]
                    elif data_x_wind[i,190] !=0 and data_x_wind[i,190] <1e8:
                        data_x_wind[i,j] = data_x_wind[i,190]
                    elif data_x_wind[i,j+3] !=0 and data_x_wind[i,j+3] <1e8:
                        data_x_wind[i,j] = data_x_wind[i,j+3]
                    elif data_x_wind[i,189] !=0 and data_x_wind[i,189] <1e8:
                        data_x_wind[i,j] = data_x_wind[i,189]
                    #looking for neighbouring lattitude points
                    elif data_x_wind[i+1,j] !=0 and data_x_wind[i+1,j] < 1e8:
                        data_x_wind[i,j] = data_x_wind[i+1,j]
                    elif data_x_wind[i,calculate_new_lon(j)] !=0 and data_x_wind[i,calculate_new_lon(j)] <1e8:
                        data_x_wind[i,j] = data_x_wind[calculate_new_lon(j),j]   

                    else:
                        print('need more conditions')
                else:
                    data_x_wind[i,j] = data_x_wind[i,j]

            elif j == 191:
                if data_x_wind[i,j] > 1e8:
                    #looking for neighbouring longitude points
                    if data_x_wind[i,0] !=0 and data_x_wind[i,0] < 1e8:
                        data_x_wind[i,j] = data_x_wind[i,0]
                    elif data_x_wind[i,j-1] !=0 and data_x_wind[i,j-1] <1e8:
                        data_x_wind[i,j] = data_x_wind[i,j-1]
                    elif data_x_wind[i,1] !=0 and data_x_wind[i,1] <1e8:
                        data_x_wind[i,j] = data_x_wind[i,1]
                    elif data_x_wind[i,j-2] !=0 and data_x_wind[i,j-2] <1e8:
                        data_x_wind[i,j] = data_x_wind[i,j-2]
                    elif data_x_wind[i,2] !=0 and data_x_wind[i,2] <1e8:
                        data_x_wind[i,j] = data_x_wind[i,2]                    
                    elif data_x_wind[i,j-3] !=0 and data_x_wind[i,j-3] <1e8:
                        data_x_wind[i,j] = data_x_wind[i,j-3]
                    #looking for neighbouring lattitude points
                    elif data_x_wind[i+1,j] !=0 and data_x_wind[i+1,j] < 1e8:
                        data_x_wind[i,j] = data_x_wind[i+1,j]
                    elif data_x_wind[i,calculate_new_lon(j)] !=0 and data_x_wind[i,calculate_new_lon(j)] <1e8:
                        data_x_wind[i,j] = data_x_wind[calculate_new_lon(j),j]  

                    else:
                        print('need more conditions')
                else:
                    data_x_wind[i,j] = data_x_wind[i,j]

            else:
                if data_x_wind[i,j] > 1e8:
                    #looking for neighbouring longitude points
                    if data_x_wind[i,j+1] !=0 and data_x_wind[i,j+1] < 1e8:
                        data_x_wind[i,j] = data_x_wind[i,j+1]
                    elif data_x_wind[i,j-1] !=0 and data_x_wind[i,j-1] <1e8:
                        data_x_wind[i,j] = data_x_wind[i,j-1]
                    elif data_x_wind[i,j+2] !=0 and data_x_wind[i,j+2] <1e8:
                        data_x_wind[i,j] = data_x_wind[i,j+2]
                    elif data_x_wind[i,j+3] !=0 and data_x_wind[i,j+3] <1e8:
                        data_x_wind[i,j] = data_x_wind[i,j+3]
                    elif data_x_wind[i,j-2] !=0 and data_x_wind[i,j-2] <1e8:
                        data_x_wind[i,j] = data_x_wind[i,j-2]
                    elif data_x_wind[i,j-3] !=0 and data_x_wind[i,j-3] <1e8:
                        data_x_wind[i,j] = data_x_wind[i,j-3]
                    #looking for neighbouring lattitude points
                    elif data_x_wind[i+1,j] !=0 and data_x_wind[i+1,j] < 1e8:
                        data_x_wind[i,j] = data_x_wind[i+1,j]
                    elif data_x_wind[i,calculate_new_lon(j)] !=0 and data_x_wind[i,calculate_new_lon(j)] <1e8:
                        data_x_wind[i,j] = data_x_wind[calculate_new_lon(j),j]  
                    else:
                        print('need more conditions')

                else:
                    data_x_wind[i,j] = data_x_wind[i,j]
        
        elif i== 143:
            if j == 0:
                if data_x_wind[i,j] > 1e8:
                    #looking for neighbouring longitude points
                    if data_x_wind[i,j+1] !=0 and data_x_wind[i,j+1] < 1e8:
                        data_x_wind[i,j] = data_x_wind[i,j+1]
                    elif data_x_wind[i,191] !=0 and data_x_wind[i,191] <1e8:
                        data_x_wind[i,j] = data_x_wind[i,191]
                    elif data_x_wind[i,j+2] !=0 and data_x_wind[i,j+2] <1e8:
                        data_x_wind[i,j] = data_x_wind[i,j+2]
                    elif data_x_wind[i,190] !=0 and data_x_wind[i,190] <1e8:
                        data_x_wind[i,j] = data_x_wind[i,190]
                    elif data_x_wind[i,j+3] !=0 and data_x_wind[i,j+3] <1e8:
                        data_x_wind[i,j] = data_x_wind[i,j+3]
                    elif data_x_wind[i,189] !=0 and data_x_wind[i,189] <1e8:
                        data_x_wind[i,j] = data_x_wind[i,189]
                     #looking for neighbouring lattitude points
                    elif data_x_wind[i,calculate_new_lon(j)] !=0 and data_x_wind[i,calculate_new_lon(j)] < 1e8:
                        data_x_wind[i,j] = data_x_wind[i,calculate_new_lon(j)]
                    elif data_x_wind[i-1,j] !=0 and data_x_wind[i-1,j] <1e8:
                        data_x_wind[i,j] = data_x_wind[i-1,j]   

                    else:
                        print('need more conditions')
                else:
                    data_x_wind[i,j] = data_x_wind[i,j]

            elif j == 191:
                if data_x_wind[i,j] > 1e8:
                    #looking for neighbouring longitude points
                    if data_x_wind[i,0] !=0 and data_x_wind[i,0] < 1e8:
                        data_x_wind[i,j] = data_x_wind[i,0]
                    elif data_x_wind[i,j-1] !=0 and data_x_wind[i,j-1] <1e8:
                        data_x_wind[i,j] = data_x_wind[i,j-1]
                    elif data_x_wind[i,1] !=0 and data_x_wind[i,1] <1e8:
                        data_x_wind[i,j] = data_x_wind[i,1]
                    elif data_x_wind[i,j-2] !=0 and data_x_wind[i,j-2] <1e8:
                        data_x_wind[i,j] = data_x_wind[i,j-2]
                    elif data_x_wind[i,2] !=0 and data_x_wind[i,2] <1e8:
                        data_x_wind[i,j] = data_x_wind[i,2]                    
                    elif data_x_wind[i,j-3] !=0 and data_x_wind[i,j-3] <1e8:
                        data_x_wind[i,j] = data_x_wind[i,j-3]
                     #looking for neighbouring lattitude points
                    elif data_x_wind[i,calculate_new_lon(j)] !=0 and data_x_wind[i,calculate_new_lon(j)] < 1e8:
                        data_x_wind[i,j] = data_x_wind[i,calculate_new_lon(j)]
                    elif data_x_wind[i-1,j] !=0 and data_x_wind[i-1,j] <1e8:
                        data_x_wind[i,j] = data_x_wind[i-1,j] 

                    else:
                        print('need more conditions')
                else:
                    data_x_wind[i,j] = data_x_wind[i,j]

            else:
                if data_x_wind[i,j] > 1e8:
                    #looking for neighbouring longitude points
                    if data_x_wind[i,j+1] !=0 and data_x_wind[i,j+1] < 1e8:
                        data_x_wind[i,j] = data_x_wind[i,j+1]
                    elif data_x_wind[i,j-1] !=0 and data_x_wind[i,j-1] <1e8:
                        data_x_wind[i,j] = data_x_wind[i,j-1]
                    elif data_x_wind[i,j+2] !=0 and data_x_wind[i,j+2] <1e8:
                        data_x_wind[i,j] = data_x_wind[i,j+2]
                    elif data_x_wind[i,j+3] !=0 and data_x_wind[i,j+3] <1e8:
                        data_x_wind[i,j] = data_x_wind[i,j+3]
                    elif data_x_wind[i,j-2] !=0 and data_x_wind[i,j-2] <1e8:
                        data_x_wind[i,j] = data_x_wind[i,j-2]
                    elif data_x_wind[i,j-3] !=0 and data_x_wind[i,j-3] <1e8:
                        data_x_wind[i,j] = data_x_wind[i,j-3]
                     #looking for neighbouring lattitude points
                    elif data_x_wind[i,calculate_new_lon(j)] !=0 and data_x_wind[i,calculate_new_lon(j)] < 1e8:
                        data_x_wind[i,j] = data_x_wind[i,calculate_new_lon(j)]
                    elif data_x_wind[i-1,j] !=0 and data_x_wind[i-1,j] <1e8:
                        data_x_wind[i,j] = data_x_wind[i-1,j] 
                    else:
                        print('need more conditions')

                else:
                    data_x_wind[i,j] = data_x_wind[i,j]



        else :
    
            if j == 0:
                if data_x_wind[i,j] > 1e8:
                    
                    #looking for neighbouring longitude points
                    if data_x_wind[i,j+1] !=0 and data_x_wind[i,j+1] < 1e8:
                        data_x_wind[i,j] = data_x_wind[i,j+1]
                    elif data_x_wind[i,191] !=0 and data_x_wind[i,191] <1e8:
                        data_x_wind[i,j] = data_x_wind[i,191]
                    elif data_x_wind[i,j+2] !=0 and data_x_wind[i,j+2] <1e8:
                        data_x_wind[i,j] = data_x_wind[i,j+2]
                    elif data_x_wind[i,190] !=0 and data_x_wind[i,190] <1e8:
                        data_x_wind[i,j] = data_x_wind[i,190]
                    elif data_x_wind[i,j+3] !=0 and data_x_wind[i,j+3] <1e8:
                        data_x_wind[i,j] = data_x_wind[i,j+3]
                    elif data_x_wind[i,189] !=0 and data_x_wind[i,189] <1e8:
                        data_x_wind[i,j] = data_x_wind[i,189]
                        
                    #looking for neighbouring lattitude points
                    elif data_x_wind[i+1,j] !=0 and data_x_wind[i+1,j] < 1e8:
                        data_x_wind[i,j] = data_x_wind[i+1,j]
                    elif data_x_wind[i-1,j] !=0 and data_x_wind[i-1,j] <1e8:
                        data_x_wind[i,j] = data_x_wind[i-1,j]                    

                    else:
                        print('need more conditions')
                else:
                    data_x_wind[i,j] = data_x_wind[i,j]

            elif j == 191:
                if data_x_wind[i,j] > 1e8:
                    #looking for neighbouring longitude points
                    if data_x_wind[i,0] !=0 and data_x_wind[i,0] < 1e8:
                        data_x_wind[i,j] = data_x_wind[i,0]
                    elif data_x_wind[i,j-1] !=0 and data_x_wind[i,j-1] <1e8:
                        data_x_wind[i,j] = data_x_wind[i,j-1]
                    elif data_x_wind[i,1] !=0 and data_x_wind[i,1] <1e8:
                        data_x_wind[i,j] = data_x_wind[i,1]
                    elif data_x_wind[i,j-2] !=0 and data_x_wind[i,j-2] <1e8:
                        data_x_wind[i,j] = data_x_wind[i,j-2]
                    elif data_x_wind[i,2] !=0 and data_x_wind[i,2] <1e8:
                        data_x_wind[i,j] = data_x_wind[i,2]                    
                    elif data_x_wind[i,j-3] !=0 and data_x_wind[i,j-3] <1e8:
                        data_x_wind[i,j] = data_x_wind[i,j-3]
                    #looking for neighbouring lattitude points
                    elif data_x_wind[i+1,j] !=0 and data_x_wind[i+1,j] < 1e8:
                        data_x_wind[i,j] = data_x_wind[i+1,j]
                    elif data_x_wind[i-1,j] !=0 and data_x_wind[i-1,j] <1e8:
                        data_x_wind[i,j] = data_x_wind[i-1,j]   

                    else:
                        print('need more conditions')
                else:
                    data_x_wind[i,j] = data_x_wind[i,j]

            else:
                if data_x_wind[i,j] > 1e8:
                    #looking for neighbouring longitude points
                    if data_x_wind[i,j+1] !=0 and data_x_wind[i,j+1] < 1e8:
                        data_x_wind[i,j] = data_x_wind[i,j+1]
                    elif data_x_wind[i,j-1] !=0 and data_x_wind[i,j-1] <1e8:
                        data_x_wind[i,j] = data_x_wind[i,j-1]
                    elif data_x_wind[i,j+2] !=0 and data_x_wind[i,j+2] <1e8:
                        data_x_wind[i,j] = data_x_wind[i,j+2]
                    elif data_x_wind[i,j+3] !=0 and data_x_wind[i,j+3] <1e8:
                        data_x_wind[i,j] = data_x_wind[i,j+3]
                    elif data_x_wind[i,j-2] !=0 and data_x_wind[i,j-2] <1e8:
                        data_x_wind[i,j] = data_x_wind[i,j-2]
                    elif data_x_wind[i,j-3] !=0 and data_x_wind[i,j-3] <1e8:
                        data_x_wind[i,j] = data_x_wind[i,j-3]
                    #looking for neighbouring lattitude points
                    elif data_x_wind[i+1,j] !=0 and data_x_wind[i+1,j] < 1e8:
                        data_x_wind[i,j] = data_x_wind[i+1,j]
                    elif data_x_wind[i-1,j] !=0 and data_x_wind[i-1,j] <1e8:
                        data_x_wind[i,j] = data_x_wind[i-1,j]   
                    else:
                        print('need more conditions')

                else:
                    data_x_wind[i,j] = data_x_wind[i,j]



# In[14]:


np.max(data_x_wind)


# In[15]:


truth = data_x_wind >= 1e8
for i in range(0,143):
    for j in range(0,192):
        if truth[i,j] == True:
            print('why?')
        else:
            pass
#for some reason there are 7 remaining points still large values. I am going to set these to zero


# In[19]:


for i in range(0,143):
    for j in range(0,192):
        if data_x_wind[i,j] >=1e8:
            data_x_wind[i,j] == 0
        else:
            data_x_wind[i,j] == data_x_wind[i,j]


# In[20]:


np.max(data_x_wind)


# In[18]:


### LITERALLY NO IDEA NOW ###


# In[ ]:




