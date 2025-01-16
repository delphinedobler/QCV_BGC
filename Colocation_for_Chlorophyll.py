#!/usr/bin/env python
# coding: utf-8

# # COLOCATION TOOL FOR IN-SITU and COPERNICUS GRIDDED PRODUCTS 

# In[ ]:


# ADAPTED TO ARGO/GLIDER NEEDS FOR CHL VALIDATION
# Author: D.Dobler (Euro-Argo ERIC)
# Date: 2024-12-20

# comments:

# Use of the copernicus marine service library

# First the subset method was tested
# Second the opendataset method is tested (was is referred to as "lazy load", even if I don't really like this term that does not tell what's behind: this is an index loading)


# <!--TABLE OF CONTENTS-->
# Contents:
# - [COLOCATION TOOL FOR IN-SITU and COPERNICUS GRIDDED PRODUCTS](#COLOCATION-TOOL-FOR-IN-SITU-and-COPERNICUS-GRIDDED-PRODUCTS)
#   - [I - Libraries imports and credentials handling](#I---Libraries-imports-and-credentials-handling)
#   - [II - Main functions](#II---Main-functions)
#     - [II.a - copernicus_marine subset function](#II.a---copernicus_marine-subset-function)
#     - [II.b - Argo data related functions - direct access](#II.b---Argo-data-related-functions---direct-access)
#     - [II.c Cerbere files related functions](#II.c-Cerbere-files-related-functions)
#   - [II.d - get all observations for one workflow](#II.d---get-all-observations-for-one-workflow)
#     - [II.e - Distance computation function](#II.e---Distance-computation-function)
#     - [II.f - In-situ observation grouping function](#II.f---In-situ-observation-grouping-function)
#   - [III - Colocation](#III---Colocation)
#     - [III.a - configuration selection](#III.a---configuration-selection)
#     - [III.b - IN-SITU data selection](#III.b---IN-SITU-data-selection)
#     - [III.b - Define needed datasets and variables for Chlorophyll-A](#III.b---Define-needed-datasets-and-variables-for-Chlorophyll-A)
#     - [III.c - spatial resolution and boundaries of the copernicus datasets](#III.c---spatial-resolution-and-boundaries-of-the-copernicus-datasets)
#     - [III.d - group extraction by geograpical criterion](#III.d---group-extraction-by-geograpical-criterion)

# ## I - Libraries imports and credentials handling

# In[5]:


import copernicusmarine
# Nota Bene: Copernicusmarine (both python and CLI) does not work when Ivanti is active. 
import requests
import xarray as xr # xarray methods can sometimes be quite long, not sure why. Changed for netCDF4 library
from netCDF4 import Dataset
import netCDF4 as nc
import pandas as pd
import numpy as np
import datetime as dt
import time
import os
import matplotlib.pyplot as plt
import traceback
import sys
import multiprocessing as mp
from multiprocessing import Pool
from dask.distributed import LocalCluster


# To know all the options from the service, uncomment the following line:
#?copernicusmarine
#?copernicusmarine.subset
#get_ipython().run_line_magic('pinfo', 'copernicusmarine.open_dataset')

# copernicusmarine.login()
# it saved credentials within:
# C:\Users\ddobler\.copernicusmarine\.copernicusmarine-credentials


# ## II - Main functions

# ### II.a - copernicus_marine subset function

def get_cms_data(did,var,lonm,lonp,latm,latp,datm,datp,zm,zp,outd,outf):
    copernicusmarine.subset(
      dataset_id=did,
      variables=var,
      minimum_longitude=lonm,
      maximum_longitude=lonp,
      minimum_latitude=latm,
      maximum_latitude=latp,
      start_datetime=datm,
      end_datetime=datp,
      minimum_depth=zm,
      maximum_depth=zp,
      output_filename = outf,
      output_directory = outd,
      force_download=True, # Important, otherwise a prompt asks for downloading confirmation.
      overwrite_output_data=True # important because if False (default value): when the output file already exists, it adds a (n) at the end. This can prevent from fetching the correct file name
    )


def get_workflow_dataset_and_var(workflow_name):

    if workflow_name == "chl":
        
        dataset_rrs='cmems_obs-oc_glo_bgc-reflectance_my_l3-multi-4km_P1D'
        rrs_var=['RRS412','RRS443','RRS490','RRS555','RRS670']
        dataset_chl='cmems_obs-oc_glo_bgc-plankton_my_l3-multi-4km_P1D'
        chl_var=['CHL']
        dataset_Kd='cmems_obs-oc_glo_bgc-transp_my_l3-multi-4km_P1D'
        Kd_var=['KD490']

        l_dataset=[dataset_chl,dataset_rrs,dataset_Kd]
        
        d_dataset_var={}
        d_dataset_var[dataset_chl]=chl_var
        d_dataset_var[dataset_rrs]=rrs_var
        d_dataset_var[dataset_Kd]=Kd_var

    return l_dataset,d_dataset_var
        

def get_resolution(workflow_name,method,clear_cache=False,cache_dir='cache_files',verbose=False):

    # intialisation (list of the datasets spatio-temporal features or stf)
    i_dataset_stf={}
    l_dataset_stf={}

    # workflow datasets and vars:
    l_dataset,d_dataset_var=get_workflow_dataset_and_var(workflow_name)

    # Test is a cache file with value is present
    if not os.path.exists(cache_dir):
        os.mkdir(cache_dir)
    cache_resolution_file = cache_dir + "/cache_datasets_for_"+workflow_name+"_workflow_spatial_resolution.txt"

    if (os.path.exists(cache_resolution_file)) & (clear_cache):
        os.remove(cache_resolution_file)
        if verbose:print("the cache file was cleared")
    
    if not os.path.exists(cache_resolution_file):

        # Initialise the cache file
        file = open(cache_resolution_file, 'w')
        line2write = "dataset_id;reso_lon_deg;reso_lat_deg;lon_min;lon_max;lat_min;lat_max;reso_time_ns;time_min;time_max"
        file.write(line2write + '\n')
        file.close()
        
        for idataset in l_dataset:
            try:
                #if method=='lazy':
                if verbose:print("\n\n\n ******* Reading spatio-temporal features for " + idataset)
                ds=copernicusmarine.open_dataset(dataset_id=idataset)#,dataset_part='default',service='arco-geo-series')
                #print(ds)
                #if method=='subset':
                #    get_cms_data(idataset,d_dataset_var[idataset],
                #             0,0,0,0,"2022-06-01T00:00:00","2022-06-01T00:00:00",0,0,"copernicus-data",
                #             idataset+"_reso.nc")
                #    ds=xr.open_dataset("copernicus-data/"+idataset+"_reso.nc")
                i_dataset_stf['reso_lon_deg']=ds.attrs['lon_step']
                i_dataset_stf['reso_lat_deg']=ds.attrs['lat_step']
                i_dataset_stf['spat_lon_min']=ds.attrs['geospatial_lon_min']
                i_dataset_stf['spat_lon_max']=ds.attrs['geospatial_lon_max']
                i_dataset_stf['spat_lat_min']=ds.attrs['geospatial_lat_min']
                i_dataset_stf['spat_lat_max']=ds.attrs['geospatial_lat_max']
                # The global acceptable min and max in time assumes to use open_dataset ("lazy" method). There
                # is no global attribute helping in retrieving this boundary when subsetting.
                # The global attributes for time min and max (time_coverage_start and time_coverage_end) 
                # are wrongly filled: we must compute min and max from the time coordinate.
                # the time_coverage_resolution global attribute is a string, tricky to parse; here again, it is
                # better to recompute from the time coordinate
                #i_dataset_stf['reso_time']=ds.attrs['time_coverage_resolution']
                #i_dataset_stf['time_min']=ds.attrs['time_coverage_start']
                #i_dataset_stf['time_max']=ds.attrs['time_coverage_end']
                i_dataset_stf['time_min']=np.array(ds['time'].min())
                i_dataset_stf['time_max']=np.array(ds['time'].max())
                i_dataset_stf['reso_time_ns']=np.timedelta64(np.mean(np.array(ds['time'][1:])-np.array(ds['time'][:-1])),'ns')
            
                l_dataset_stf[idataset]=i_dataset_stf
                ds.close()
            except:
                print("ERROR: while downloading data from cmems with the " + method + " method")
                print(traceback.format_exc())

            try:
                file = open(cache_resolution_file, 'a')
                reso_time_str=str(np.timedelta64(i_dataset_stf['reso_time_ns'],'ns')/np.timedelta64(1,'ns'))
                line2write = idataset + ";" + str(i_dataset_stf['reso_lon_deg']) + ";" + str(i_dataset_stf['reso_lat_deg']) + ";" + \
                             str(i_dataset_stf['spat_lon_min']) + ";" + str(i_dataset_stf['spat_lon_max']) + ";" + \
                             str(i_dataset_stf['spat_lat_min']) + ";" + str(i_dataset_stf['spat_lat_max']) + ";" + \
                             reso_time_str + ";" + str(i_dataset_stf['time_min']) + ";" + str(i_dataset_stf['time_max'])
                file.write(line2write + '\n')
                
                file.close()
            except:
                print("ERROR: while writing the cache_resolution file: " + cache_resolution_file)
                file.close()
                if os.path.exists(cache_resolution_file):
                    os.remove(cache_resolution_file)
                print(traceback.format_exc())
    else:
        if verbose: print("Reading spatial resolution and boundaries from cache file")
        Reso_index=pd.read_csv(cache_resolution_file,sep=";")
        for idataset in l_dataset:
            i_dataset_stf['reso_lon_deg']=Reso_index['reso_lon_deg'][np.where(Reso_index['dataset_id']==idataset)[0][0]]
            i_dataset_stf['reso_lat_deg']=Reso_index['reso_lat_deg'][np.where(Reso_index['dataset_id']==idataset)[0][0]]
            i_dataset_stf['spat_lon_min']=Reso_index['lon_min'][np.where(Reso_index['dataset_id']==idataset)[0][0]]
            i_dataset_stf['spat_lon_max']=Reso_index['lon_max'][np.where(Reso_index['dataset_id']==idataset)[0][0]]
            i_dataset_stf['spat_lat_min']=Reso_index['lat_min'][np.where(Reso_index['dataset_id']==idataset)[0][0]]
            i_dataset_stf['spat_lat_max']=Reso_index['lat_max'][np.where(Reso_index['dataset_id']==idataset)[0][0]]
            reso_time_ns_int64=(Reso_index['reso_time_ns'][np.where(Reso_index['dataset_id']==idataset)[0][0]]).astype('int64')
            i_dataset_stf['reso_time_ns']=np.timedelta64(reso_time_ns_int64,'ns')
            i_dataset_stf['time_min']=np.datetime64(Reso_index['time_min'][np.where(Reso_index['dataset_id']==idataset)[0][0]])
            i_dataset_stf['time_max']=np.datetime64(Reso_index['time_max'][np.where(Reso_index['dataset_id']==idataset)[0][0]])
            l_dataset_stf[idataset]=i_dataset_stf

    if verbose:
        for idataset in l_dataset:
            print("\n For copernicus dataset: " + idataset)
            print("reso_longitude = ",i_dataset_stf['reso_lon_deg']," deg, reso_latitude = ",i_dataset_stf['reso_lat_deg'],"deg")
            print("spat_lon_min   = ",i_dataset_stf['spat_lon_min']," deg, spat_lon_max  = ",i_dataset_stf['spat_lon_max'],"deg")
            print("spat_lat_min   = ",i_dataset_stf['spat_lat_min']," deg, spat_lat_max  = ",i_dataset_stf['spat_lat_max'],"deg")
            print("reso_time_ns = ",i_dataset_stf['reso_time_ns'])
            print("time_min   = ",i_dataset_stf['time_min']," , time_max  = ",i_dataset_stf['time_max'])
            print(i_dataset_stf['reso_time_ns'].dtype)
            print(i_dataset_stf['time_min'].dtype)


    return l_dataset_stf


# ### II.b - Argo data related functions - direct access

def qc_from_char2int(val):
    shape_ini=val.shape
    tmp=((np.array(val)).astype('|S1')).tobytes().decode()
    tmp=tmp.replace(' ','6') # beware: this is only for computational issue, QC6 is unused usually
    tmp=list(tmp)
    out_val=np.array(tmp,dtype='int')
    out_val=np.reshape(out_val,shape_ini)
    return out_val


def get_file_from_url(URL,DL_FILE):
    response = requests.get(URL)
    if response.status_code == 404:
        print("No " + URL + " found in the gdac")
    else:
        open(DL_FILE, "wb").write(response.content)
        print(URL + " found in the gdac, locally downloaded in " + DL_FILE)


def get_and_read_Argo_meta_index():
    META_index_file="ar_index_global_meta.csv"
    URL = "https://data-argo.ifremer.fr/ar_index_global_meta.txt"
    get_file_from_url(URL,META_index_file)
    META_index=pd.read_csv(META_index_file,header=9,sep="/",names=['dac','wmo','remaining'],dtype={'dac': 'str', 'wmo': 'str', 'remaining': 'str'})
    return META_index


def get_dac_from_meta_index(wmo):
    META_index=get_and_read_Argo_meta_index()
    dac=META_index['dac'][np.where(META_index['wmo']==wmo)[0][0]]
    return dac


def get_argo_data_from_direct_access(wmo,workflow_name):
    SPROF_FILE=wmo + "_Sprof.nc"
    dac = get_dac_from_meta_index(wmo)
    URL = "https://data-argo.ifremer.fr/dac/"+dac+"/"+wmo+"/"+wmo+"_Sprof.nc"
    get_file_from_url(URL,SPROF_FILE)

    # xarray was sometimes taking several seconds for an unknown reason
    # As there is no challenge here in terms loading capacity
    # the NetCDF.Dataset was used: it never showed the delay experienced with xarray, thus kept
    ds=Dataset(SPROF_FILE,'r')
    ds.set_auto_mask(False) # to avoid masked array, a little bit more tricky to manage
    longitudes=ds.variables['LONGITUDE'][:]
    latitudes=ds.variables['LATITUDE'][:]
    position_qc=qc_from_char2int(ds.variables['POSITION_QC'][:])
    JULD=ds.variables['JULD'][:]
    dates_qc=qc_from_char2int(ds.variables['JULD_QC'][:])
    cycles=ds.variables['CYCLE_NUMBER'][:]
    PRES=ds.variables['PRES'][:]
    if workflow_name == 'chl':
        CHLA=ds.variables['CHLA'][:]
        CHLA_QC=ds.variables['CHLA_QC'][:]
        
    #print(ds)
    ds.close()
    
    ref_date=np.datetime64("1950-01-01T00:00:00")
    dates=JULD*86400*np.timedelta64(1, 's')+ref_date
    
    # output for colocation computation (1-D)
    df=pd.DataFrame({'CYCLE':cycles,'DATE': dates, 'LAT': latitudes, 'LON': longitudes, 'DATE_QC': dates_qc, 'POSITION_QC': position_qc})


    # output for colocation display
    n_prof,n_levels=PRES.shape
    prof=np.arange(n_prof)
    levels=np.arange(n_levels)
    
    ds = xr.Dataset(
        data_vars=dict(
            DATE=(["prof"], dates),
            LAT=(["prof"], latitudes),
            LON=(["prof"], longitudes),
            CYCLE=(["prof"], cycles),
            DATE_QC=(["prof"], dates_qc),
            POSITION_QC=(["prof"], position_qc),
            PRES=(["prof", "levels"], PRES),
            CHLA=(["prof", "levels"], CHLA),
            CHLA_QC=(["prof", "levels"], CHLA_QC),            
        ),
        coords=dict(
            prof=prof,
            n_levels=levels,
        ),
        attrs=dict(description="Observation related data"),
    )

    return df,ds


# ### II.c Cerbere files related functions


def get_argo_data_from_cerbere_access(cerbere_dir,wmo,workflow_name):
    
    cerbere_file = cerbere_dir + "gdac_" + wmo + "_202212_harm_agg.nc"
    ds=Dataset(cerbere_file,'r')
    ds.set_auto_mask(False) # to avoid masked array, a little bit more tricky to manage
    longitudes=ds.variables['lon'][:]
    latitudes=ds.variables['lat'][:]
    position_qc=qc_from_char2int(ds.variables['pos_qc'][:])
    JULD=ds.variables['time'][:]
    dates_qc=qc_from_char2int(ds.variables['time_qc'][:])
    cycles=ds.variables['cycle_number'][:]
    PRES=ds.variables['pressure_raw'][:]
    if workflow_name == 'chl':
        CHLA=ds.variables['chlorophylle_raw'][:]
        CHLA_QC=ds.variables['chlorophylle_raw_qc'][:]

    print(ds)
    ds.close()
    
    ref_date=np.datetime64("1950-01-01T00:00:00")
    dates=JULD*np.timedelta64(1, 's')+ref_date # In cerbere format, time units are seconds from ref-date
    

    # output for colocation computation (1-D)
    df=pd.DataFrame({'CYCLE':cycles,'DATE': dates, 'LAT': latitudes, 'LON': longitudes, 'DATE_QC': dates_qc, 'POSITION_QC': position_qc,})

    
    # output for colocation display
    n_prof,n_levels=PRES.shape
    prof=np.arange(n_prof)
    levels=np.arange(n_levels)
    
    ds = xr.Dataset(
        data_vars=dict(
            DATE=(["prof"], dates),
            LAT=(["prof"], latitudes),
            LON=(["prof"], longitudes),
            CYCLE=(["prof"], cycles),
            DATE_QC=(["prof"], dates_qc),
            POSITION_QC=(["prof"], position_qc),
            PRES=(["prof", "levels"], PRES),
            CHLA=(["prof", "levels"], CHLA),
            CHLA_QC=(["prof", "levels"], CHLA_QC),            
        ),
        coords=dict(
            prof=prof,
            n_levels=levels,
        ),
        attrs=dict(description="Observation related data"),
    )

    return df,ds


# ## II.d - get all observations for one workflow


def get_argo_data_from_index(workflow_name):
    
    BIO_Index_file="argo_bio-profile_index.csv"
    URL = "https://data-argo.ifremer.fr/argo_bio-profile_index.txt"
    get_file_from_url(URL,BIO_Index_file)
    BIO_Index=pd.read_csv(BIO_Index_file,header=8,sep=",")

    # Removing lines with incomplete coordinates:
    #print("intial size : ",BIO_Index.shape)
    BIO_Index.drop(BIO_Index.index[BIO_Index['date'].isnull()],inplace=True)
    #print("after removing null dates : ",BIO_Index.shape)
    BIO_Index.drop(BIO_Index.index[BIO_Index['latitude'].isnull()],inplace=True)
    #print("after removing null latitudes : ",BIO_Index.shape)
    BIO_Index.drop(BIO_Index.index[BIO_Index['longitude'].isnull()],inplace=True)
    #print("after removing null longitudes : ",BIO_Index.shape)

    # Keeping lines including the worklow parameter:
    BIO_Index=BIO_Index[BIO_Index['parameters'].str.contains("CHLA")]
    #print("after selecting chla lines : ",BIO_Index.shape)
  
   
    latitudes=np.array(BIO_Index['latitude'])
    longitudes=np.array(BIO_Index['longitude'])
    dates=BIO_Index['date'].astype(str)
    dates_str_iso_8601=dates.str[:4]+"-"+dates.str[4:6]+"-"+dates.str[6:8]+"T"+dates.str[8:10]+":"+dates.str[10:12]+":"+dates.str[12:14]
    dates_dt64=np.array(dates_str_iso_8601,dtype='datetime64')
    #print(dates_dt64)
    
    
    dates_qc=np.ones(dates_dt64.shape)
    position_qc=np.ones(dates_dt64.shape)
    
    df=pd.DataFrame({'DATE': dates_dt64, 'LAT': latitudes, 'LON': longitudes, 'DATE_QC': dates_qc, 'POSITION_QC': position_qc})

    return df #,ds


# ### II.e - Distance computation function

def compute_earth_radius_elliptical(lat_deg):
    
    # This function returns the earth radius at a given latitude, assuming an
    # elliptical earth model.
    
    if type(lat_deg)==np.ndarray:
        lat_deg=lat_deg.astype('float64')
        
    a=6378137 # equatorial radius
    b=6356750 # polar radius
    e=np.sqrt(1-(b/a)**2)
    lat_rad=lat_deg*np.pi/180
    earth_radius_m=a*np.sqrt(1-e**2*(np.sin(lat_rad))**2)
    
    return earth_radius_m
    
def compute_distance(lonA=0,latA=0,lonB=1,latB=0,verbose=False):
    
    # force float64 for input data to deal with default
    # ndarray dtype which is float32
    # and in this case, the computation is done in float32 
    # which can lead to up to 8% relative
    # error a distance of 1/12 deg (8/9 km).

    lonA=np.array(lonA).astype('float64')
    latA=np.array(latA).astype('float64')
    lonB=np.array(lonB).astype('float64')    
    latB=np.array(latB).astype('float64')
    
    
    #then compute earth_radius median
    #Earth_radius=6376*10**3 # in [m]
    lat_med=0.5*(latA+latB)
    Earth_radius=compute_earth_radius_elliptical(lat_med)
    #print(lat_med,Earth_radius)
    
    # first, put them in radians
    lonA_rad=lonA*np.pi/180
    latA_rad=latA*np.pi/180
    lonB_rad=lonB*np.pi/180
    latB_rad=latB*np.pi/180
    #print(lonA_rad,latA_rad,lonB_rad,latB_rad)
    
    if ((len(lonA.shape)!=0) & (len(lonB.shape) !=0)):
        if (len(lonA) > len(lonB)): distance=np.zeros(lonA.shape)
        else: distance=np.zeros(lonB.shape)
    if ((len(lonA.shape)==0) & (len(lonB.shape) !=0)):
        distance=np.zeros(lonB.shape)
    if ((len(lonA.shape)!=0) & (len(lonB.shape) ==0)):
        distance=np.zeros(lonA.shape)
    if ((len(lonA.shape)==0) & (len(lonB.shape) ==0)): 
        distance=0.0
    
    eps=1e-13
    is_A_an_array=False
    is_B_an_array=False
    try: 
        if np.size(lonA) > 1:is_A_an_array=True
    except: 
        eps=eps
    try: 
        if np.size(lonB) > 1:is_B_an_array=True
    except: 
        eps=eps
    
    if verbose:
        print("is_A_an_array,is_B_an_array:")
        print(is_A_an_array,is_B_an_array)
    
    
    if ((is_A_an_array==True) & (is_B_an_array==True)):
        #check where there is equality:
        i_neq=np.where((abs(lonA-lonB)>eps) | (abs(latA-latB) > eps))

        # then compute distance in [m]
        distance[i_neq]=Earth_radius[i_neq]*np.arccos(np.sin(latA_rad[i_neq])*np.sin(latB_rad[i_neq]) + \
                                 np.cos(latA_rad[i_neq])*np.cos(latB_rad[i_neq])*np.cos(lonB_rad[i_neq]-lonA_rad[i_neq]))
    
        
    if ( (is_A_an_array==False) & (is_B_an_array==True)):
        #check where there is equality:
        i_neq=np.where((abs(lonA-lonB)>eps) | (abs(latA-latB) > eps))
        # then compute distance in [m]
        AA=np.sin(latA_rad)*np.sin(latB_rad[i_neq])
        BB=np.cos(latA_rad)*np.cos(latB_rad[i_neq])*np.cos(lonB_rad[i_neq]-lonA_rad)
        cos_val=AA+BB
        distance[i_neq]=Earth_radius[i_neq]*np.arccos(cos_val)
    
    
    if ((is_A_an_array==True) & (is_B_an_array==False)):
        #check where there is equality:
        i_neq=np.where((abs(lonA-lonB)>eps) | (abs(latA-latB) > eps))

        # then compute distance in [m]
        distance[i_neq]=Earth_radius[i_neq]*np.arccos(np.sin(latA_rad[i_neq])*np.sin(latB_rad) + \
                                 np.cos(latA_rad[i_neq])*np.cos(latB_rad)*np.cos(lonB_rad-lonA_rad[i_neq]))
    
    
    if ((is_A_an_array==False) & (is_B_an_array==False)):
        if (abs(lonA-lonB)>eps) | (abs(latA-latB) > eps):
            distance=Earth_radius*np.arccos(np.sin(latA_rad)*np.sin(latB_rad)+ 
                                 np.cos(latA_rad)*np.cos(latB_rad)*np.cos(lonB_rad-lonA_rad))
    
            
    
    
    return distance


# ### II.f - In-situ observation grouping function

def create_obs_groups(gp_crit,i_dataset,i_dataset_stf,df_in_situ,verbose=False):
    
    # create groups by spatio-temporal criterion (will be referred to as medium cube)

    # transform capacity criterion into physical values
    # degree are kept because the copernicus grids are regular in degrees. This is also the reason why below, 
    # distances computation are converted into equivalent degree at the equator.For the time resolution, it is a little bit trickier. 
    # the global attributes is a string 'P1D' ... 
    gp_max_x_deg=gp_crit['gp_max_x_n']*i_dataset_stf['reso_lon_deg']
    gp_max_y_deg=gp_crit['gp_max_y_n']*i_dataset_stf['reso_lat_deg']
    gp_max_t_ns=gp_crit['gp_max_t_n']*i_dataset_stf['reso_time_ns']
    
    
    if verbose:
        print("gp_max_x_deg = {0:.3f}, gp_max_y_deg = {1:.3f}, gp_max_t_ns = {2:d}".format(gp_max_x_deg,gp_max_y_deg,gp_max_t_ns))
    
    # cast gp_max_t_ns in timedelta64 type (no more need to cast with new)
    #gp_max_t_days_dt64=np.timedelta64(gp_max_t_days,'D')
    
    # first create a "fictive" observation id list:
    list_obs_id = np.arange(0,df_in_situ.shape[0])
    print("Initial number of observation:", len(list_obs_id))
    
    i_group=0
    
    # create a dicionnary with the indexes of the various groups
    group_of_obs={}

    # First discard observations that are too old for the copernicus dataset:
    print("Discarding observations that are too old")
    group_of_obs_too_old={}
    i_too_old=np.where( (df_in_situ['DATE'] < i_dataset_stf['time_min']) )[0]
    i_obs_group_n=list_obs_id[i_too_old]
    group_of_obs_too_old=i_obs_group_n
    list_obs_id=np.delete(list_obs_id,i_too_old)
    print("Left number of observation:", len(list_obs_id))

    # Second discard observations that are too recent for the copernicus dataset:
    print("Discarding observations that are too recent")
    group_of_obs_too_recent={}
    i_too_recent=np.where( (df_in_situ['DATE'] > i_dataset_stf['time_max']) )[0]
    i_obs_group_n=list_obs_id[i_too_recent]
    group_of_obs_too_recent=i_obs_group_n
    list_obs_id=np.delete(list_obs_id,i_too_recent)
    print("Left number of observation:", len(list_obs_id))

    earth_radius_eq=compute_earth_radius_elliptical(0)
    
    while (len(list_obs_id) > 0) & (i_group<=df_in_situ.shape[0]) :
    #while (len(list_obs_id) > 0) & (i_group<=600) :
    
        lon = df_in_situ['LON'][list_obs_id[0]]
        lat = df_in_situ['LAT'][list_obs_id[0]]
        dat = df_in_situ['DATE'][list_obs_id[0]]
    
        lon_obs_left=df_in_situ['LON'][list_obs_id]
        lat_obs_left=df_in_situ['LAT'][list_obs_id]
        dat_obs_left=df_in_situ['DATE'][list_obs_id]
    
        #compute equatorial equivalent distances
        dist_m_2_deg_at_equat=(180/(np.pi*earth_radius_eq))
        dist_lon_deg=compute_distance(lon,0,lon_obs_left,np.zeros(lon_obs_left.shape)) * dist_m_2_deg_at_equat
        dist_lat_deg=compute_distance(0,lat,np.zeros(lat_obs_left.shape),lat_obs_left) * dist_m_2_deg_at_equat
        dist_time_dt64=abs(dat-dat_obs_left)
    
        if verbose:
            print("First observation position: {0:}, {1:.3f}°N {2:.3f}°E".format(dat,lat,lon))
            print(dist_lon_deg[:5])
            print((lon-df_in_situ['LON'])[:5])
            print(dist_lat_deg[:5])
            print((lat-df_in_situ['LAT'])[:5])
            print(dist_time_dt64[:5])

        # Select close-by observations that are within the dataset time boundaries
        # This time boundary filter can also be done before calling the function.
        i_close_by=np.where((dist_lon_deg<=gp_max_x_deg) & \
                            (dist_lat_deg<=gp_max_y_deg) & \
                            (dist_time_dt64 <= gp_max_t_ns) )[0]
        i_obs_group_n=list_obs_id[i_close_by]
        group_of_obs[i_group]=i_obs_group_n
        list_obs_id=np.delete(list_obs_id,i_close_by)
        
        if verbose:
            print("i_close_by=\n", i_close_by)
            print("i_obs_group_n=\n",i_obs_group_n)
            print("list_obs_id after deleting")
            print(list_obs_id)
            print("np.min(dat_obs_left[i_obs_group_n]),np.max(dat_obs_left[i_obs_group_n])")
            print(np.min(dat_obs_left[i_obs_group_n]),np.max(dat_obs_left[i_obs_group_n]))
            print("np.min(dist_time_dt64[i_obs_group_n]),np.max(dist_time_dt64[i_obs_group_n])")
            print(np.min(dist_time_dt64[i_obs_group_n]),np.max(dist_time_dt64[i_obs_group_n]))
            #print("[dist_time_dt64[i_obs_group_n] dat_obs_left[i_obs_group_n]]")
            #print([dist_time_dt64[i_obs_group_n] dat_obs_left[i_obs_group_n]])
            print('\n')
            
        print("i_group={0:d};nb_elt_group={1:d};n_elt_left_to_group={2:d}".format(i_group+1,len(i_close_by),len(list_obs_id)))
            
        i_group = i_group + 1  
    if verbose: print(group_of_obs)
    return group_of_obs,group_of_obs_too_old,group_of_obs_too_recent


# ## II.g - get_copernicus_data_for_a_group_of_obs

def get_copernicus_data_for_a_group_of_obs(dataset_id,d_dataset_var,group_of_obs,df_in_situ,i_dataset_stf,delta_px,outfile_dir,cache_index_file,
                                           copernicus_method,indexation_method,record_format,log_file,analysis_date,location,gp_crit,i_obs_group,verbose=True):

    stime=time.time()
    spatial_extension_square_deg_max=0
    temporal_extension_days_max=0
    analysis_date=str(np.datetime64('today','D'))

    # extract the grouped observations from df_in_situ:
    idf=df_in_situ.iloc[group_of_obs]
    dates_qc=np.array(idf['DATE_QC'])
    position_qc=np.array(idf['POSITION_QC'])
    latitudes=np.array(idf['LAT'])
    longitudes=np.array(idf['LON'])
    dates=np.array(idf['DATE'])
    
    # compute medium box boundaries, accounting for 180 crossing:
    i_good_position=np.where(((position_qc==1) | (position_qc==2) | (position_qc==5) | (position_qc==8)) )
    i_good_dates=np.where(((dates_qc==1) | (dates_qc==2) | (dates_qc==5) | (dates_qc==8)) )
    
    delta_lon=i_dataset_stf['reso_lon_deg']*delta_px['x']
    delta_lat=i_dataset_stf['reso_lat_deg']*delta_px['y']

    # compute latitude boundaries
    bbox_lat_min=max(i_dataset_stf['spat_lat_min'],np.min(latitudes[i_good_position])-delta_lat)
    bbox_lat_max=min(i_dataset_stf['spat_lat_max'],np.max(latitudes[i_good_position])+delta_lat)
    dlat=bbox_lat_max-bbox_lat_min

    # compute longitude boundaries and check whether 180 was crossed
    # 2024/12/18 update: let's do that much simplier
    # Is the 180° crossed (assumption: the medium cube will never reach 180° large):
    if np.max(longitudes[i_good_position])-np.min(longitudes[i_good_position]) < 180:
        cross_180=0
        bbox_lon_min=np.min(longitudes[i_good_position])-delta_lon
        bbox_lon_max=np.max(longitudes[i_good_position])+delta_lon
    else:
        cross_180=1
        bbox_lon_min=np.max(longitudes[i_good_position])-delta_lon
        bbox_lon_max=np.min(longitudes[i_good_position])+delta_lon
    dlon=bbox_lon_max-bbox_lon_min
        
    # compute time boundaries
    bbox_dates_min=str(np.datetime_as_string(np.min(dates[i_good_dates])-np.timedelta64(delta_px['t'],'D'),'s'))
    bbox_dates_max=str(np.datetime_as_string(np.max(dates[i_good_dates])+np.timedelta64(delta_px['t'],'D'),'s'))

    # compute spatio-temporal extensions:
    spatial_extenstion_square_deg=dlon*dlat
    temporal_extension_days=np.timedelta64(np.max(dates[i_good_dates])-np.min(dates[i_good_dates]),'D') / np.timedelta64(1, 'D')
    spatial_extenstion_square_deg_max=max(spatial_extenstion_square_deg,spatial_extension_square_deg_max)
    temporal_extension_days_max=max(temporal_extension_days_max,temporal_extension_days)

    # some debug printing:
    if verbose:
        print("obs_lat_min     = {0:.2f}\t\t, obs_lat_max     = {1:.2f}".format(np.min(latitudes[i_good_position]),np.max(latitudes[i_good_position])))
        print("bbox_lat_min    = {0:.2f}\t\t, bbox_lat_max    = {1:.2f}\t dlat = {2:.2f}".format(bbox_lat_min,bbox_lat_max,dlat))
        print("obs_lon_min     = {0:.2f}\t\t, obs_lon_max     = {1:.2f}".format(np.min(longitudes[i_good_position]),np.max(longitudes[i_good_position])))
        print("bbox_lon_min    = {0:.2f}\t\t, bbox_lon_max    = {1:.2f}\t dlon = {2:.2f}".format(bbox_lon_min,bbox_lon_max,dlon))
        print("bbox_dat_min = ",bbox_dates_min,"\t, bbox_dates_max = ",bbox_dates_max)
        print("spatial_extension = {0:.2f} square_degrees temporal_extension = {1:.1f} days".format(spatial_extenstion_square_deg,temporal_extension_days))

    
    # Define cache file name
    outfile_name="{0:s}_{1:s}_{2:s}_{3:.1f}_{4:.1f}_{5:.1f}_{6:.1f}.nc".format(dataset_id,bbox_dates_min[:10],bbox_dates_max[:10],bbox_lat_min,
                                                                        bbox_lat_max,bbox_lon_min,bbox_lon_max)

    if verbose: print("outfile_name=",outfile_name)

    # Test whether the information are not yet downloaded:
    # TODO: move that part when the group of observations are built
    cache_index=pd.read_csv(cache_index_file,sep=";",dtype={'dataset_id': 'str', 'date_min': 'str', 'date_max': 'str',
                                                                'lat_min' : 'float','lat_max':'float','lon_min' : 'float','lon_max':'float',
                                                                'cross_180' : 'int','file_name':'str'})
    index_line_exist=False
    if verbose: print("Number of lines in the cache index:",cache_index['dataset_id'].size)
    test_presence=np.array([])
    if cache_index['dataset_id'].size > 0:
        # dataset_id;date_min;date_max;lat_min;lat_max;lon_min;lon_max;cross_180;file_name
        test_presence=np.where((cache_index['dataset_id'] == dataset_id) &
                               (cache_index['date_min']   <= bbox_dates_min) &
                               (cache_index['date_max']   >= bbox_dates_max) &
                               (cache_index['lat_min']    <= bbox_lat_min) &
                               (cache_index['lat_max']    >= bbox_lat_max) &
                               (cache_index['lon_min']    <= bbox_lon_min) &
                               (cache_index['lon_max']    >= bbox_lon_max) & 
                               (cache_index['cross_180']  == cross_180))[0]
        print("test_presence=",test_presence)
        if (test_presence.size > 0): 
            index_line_exist=True
            print("the information to download already exists in a file, no need to download again from copernicus")


    if not index_line_exist:

        if copernicus_method == 'lazy':
            print("Subsetting data with the 'lazy' load")
            ii_dat=np.where((dat_cop > np.datetime64(bbox_dates_min)) & (dat_cop < np.datetime64(bbox_dates_max)))
            ii_lat=np.where((lat_cop > bbox_lat_min) & (lat_cop < bbox_lat_max))
            if cross_180 == 0:
                ii_lon=np.where((lon_cop > bbox_lon_min) & (lon_cop < bbox_lon_max))
            if cross_180 == 1:
                ii_lon=np.where((lon_cop < bbox_lon_max) | (lon_cop > bbox_lon_min))

            # Bench mark the different ways of subsetting an xarray dataset:
            # the in-between solution:
            if indexation_method == 'sel':
                ds_cop_group=ds_cop[d_dataset_var[dataset_id]].sel(time=dat_cop[ii_dat],
                                                          latitude=lat_cop[ii_lat],
                                                          longitude=lon_cop[ii_lon], 
                                                          method="nearest")
            #print(ii_dat[0])

            # The longest (even if counter-intuitive)
            if indexation_method == 'isel':
                ds_cop_group=ds_cop[d_dataset_var[dataset_id]].isel(time=ii_dat[0],
                                                      latitude=ii_lat[0],
                                                      longitude=ii_lon[0])
            
            # The fastest is direct indexing but can be done only by parameter, not the entire dataset.
            if indexation_method == 'direct':
                ds_cop_group=ds_cop[d_dataset_var[dataset_id][0]][np.min(ii_dat[0]):np.max(ii_dat[0]),
                                                               np.min(ii_lat[0]):np.max(ii_lat[0]),
                                                               np.min(ii_lon[0]):np.max(ii_lon[0])]

            
            print("lazy indexing ended")
            if record_format == 'NetCDF': # 'values' or 'NetCDF'
                ds_cop_group.to_netcdf(outfile_dir+"/"+outfile_name)
                print("to_netcdf ended")
            if record_format == 'values': # 'values' or 'NetCDF'
                for iVAR in d_dataset_var[dataset_id]:
                    print("Get variable ",iVAR)
                    if indexation_method == 'direct':
                        tmp=ds_cop_group.values
                    else:
                        tmp=ds_cop_group[iVAR].values                            
                    print(tmp)
                    print("to_values ended")
            if record_format == 'computation':
                print("entering " + record_format + " record format")
                #print(ds_cop_group)
                ds_average=ds_cop_group.mean()
                print(ds_average)
                print("exiting " + record_format + " record format")
                
        

        if copernicus_method=='subset':
            print("Subsetting data with the subset method")
            if cross_180 == 1 :
                # There is a need to split the request in two:
                
                outfile_name_1=outfile_name + "_1.nc"
                xmin=i_dataset_stf['spat_lon_min']
                xmax=min(bbox_lon_max,i_dataset_stf['spat_lon_max'])
                if verbose:
                    print("180 was crossed, splitting the subset request")
                    print("First request arguments:")
                    print(dataset_id,d_dataset_var[dataset_id],
                         xmin,xmax,bbox_lat_min,bbox_lat_max,bbox_dates_min,bbox_dates_max,0,0,
                         outfile_dir,outfile_name_1)
                get_cms_data(dataset_id,d_dataset_var[dataset_id],
                     xmin,xmax,bbox_lat_min,bbox_lat_max,bbox_dates_min,bbox_dates_max,0,0,
                     outfile_dir,outfile_name_1)
                
                outfile_name_2=outfile_name + "_2.nc"
                xmin=max(bbox_lon_min,i_dataset_stf['spat_lon_min'])
                xmax=i_dataset_stf['spat_lon_max']
                if verbose:
                    print("second request arguments:")
                    print(dataset_id,d_dataset_var[dataset_id],
                         xmin,xmax,bbox_lat_min,bbox_lat_max,bbox_dates_min,bbox_dates_max,0,0,
                         outfile_dir,outfile_name_1)
                get_cms_data(dataset_id,d_dataset_var[dataset_id],
                     xmin,xmax,bbox_lat_min,bbox_lat_max,bbox_dates_min,bbox_dates_max,0,0,
                     outfile_dir,outfile_name_2)

                # merge results
                print("merging results into " + outfile_dir + "/" + outfile_name)
                ds = xr.open_mfdataset([outfile_dir + "/" + outfile_name_1,outfile_dir + "/" + outfile_name_2], concat_dim=['longitude'], combine= "nested")
                ds.to_netcdf(outfile_dir + "/" + outfile_name)
                ds.close()
            else:
                xmin=max(bbox_lon_min,i_dataset_stf['spat_lon_min'])
                xmax=min(bbox_lon_max,i_dataset_stf['spat_lon_max'])
                get_cms_data(dataset_id,d_dataset_var[dataset_id],
                     xmin,xmax,bbox_lat_min,bbox_lat_max,bbox_dates_min,bbox_dates_max,0,0,
                     outfile_dir,outfile_name)

        # Saving in cache the 
        file = open(cache_index_file, 'a')
        #"dataset_id;date_min;date_max;lat_min;lat_max;lon_min;lon_max;cross_180;file_name"
        line2write = "{0:s};{1:s};{2:s};{3:.6f};{4:.6f};{5:.6f};{6:.6f};{7:d};{8:s}".format(dataset_id,bbox_dates_min,bbox_dates_max,
                                                                                            bbox_lat_min,bbox_lat_max,bbox_lon_min,
                                                                                            bbox_lon_max,cross_180,outfile_name)
        file.write(line2write + '\n')
        file.close()

    file = open(log_file, 'a')
    if (test_presence.size == 0):
        
        if (copernicus_method == 'lazy') & (record_format == 'values') : 
            file_size=0
        else:
            file_size=os.path.getsize(outfile_dir+"/"+outfile_name)

        if copernicus_method == 'lazy':
            str_method=copernicus_method + "_" + indexation_method
        else:
            str_method=copernicus_method
        
        line2write_fmt="{0:s};{1:s};{2:d}_{3:d}_{4:d};{5:s};{6:s};{7:s};{8:.0f};{9:.5f};{10:.2f};{11:.0f};{12:.0f}"
        line2write=line2write_fmt.format(analysis_date,location,gp_crit['gp_max_x_n'],gp_crit['gp_max_y_n'],gp_crit['gp_max_t_n'],dataset_id,str_method,
                                         record_format,i_obs_group,time.time()-stime,spatial_extenstion_square_deg_max,temporal_extension_days_max,
                                         file_size)
        
        print(line2write)
        file.write(line2write + '\n')
    file.close()


# ## III - Colocation

# ### III.a - configuration selection

if __name__ == '__main__':


    # TODO: use a config file for the following parameters
    # choose the input depending on your needs (the output can be tuned
    access_type='ARGO_INDEX' # 'ARGO_DIRECT' or 'ARGO_CERBERE' or 'ARGO_INDEX' for the moment. this parameter will be used afterwards for plugging cerberized data
    cerbere_dir="C:/Users/ddobler/Documents/08_DD_scripts/09_FAIR-EASE/cerbere-data/"

    #wmo='6901578' # long journey float
    wmo='6903024' # crosses 180 line (cycles 139 to 145 are on the West side of the line, the others on the East side)
    workflow_name='chl'

    clear_cache=True
    copernicus_method='subset' # 'lazy' or 'subset' : I kept both, can be tuned
    indexation_method='sel' # 'sel' or 'isel' or 'index' (in case of lazy access)
    record_format='NetCDF' # 'values' or 'NetCDF' or 'computation': either data are get (.values) or locally saved in a NetCDF file. Used for performance assessments
    verbose=False # the copernicus library can not yet be turned into quiet mode (but this works for informative prints)
    extract_data=True 

    # Depending on your capacity, tune the grouping options
    gp_max_x_n=25#15#25#50#100 # i.e. within gp_max_x_n*reso_lon_deg, e.g. 200*0.04 = 8 deg
    gp_max_y_n=25#15#25#50#100
    gp_max_t_n=50#30#50#100#200



    if copernicus_method == 'subset':
        record_format="NetCDF"
        indexation_method=""

    gp_crit={}
    gp_crit['gp_max_x_n']=gp_max_x_n
    gp_crit['gp_max_y_n']=gp_max_y_n
    gp_crit['gp_max_t_n']=gp_max_t_n

    if verbose:
        print("Estimate of the number of copernicus points to fetch: {:d}".format( gp_max_x_n*gp_max_y_n*gp_max_t_n))


    # ### III.b - IN-SITU data selection

    # In[ ]:


    if access_type == 'ARGO_DIRECT':
        df_in_situ,ds_in_situ=get_argo_data_from_direct_access(wmo,workflow_name)
    if access_type == 'ARGO_CERBERE':
        df_in_situ,ds_in_situ=get_argo_data_from_cerbere_access(cerbere_dir,wmo,workflow_name)
    if access_type == 'ARGO_INDEX':
        df_in_situ=get_argo_data_from_index(workflow_name)


    # ### III.b - Define needed datasets and variables for Chlorophyll-A

    # In[ ]:


    # Retrieve the copernicus dataset names and variables associated to the workflow
    l_dataset,d_dataset_var=get_workflow_dataset_and_var(workflow_name)


    # ### III.c - spatial resolution and boundaries of the copernicus datasets

    # Retrieve the spatial resolution and boundaries of the copernicus datasets
    stime=time.time()
    #l_dataset_stf=get_resolution(workflow_name,method=copernicus_method,clear_cache=clear_cache,verbose=verbose)
    l_dataset_stf=get_resolution(workflow_name,method=copernicus_method,clear_cache=False,verbose=verbose)
    print('Execution time: {0:.1f} s'.format(time.time()-stime))
    # Performance from Ifremer site: 17 s using subset method vs 10s in lazy load.


    # ### III.d - group extraction by geograpical criterion


    stime=time.time()
    group_of_obs,group_of_obs_too_old,group_of_obs_too_recent=create_obs_groups(gp_crit,l_dataset[0],l_dataset_stf[l_dataset[0]],df_in_situ[:2000],verbose=verbose)
    print('Execution time: {0:.1f} s'.format(time.time()-stime))


    # this will be a function
    # get_copernicus_data(wmo,dataset_id,cycle_step,delta_x_px,delta_y_px,delta_t_days)
    # print lines will be commented.
    delta_px={}
    delta_px['x']=5
    delta_px['y']=5
    delta_px['t']=5

    log_file='perfo.log'
    file = open(log_file, 'a')
    line2write="date;location;group_crit;dataset_id;copernicus_method;record_format;cycle_step;" +\
               "execution_time[s];spatial_extension[square_degrees];temporal_extension[days];cache file size[B]"
    print(line2write)
    file.write(line2write + '\n')
    file.close()

    analysis_date=str(np.datetime64('today','D'))
    location="office"

    # These lines initiate or read the cache file
    cache_dir="cache_files"
    if not os.path.exists(cache_dir):
        os.mkdir(cache_dir)
    cache_index_file = cache_dir + "/cache_dowloaded_data_index.txt"

    if (os.path.exists(cache_index_file)) & (clear_cache):
        os.remove(cache_index_file)
        if verbose:print("the cache file was cleared")

    if not os.path.exists(cache_index_file):
        file = open(cache_index_file, 'w')
        line2write = "dataset_id;date_min;date_max;lat_min;lat_max;lon_min;lon_max;cross_180;file_name"
        file.write(line2write + '\n')
        file.close()
        
    

    #for dataset_id in l_dataset:
    for dataset_id in [l_dataset[0]]:


        n_obs_group=len(group_of_obs)

        print("\n\n Workflow {0:s}; dataset {1:s} ".format(workflow_name,dataset_id))
        print("Variables to extract: ",d_dataset_var[dataset_id])

        # If method is lazy, pre-load index for the dataset
        if copernicus_method=='lazy':
            ds_cop=copernicusmarine.open_dataset(dataset_id=dataset_id)
            lat_cop=ds_cop['latitude']
            lon_cop=ds_cop['longitude']
            dat_cop=ds_cop['time']

        outfile_dir="copernicus-data/worflow_{0:s}_xn_{1:03d}_yn_{2:03d}_tn_{3:03d}".format(workflow_name,gp_crit['gp_max_x_n'],
                                                                                                gp_crit['gp_max_x_n'],
                                                                                                gp_crit['gp_max_t_n'])
        if not os.path.exists(outfile_dir):
            os.mkdir(outfile_dir)
        
        start = time.perf_counter()
        
        
        n_gr = int(sys.argv[1])
        n_gr = n_obs_group
        parallelisation = sys.argv[2]
        
        if parallelisation == 'mpProcess':
            pr={}
        if parallelisation == 'mpAsync':
            pool = Pool()
            res={}
        if parallelisation == 'dask':
            client = LocalCluster().get_client()
            res={}
            
        #for i_obs_group in range(n_obs_group):
        for i_obs_group in range(n_gr):

            #ii=input('start?')

            #print("\n i_obs_group/n_obs_group = ",i_obs_group+1,"/",n_obs_group)
            
            if parallelisation == 'no':
                print("\n i_obs_group/n_obs_group = ",i_obs_group+1,"/",n_obs_group, " no parallelisation")
                get_copernicus_data_for_a_group_of_obs(dataset_id,d_dataset_var,group_of_obs[i_obs_group],df_in_situ,
                                           l_dataset_stf[dataset_id],delta_px,outfile_dir,cache_index_file,
                                           copernicus_method,indexation_method,record_format,log_file,analysis_date,
                                           location,gp_crit,i_obs_group,verbose=verbose)
            if parallelisation == 'mpProcess':
                print("\n i_obs_group/n_obs_group = ",i_obs_group+1,"/",n_obs_group, " mp parallelisation active - method process")
                pr[i_obs_group]=mp.Process(target=get_copernicus_data_for_a_group_of_obs,
                                       args=(dataset_id,d_dataset_var,group_of_obs[i_obs_group],df_in_situ,
                                       l_dataset_stf[dataset_id],delta_px,outfile_dir,cache_index_file,
                                       copernicus_method,indexation_method,record_format,log_file,analysis_date,
                                       location,gp_crit,i_obs_group,verbose))

                pr[i_obs_group].start()
                
            if parallelisation == 'mpAsync':
                print("\n i_obs_group/n_obs_group = ",i_obs_group+1,"/",n_obs_group, " mp parallelisation active - method async")
                res[i_obs_group] = pool.apply_async(get_copernicus_data_for_a_group_of_obs, [dataset_id,d_dataset_var,group_of_obs[i_obs_group],df_in_situ,
                                   l_dataset_stf[dataset_id],delta_px,outfile_dir,cache_index_file,
                                   copernicus_method,indexation_method,record_format,log_file,analysis_date,
                                   location,gp_crit,i_obs_group,verbose])
                
            if parallelisation == 'dask':
                print("\n i_obs_group/n_obs_group = ",i_obs_group+1,"/",n_obs_group, " dask parallelisation active")
                res=client.submit(get_copernicus_data_for_a_group_of_obs,dataset_id,d_dataset_var,group_of_obs[i_obs_group],df_in_situ,
                                   l_dataset_stf[dataset_id],delta_px,outfile_dir,cache_index_file,
                                   copernicus_method,indexation_method,record_format,log_file,analysis_date,
                                   location,gp_crit,i_obs_group,verbose)
        
        for i_obs_group in range(n_gr):
            
            if parallelisation == 'mpProcess':
                pr[i_obs_group].join()

            if parallelisation == 'mpAsync':
                ans = res[i_obs_group].get(timeout=60)
        
        #if parallelisation == 'dask':
        #    res = client.gather(res)
            
        finish = time.perf_counter()
        print(f'It took {finish-start:.3f} second(s) to finish')
            
