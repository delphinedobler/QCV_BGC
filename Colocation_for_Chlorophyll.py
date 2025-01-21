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
import portalocker
import pickle


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
        line2write = "dataset_id;reso_lon_deg;reso_lat_deg;lon_min;lon_max;lat_min;lat_max;reso_time_ns;time_min;time_max"
        file = open(cache_resolution_file, 'w')
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
                
                reso_time_str=str(np.timedelta64(i_dataset_stf['reso_time_ns'],'ns')/np.timedelta64(1,'ns'))
                line2write = idataset + ";" + str(i_dataset_stf['reso_lon_deg']) + ";" + str(i_dataset_stf['reso_lat_deg']) + ";" + \
                             str(i_dataset_stf['spat_lon_min']) + ";" + str(i_dataset_stf['spat_lon_max']) + ";" + \
                             str(i_dataset_stf['spat_lat_min']) + ";" + str(i_dataset_stf['spat_lat_max']) + ";" + \
                             reso_time_str + ";" + str(i_dataset_stf['time_min']) + ";" + str(i_dataset_stf['time_max'])
                             
                file = open(cache_resolution_file, 'a')
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


def get_argo_data_from_direct_access(argo_dir,wmo,workflow_name,dl=True):
    SPROF_FILE = argo_dir + "/" + wmo + "_Sprof.nc"
    
    if dl:
        
        dac = get_dac_from_meta_index(wmo)
        URL = "https://data-argo.ifremer.fr/dac/"+dac+"/"+wmo+"/"+wmo+"_Sprof.nc"
        print("Downloading in-situ data from ", URL)
        get_file_from_url(URL,SPROF_FILE)

    # xarray was sometimes taking several seconds for an unknown reason
    # As there is no challenge here in terms loading capacity
    # the NetCDF.Dataset was used: it never showed the delay experienced with xarray, thus kept
    print("Reading in-situ data from ", SPROF_FILE)
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


def get_argo_data_from_cerbere_access(cerbere_dir,wmo,workflow_name,dl=False):
    
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


def get_argo_data_from_index(argo_dir,workflow_name,dl=True):
    
    BIO_Index_file=argo_dir + "/argo_bio-profile_index.csv"
    
    if dl:
        URL = "https://data-argo.ifremer.fr/argo_bio-profile_index.txt"
        get_file_from_url(URL,BIO_Index_file)
        print("Downloading in-situ data from ", URL)
        
    
    print("Reading in-situ data from ", BIO_Index_file)
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


def get_bbox_from_df(df,i_dataset_stf,delta_px,verbose=False):
    #verbose=True
    spatial_extension_square_deg_max=0
    temporal_extension_days_max=0
    
    # extract the grouped observations from df:
    dates_qc=np.array(df['DATE_QC'])
    position_qc=np.array(df['POSITION_QC'])
    latitudes=np.array(df['LAT'])
    longitudes=np.array(df['LON'])
    dates=np.array(df['DATE'],dtype='datetime64')
    
    # compute boundaries
    if np.isscalar(df['DATE_QC']):
        latm=latitudes
        latp=latitudes
        lonm=longitudes
        lonp=longitudes
        datm=dates
        datp=dates
    else:
        i_good_position=np.where(((position_qc==1) | (position_qc==2) | (position_qc==5) | (position_qc==8)) )
        i_good_dates=np.where(((dates_qc==1) | (dates_qc==2) | (dates_qc==5) | (dates_qc==8)) )
        lat_ok=latitudes[i_good_position]
        lon_ok=longitudes[i_good_position]
        dat_ok=dates[i_good_dates]
        latm=np.min(lat_ok)
        latp=np.max(lat_ok)
        lonm=np.min(lon_ok)
        lonp=np.max(lon_ok)
        datm=np.min(dat_ok)
        datp=np.max(dat_ok)

    delta_lon=i_dataset_stf['reso_lon_deg']*delta_px['x']
    delta_lat=i_dataset_stf['reso_lat_deg']*delta_px['y']
        
    bbox_lat_min=max(i_dataset_stf['spat_lat_min'],latm-delta_lat)
    bbox_lat_max=min(i_dataset_stf['spat_lat_max'],latp+delta_lat)
    dlat=bbox_lat_max-bbox_lat_min

    # compute longitude boundaries and check whether 180 was crossed
    # 2024/12/18 update: let's do that much simplier
    # Is the 180° crossed (assumption: the medium cube will never reach 180° large):
    if lonp-lonm < 180:
        cross_180=0
        bbox_lon_min=lonm-delta_lon
        bbox_lon_max=lonp+delta_lon

        if bbox_lon_min < -180:
            cross_180 = 1
            bbox_lon_min = bbox_lon_min + 360

        if bbox_lon_max > 180:
            cross_180 = 1
            bbox_lon_max = bbox_lon_max - 360
            
    else:
        cross_180=1
        lonp=np.min(lon_ok[lon_ok>0])
        lonm=np.max(lon_ok[lon_ok<0])
        bbox_lon_min=lonp-delta_lon
        bbox_lon_max=lonm+delta_lon
    
    
    
    dlon=bbox_lon_max-bbox_lon_min
        
    # compute time boundaries
    bbox_dates_min=str(np.datetime_as_string(datm-np.timedelta64(delta_px['t'],'D'),'s'))
    bbox_dates_max=str(np.datetime_as_string(datp+np.timedelta64(delta_px['t'],'D'),'s'))

    # compute spatio-temporal extensions:
    spatial_extension_square_deg=dlon*dlat
    temporal_extension_days=np.timedelta64(datp-datm,'D') / np.timedelta64(1, 'D')
    spatial_extension_square_deg_max=max(spatial_extension_square_deg,spatial_extension_square_deg_max)
    temporal_extension_days_max=max(temporal_extension_days_max,temporal_extension_days)

    # some debug printing:
    if verbose:
        print("obs_lat_min     = {0:.2f}\t\t, obs_lat_max     = {1:.2f}".format(latm,latp))
        print("bbox_lat_min    = {0:.2f}\t\t, bbox_lat_max    = {1:.2f}\t dlat = {2:.2f}".format(bbox_lat_min,bbox_lat_max,dlat))
        print("obs_lon_min     = {0:.2f}\t\t, obs_lon_max     = {1:.2f}".format(lonm,lonp))
        print("bbox_lon_min    = {0:.2f}\t\t, bbox_lon_max    = {1:.2f}\t dlon = {2:.2f}".format(bbox_lon_min,bbox_lon_max,dlon))
        print("bbox_dat_min = ",bbox_dates_min,"\t, bbox_dates_max = ",bbox_dates_max)
        print("spatial_extension = {0:.2f} square_degrees temporal_extension = {1:.1f} days".format(spatial_extension_square_deg,temporal_extension_days))
        print(dates)
        print(latitudes)
        print(longitudes)
        
    bbox={}
    bbox['bbox_dates_min']=bbox_dates_min
    bbox['bbox_dates_max']=bbox_dates_max
    bbox['bbox_lat_min']=bbox_lat_min
    bbox['bbox_lat_max']=bbox_lat_max
    bbox['bbox_lon_min']=bbox_lon_min
    bbox['bbox_lon_max']=bbox_lon_max
    bbox['cross_180']=cross_180
    bbox['spatial_extension_square_deg']=spatial_extension_square_deg
    bbox['temporal_extension_days']=temporal_extension_days

    return bbox

def get_data_to_colocate(df,dataset_id,i_dataset_stf,delta_px,cache_copernicus_downloaded_data_index,verbose=False,log4debug=False,log_file_col_1="",log_file_col_2=""):
    
    # Test whether the information are not yet downloaded:
    cache_index=pd.read_csv(cache_copernicus_downloaded_data_index,sep=";",dtype={'dataset_id': 'str', 'date_min': 'str', 'date_max': 'str',
                                                                'lat_min' : 'float','lat_max':'float','lon_min' : 'float','lon_max':'float',
                                                                'cross_180' : 'int','file_name':'str'})

    
    # extract the observations from df:
    dates_qc=np.array(df['DATE_QC'])
    position_qc=np.array(df['POSITION_QC'])
    latitudes=np.array(df['LAT'])
    longitudes=np.array(df['LON'])
    dates=np.array(df['DATE'])
    
    if log4debug == True:
        file=open(log_file_col_1,'w')
        file.write("i_obs;date;lat;lon;date_qc;pos_qc\n")
        file.close()

        n_obs=np.size(dates_qc)
        for i_obs in range(n_obs):
            df_iobs=df.iloc[i_obs]
            line2write="{0:d};{1:s};{2:.4f};{3:.4f};{4:.0f};{5:.0f}".format(i_obs,str(df_iobs['DATE']),df_iobs['LAT'],df_iobs['LON'],df_iobs['DATE_QC'],df_iobs['POSITION_QC'])
            print(line2write)
            file=open(log_file_col_1,'a')
            file.write(line2write + '\n')
            file.close()
    
    # keep only correctly located observations
    i_good_obs=np.where(   ((position_qc==1) | (position_qc==2) | (position_qc==5) | (position_qc==8)) & \
                              ((dates_qc==1) | (dates_qc==2) | (dates_qc==5) | (dates_qc==8))  )[0]
    
    df_ok=df.iloc[i_good_obs]
    
    # test if the cache index contains lines
    if cache_index['dataset_id'].size > 0:
        
        if log4debug:
            file=open(log_file_col_2,'w')
            file.write("i_good_obs;date;lat;lon;bbox_datm;bbox_datp;bbox_latm;bbox_latp;bbox_lonm;bbox_lonp;cross_180;in_cache;cache_file_name\n")
            file.close()
        
        n_obs=np.size(i_good_obs)
        print("Searching data to colocate within ",n_obs," observations")
        
        i_obs_to_colocate=[]
        eps=1e-6 # must be consistent with cache precision
        for i_obs in range(n_obs):
            
            df_iobs=df_ok.iloc[i_obs]
            bbox=get_bbox_from_df(df_iobs,i_dataset_stf,delta_px)
 
            # dataset_id;date_min;date_max;lat_min;lat_max;lon_min;lon_max;cross_180;file_name

            test_presence=np.where((cache_index['dataset_id'] == dataset_id) &
                               (cache_index['date_min']   <= bbox['bbox_dates_min']) &
                               (cache_index['date_max']   >= bbox['bbox_dates_max']) &
                               (cache_index['lat_min']-eps    <= bbox['bbox_lat_min']) &
                               (cache_index['lat_max']+eps    >= bbox['bbox_lat_max']) &
                               ( ((cache_index['lon_min']-eps    <= bbox['bbox_lon_min']) &
                                  (cache_index['lon_max']+eps    >= bbox['bbox_lon_max']) ) |
                                 ((cache_index['lon_min']-eps    <= bbox['bbox_lon_min']) &
                                  (cache_index['lon_min']-eps    <= bbox['bbox_lon_max']) & 
                                  (cache_index['cross_180']  == 1)) |
                                 ((cache_index['lon_max']+eps    >= bbox['bbox_lon_min']) &
                                  (cache_index['lon_max']+eps    >= bbox['bbox_lon_max']) & 
                                  (cache_index['cross_180']  == 1)) )
                               )[0]

            
            #print(str(bbox['bbox_dates_min']))
            if test_presence.size == 0:
                in_cache="no"
                cache_file=""
            else:
                in_cache=str(test_presence[0])
                cache_file=cache_index['file_name'][test_presence[0]]
            
            line2write="{0:d};{1:s};{2:.4f};{3:.4f};{4:s};{5:s};{6:.4f};{7:.4f};{8:.4f};{9:.4f};{10:d};{11:s};{12:s}\n".format(i_obs,str(df_iobs['DATE']),df_iobs['LAT'],df_iobs['LON'],str(bbox['bbox_dates_min']),str(bbox['bbox_dates_max']),bbox['bbox_lat_min'],bbox['bbox_lat_max'],bbox['bbox_lon_min'],bbox['bbox_lon_max'],bbox['cross_180'],in_cache,cache_file)
            
            if log4debug:
                file=open(log_file_col_2,'a')
                file.write(line2write)
                file.close()
            
            if verbose & ~log4debug:print(line2write)
                                   
            #print("test_presence=",test_presence)
            #if (test_presence.size > 0): 
            #    print("i_obs ",i_obs," already colocated")
            if (test_presence.size == 0):
                if len(i_obs_to_colocate)==0:
                    i_obs_to_colocate=[i_obs]
                else:
                    i_obs_to_colocate.append(i_obs)
        
        print("Number of points to colocate:",len(i_obs_to_colocate))
        df_to_colocate=df_ok.iloc[i_obs_to_colocate]
        
    else:
        
        df_to_colocate=df_ok
        

    return df_to_colocate

def create_obs_groups(gp_crit,i_dataset,i_dataset_stf,df_in_situ_ini,verbose=False,log4debug=False,log_file_grp=""):
    
    # create groups by spatio-temporal criterion (will be referred to as medium cube)

    # transform capacity criterion into physical values
    # degree are kept because the copernicus grids are regular in degrees. This is also the reason why below, 
    # distances computation are converted into equivalent degree at the equator.For the time resolution, it is a little bit trickier. 
    # the global attributes is a string 'P1D' ... 
    gp_max_x_deg=gp_crit['gp_max_x_n']*i_dataset_stf['reso_lon_deg']
    gp_max_y_deg=gp_crit['gp_max_y_n']*i_dataset_stf['reso_lat_deg']
    gp_max_t_ns=gp_crit['gp_max_t_n']*i_dataset_stf['reso_time_ns']
    
    
    if verbose:
        print("gp_max_x_deg = {0:.3f}, gp_max_y_deg = {1:.3f}, gp_max_t_ns = {2:s}".format(gp_max_x_deg,gp_max_y_deg,gp_max_t_ns))
    
    # cast gp_max_t_ns in timedelta64 type (no more need to cast with new)
    #gp_max_t_days_dt64=np.timedelta64(gp_max_t_days,'D')
    
    # first create a "fictive" observation id list:
    list_obs_id = np.arange(0,df_in_situ_ini.shape[0])
    print("Initial number of observation:", len(list_obs_id))
    
    #initialise log file for investigations
    if log4debug:
        file=open(log_file_grp,'w')
        file.write("i_obs;date;lat;lon;date_qc;pos_qc;id_group;comment;nb_elt_group\n")
        file.close()
    
    i_group=0
    
    # create a dicionnary with the indexes of the various groups
    group_of_obs={}

    # First discard observations that are too old for the copernicus dataset:
    print("Discarding observations that are too old")
    group_of_obs_too_old={}
    i_too_old=np.where( (df_in_situ_ini['DATE'] < i_dataset_stf['time_min']) )[0]
    i_obs_group_n=list_obs_id[i_too_old]
    group_of_obs_too_old=i_obs_group_n
    list_obs_id=np.delete(list_obs_id,i_too_old)
    print("Number of too old obs = {0:d}; Remaining obs to group = {1:d}".format(len(i_obs_group_n),len(list_obs_id)))
    if log4debug:
        for i_obs_too_old in i_obs_group_n:
            line2write= "{0:d};{1:s};{2:.4f};{3:.4f};{4:.0f};{5:.0f};no_group;too_old\n".format(i_obs_too_old,str(df_in_situ_ini['DATE'][i_obs_too_old]),df_in_situ_ini['LAT'][i_obs_too_old], \
            df_in_situ_ini['LON'][i_obs_too_old],df_in_situ_ini['DATE_QC'][i_obs_too_old],df_in_situ_ini['POSITION_QC'][i_obs_too_old])
            file=open(log_file_grp,'a')
            file.write(line2write)
            file.close()
        

    # Second discard observations that are too recent for the copernicus dataset:
    print("Discarding observations that are too recent")
    group_of_obs_too_recent={}
    i_too_recent=np.where( (df_in_situ_ini['DATE'] > i_dataset_stf['time_max']) )[0]
    i_obs_group_n=list_obs_id[i_too_recent]
    group_of_obs_too_recent=i_obs_group_n
    list_obs_id=np.delete(list_obs_id,i_too_recent)
    print("Number of too recent obs = {0:d}; Remaining obs to group = {1:d}".format(len(i_obs_group_n),len(list_obs_id)))
    if log4debug:
        for i_obs_too_recent in i_obs_group_n:
            line2write= "{0:d};{1:s};{2:.4f};{3:.4f};{4:.0f};{5:.0f};no_group;too_recent\n".format(i_obs_too_recent,str(df_in_situ_ini['DATE'][i_obs_too_recent]),df_in_situ_ini['LAT'][i_obs_too_recent], \
            df_in_situ_ini['LON'][i_obs_too_recent],df_in_situ_ini['DATE_QC'][i_obs_too_recent],df_in_situ_ini['POSITION_QC'][i_obs_too_recent])
            file=open(log_file_grp,'a')
            file.write(line2write)
            file.close()

    earth_radius_eq=compute_earth_radius_elliptical(0)
    
    while (len(list_obs_id) > 0) & (i_group<=df_in_situ_ini.shape[0]) :
    #while (len(list_obs_id) > 0) & (i_group<=600) :
    
        lon = df_in_situ_ini['LON'][list_obs_id[0]]
        lat = df_in_situ_ini['LAT'][list_obs_id[0]]
        dat = df_in_situ_ini['DATE'][list_obs_id[0]]
    
        lon_obs_left=df_in_situ_ini['LON'][list_obs_id]
        lat_obs_left=df_in_situ_ini['LAT'][list_obs_id]
        dat_obs_left=df_in_situ_ini['DATE'][list_obs_id]
    
        #compute equatorial equivalent distances
        dist_m_2_deg_at_equat=(180/(np.pi*earth_radius_eq))
        dist_lon_deg=compute_distance(lon,0,lon_obs_left,np.zeros(lon_obs_left.shape)) * dist_m_2_deg_at_equat
        dist_lat_deg=compute_distance(0,lat,np.zeros(lat_obs_left.shape),lat_obs_left) * dist_m_2_deg_at_equat
        dist_time_dt64=abs(dat-dat_obs_left)
    
        if verbose:
            print("First observation position: {0:}, {1:.3f}°N {2:.3f}°E".format(dat,lat,lon))
            print(dist_lon_deg[:5])
            print((lon-df_in_situ_ini['LON'])[:5])
            print(dist_lat_deg[:5])
            print((lat-df_in_situ_ini['LAT'])[:5])
            print(dist_time_dt64[:5])

        # Select close-by observations that are within the dataset time boundaries
        # This time boundary filter can also be done before calling the function.
        i_close_by=np.where((dist_lon_deg<=gp_max_x_deg) & \
                            (dist_lat_deg<=gp_max_y_deg) & \
                            (dist_time_dt64 <= gp_max_t_ns) )[0]
        i_obs_group_n=list_obs_id[i_close_by]
        group_of_obs[i_group]=i_obs_group_n
        list_obs_id=np.delete(list_obs_id,i_close_by)
        
        if log4debug:
            for i_obs in i_obs_group_n:
                line2write= "{0:d};{1:s};{2:.4f};{3:.4f};{4:.0f};{5:.0f};{6:d};;{7:d}\n".format(i_obs,str(df_in_situ_ini['DATE'][i_obs]),df_in_situ_ini['LAT'][i_obs], \
                df_in_situ_ini['LON'][i_obs],df_in_situ_ini['DATE_QC'][i_obs],df_in_situ_ini['POSITION_QC'][i_obs],i_group+1,len(i_close_by))
                file=open(log_file_grp,'a')
                file.write(line2write)
                file.close()
        
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
            
        if(i_group%100==0):print("i_group={0:d};nb_elt_group={1:d};n_elt_left_to_group={2:d}".format(i_group+1,len(i_close_by),len(list_obs_id)))
            
        i_group = i_group + 1  
    if verbose: print(group_of_obs)
    return group_of_obs,group_of_obs_too_old,group_of_obs_too_recent


# ## II.g - get_copernicus_data_for_a_group_of_obs


def get_copernicus_data_for_a_group_of_obs(dataset_id,d_dataset_var,i_obs_group,group_of_obs,df_in_situ_ini,i_dataset_stf,delta_px,outfile_dir,cache_copernicus_downloaded_data_index,
                                           copernicus_method,indexation_method,record_format,log_file_cop,analysis_date,location,gp_crit,verbose=True):

    stime=time.time()
    
    analysis_date=str(np.datetime64('today','D'))

    if verbose:
        print("################################")
        print("Getting bbox for group number ",i_obs_group)
        print("It contains the following i_obs:")
        print(group_of_obs)
        
    bbox=get_bbox_from_df(df_in_situ_ini.iloc[group_of_obs],i_dataset_stf,delta_px,verbose=verbose)
    
    bbox_dates_min=bbox['bbox_dates_min']
    bbox_dates_max=bbox['bbox_dates_max']
    bbox_lat_min=bbox['bbox_lat_min']
    bbox_lat_max=bbox['bbox_lat_max']
    bbox_lon_min=bbox['bbox_lon_min']
    bbox_lon_max=bbox['bbox_lon_max']
    cross_180=bbox['cross_180']
    spatial_extension_square_deg=bbox['spatial_extension_square_deg']
    temporal_extension_days=bbox['temporal_extension_days']


    if copernicus_method == 'lazy':
        print("Subsetting data with the 'lazy' load")
        outfile_name="{0:s}_{1:s}_{2:s}_{3:.1f}_{4:.1f}_{5:.1f}_{6:.1f}.nc".format(dataset_id,bbox_dates_min[:10],bbox_dates_max[:10],bbox_lat_min,
                                                                        bbox_lat_max,bbox_lon_min,bbox_lon_max)
        print("outfile_name=",outfile_name)
        
        
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
        tmin=np.datetime_as_string(max(np.datetime64(bbox_dates_min),np.datetime64(i_dataset_stf['time_min'])),unit='s')
        tmax=np.datetime_as_string(min(np.datetime64(bbox_dates_max),np.datetime64(i_dataset_stf['time_max'])),unit='s')
        ymin=max(bbox_lat_min,i_dataset_stf['spat_lat_min'])
        ymax=min(bbox_lat_max,i_dataset_stf['spat_lat_max'])
        
        outfile_name="{0:s}_{1:s}_{2:s}_{3:.1f}_{4:.1f}_{5:.1f}_{6:.1f}.nc".format(dataset_id,tmin[:10],tmax[:10],ymin,ymax,bbox_lon_min,bbox_lon_max)
        print("outfile_name=",outfile_name)
        if cross_180 == 1 :
            # There is a need to split the request in two:
            
            outfile_name_1=outfile_name + "_1.nc"
            xmin=i_dataset_stf['spat_lon_min']
            xmax=min(bbox_lon_max,i_dataset_stf['spat_lon_max'])
            if verbose:
                print("180 was crossed, splitting the subset request")
                print("First request arguments:")
                print(dataset_id,d_dataset_var[dataset_id],
                     xmin,xmax,ymin,ymax,tmin,tmax,0,0,
                     outfile_dir,outfile_name_1)
            get_cms_data(dataset_id,d_dataset_var[dataset_id],xmin,xmax,ymin,ymax,tmin,tmax,0,0,outfile_dir,outfile_name_1)
            
            outfile_name_2=outfile_name + "_2.nc"
            xmin=max(bbox_lon_min,i_dataset_stf['spat_lon_min'])
            xmax=i_dataset_stf['spat_lon_max']
            if verbose:
                print("second request arguments:")
                print(dataset_id,d_dataset_var[dataset_id],
                     xmin,xmax,ymin,ymax,tmin,tmax,0,0,
                     outfile_dir,outfile_name_1)
            get_cms_data(dataset_id,d_dataset_var[dataset_id],xmin,xmax,ymin,ymax,tmin,tmax,0,0,outfile_dir,outfile_name_2)

            # merge results
            print("merging results into " + outfile_dir + "/" + outfile_name)
            ds = xr.open_mfdataset([outfile_dir + "/" + outfile_name_1,outfile_dir + "/" + outfile_name_2], concat_dim=['longitude'], combine= "nested")
            ds.to_netcdf(outfile_dir + "/" + outfile_name)
            ds.close()
        else:
            xmin=max(bbox_lon_min,i_dataset_stf['spat_lon_min'])
            xmax=min(bbox_lon_max,i_dataset_stf['spat_lon_max'])
            get_cms_data(dataset_id,d_dataset_var[dataset_id],xmin,xmax,ymin,ymax,tmin,tmax,0,0,outfile_dir,outfile_name)

    
    # Saving in the cache file 
    #"dataset_id;date_min;date_max;lat_min;lat_max;lon_min;lon_max;cross_180;file_name"
    line2write = "{0:s};{1:s};{2:s};{3:.6f};{4:.6f};{5:.6f};{6:.6f};{7:d};{8:s};{9:d}".format(dataset_id,tmin,tmax,ymin,ymax,bbox_lon_min,bbox_lon_max,cross_180,outfile_name,i_obs_group)
    locked=True
    nb_tries=0
    while (locked == True) and (nb_tries < 100):
        try:
            file = open(cache_copernicus_downloaded_data_index, 'a')
            portalocker.lock(file, portalocker.LockFlags.EXCLUSIVE)
            file.write(line2write + '\n')
            portalocker.unlock(file)
            file.close()
            locked=False
        except:
            locked=True
        nb_tries=nb_tries+1

    
    # logging information for performance and debug purpose
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
                                     record_format,i_obs_group,time.time()-stime,spatial_extension_square_deg,temporal_extension_days,
                                     file_size)
    
    print(line2write)
    locked=True
    nb_tries=0
    while (locked == True) and (nb_tries < 100):
        try:
            file = open(log_file_cop, 'a')
            portalocker.lock(file, portalocker.LockFlags.EXCLUSIVE)
            file.write(line2write + '\n')
            portalocker.unlock(file)
            file.close()
            locked=False
        except:
            locked=True
        nb_tries=nb_tries+1



# ## III - Colocation

# ### III.a - configuration selection

if __name__ == '__main__':

    
    import Colocation_cfg as cf
    
    access_type=cf.access_type
    cerbere_dir=cf.cerbere_dir
    argo_dir=cf.argo_dir
    wmo=cf.wmo
    workflow_name=cf.workflow_name
    clear_copernicus_resolution_cache=cf.clear_copernicus_resolution_cache
    clear_copernicus_downloaded_data_index_cache=cf.clear_copernicus_downloaded_data_index_cache
    copernicus_method=cf.copernicus_method
    indexation_method=cf.indexation_method
    record_format=cf.record_format
    verbose=cf.verbose
    gp_crit=cf.gp_crit
    delta_px=cf.delta_px
    cache_dir=cf.cache_dir
    cache_copernicus_downloaded_data_index=cf.cache_copernicus_downloaded_data_index
    log_file_cop=cf.log_file_cop
    log_file_col_1=cf.log_file_col_1
    log_file_col_2=cf.log_file_col_2
    log_file_grp=cf.log_file_grp
    location=cf.location
    outfile_dir=cf.outfile_dir

    if copernicus_method == 'subset':
        record_format="NetCDF"
        indexation_method=""

    if verbose:
        print("Estimate of the number of copernicus points to fetch: {:d}".format(gp_crit['gp_max_x_n']*gp_crit['gp_max_y_n']*gp_crit['gp_max_t_n']))
    
    # Debug parameterization
    # in steps to run: beware of choosing one of the possibility
    steps_2_run=[""]
    #steps_2_run=["get_copernicus_dataset_resolution"]
    #steps_2_run=["get_insitu_data","get_copernicus_dataset_resolution","compute_group_of_obs_from_in_situ_data"]
    #steps_2_run=["get_insitu_data","get_copernicus_dataset_resolution","get_remaining_data_2_colocate_from_cache","compute_group_of_obs_from_in_situ_data"]
    i0,i1=int(sys.argv[2]),int(sys.argv[3])
    parallelisation = sys.argv[1]

    # ### III.b - IN-SITU data selection
    print("\n#STEP 1: GET IN SITU DATA FROM:",access_type,"...")
    if "get_insitu_data" in steps_2_run:
        dl=True
    else:
        dl=False
        
    if access_type == 'ARGO_DIRECT':
        df_in_situ_ini,ds_in_situ=get_argo_data_from_direct_access(argo_dir,wmo,workflow_name,dl=dl)
    if access_type == 'ARGO_CERBERE':
        df_in_situ_ini,ds_in_situ=get_argo_data_from_cerbere_access(cerbere_dir,wmo,workflow_name,dl=dl)
    if access_type == 'ARGO_INDEX':
        df_in_situ_ini=get_argo_data_from_index(argo_dir,workflow_name,dl=dl)
        
    if "get_insitu_data" in steps_2_run:
        print("...completed")
    else:
        print("...skipped, data read from local repository")


    # ### III.b - Define needed datasets and variables for Chlorophyll-A
    print("\n#STEP 2: GET WORKFLOW COPERNICUS DATASETS AND VAR for WORFLOW :",workflow_name,"...")
    l_dataset,d_dataset_var=get_workflow_dataset_and_var(workflow_name)
    print("...completed")

    # ### III.c - spatial resolution and boundaries of the copernicus datasets
    print("\n#STEP 3: GET COPERNICUS DATASETS SPATIO-TEMPORAL RESOLUTION ...")
    if "get_copernicus_dataset_resolution" in steps_2_run:
        l_dataset_stf=get_resolution(workflow_name,method=copernicus_method,clear_cache=True,verbose=verbose)
        print("...completed")
    else:
        l_dataset_stf=get_resolution(workflow_name,method=copernicus_method,clear_cache=False,verbose=verbose)
        print("...skipped, resolution read from cache file")


    # Initialise the log file header
    line2write="date;location;group_crit;dataset_id;copernicus_method;record_format;cycle_step;" +\
               "execution_time[s];spatial_extension[square_degrees];temporal_extension[days];cache file size[B]"
    print(line2write)
    file = open(log_file_cop, 'a')
    file.write(line2write + '\n')
    file.close()
    analysis_date=str(np.datetime64('today','D'))

    # If asked, clear and reinitialise the cache file
    if not os.path.exists(cache_dir):
        os.mkdir(cache_dir)

    if (os.path.exists(cache_copernicus_downloaded_data_index)) & (clear_copernicus_downloaded_data_index_cache):
        os.remove(cache_copernicus_downloaded_data_index)
        if verbose:print("the copernicus_downloaded_data_index_cache file was cleared")

    if not os.path.exists(cache_copernicus_downloaded_data_index):
        line2write = "dataset_id;date_min;date_max;lat_min;lat_max;lon_min;lon_max;cross_180;file_name;i_group"
        file = open(cache_copernicus_downloaded_data_index, 'w')
        file.write(line2write + '\n')
        file.close()

    # ### III.d - group extraction by geographical criterion

    #for dataset_id in l_dataset:
    for dataset_id in [l_dataset[0]]:
        
        print("\n#STEP 4: GET REMAINING IN-SITU DATA TO COLOCATE FROM CACHE INDEX OF ALREADY LOCALLY DOWNLOADED COPERNICUS DATA ...")
        # Test the existence of an index-cache file and if it exists, assess the existence of already downloaded data
        if "get_remaining_data_2_colocate_from_cache" in steps_2_run:
            df_to_colocate=get_data_to_colocate(df_in_situ_ini,dataset_id,l_dataset_stf[dataset_id],delta_px,cache_copernicus_downloaded_data_index,verbose=verbose,log4debug=True,log_file_col_1=log_file_col_1,log_file_col_2=log_file_col_2)
            print("...completed")
        else:
            print("...skipped")
        
        # group observation to colocate in spatio-temporal medium cubes
        print("\n#STEP 5: CREATE GROUPS OF IN-SITU OBSERVATIONS USING CLOSE-BY IN SPACE AND TIME CRITERIA ...")
        if "compute_group_of_obs_from_in_situ_data" in steps_2_run:
            stime=time.time()
            #group_of_obs,group_of_obs_too_old,group_of_obs_too_recent=create_obs_groups(gp_crit,dataset_id,l_dataset_stf[dataset_id],df_to_colocate,verbose=verbose)
            group_of_obs,group_of_obs_too_old,group_of_obs_too_recent=create_obs_groups(gp_crit,dataset_id,l_dataset_stf[dataset_id],df_in_situ_ini,verbose=verbose,log4debug=True,log_file_grp=log_file_grp)
            print('Execution time: {0:.1f} s'.format(time.time()-stime))
            # For debug purpose: save variable 
            with open("sav_group_of_obs.pkl", 'wb') as file:
                pickle.dump(group_of_obs, file)
            print("...completed")
        else:
            print("...skipped, reading from saved variable ...")
            with open('sav_group_of_obs.pkl', 'rb') as file:
                group_of_obs = pickle.load(file)
            print("...completed")
        
        
        
        print("\n#STEP 6: DOWNLOAD COPERNICUS DATA USING ",parallelisation, " parallelisation method.")
        
        
        

        print("\n\n Workflow {0:s}; dataset {1:s} ".format(workflow_name,dataset_id))
        print("Variables to extract: ",d_dataset_var[dataset_id])

        # If method is lazy, pre-load index for the dataset
        if copernicus_method=='lazy':
            ds_cop=copernicusmarine.open_dataset(dataset_id=dataset_id)
            lat_cop=ds_cop['latitude']
            lon_cop=ds_cop['longitude']
            dat_cop=ds_cop['time']

        
        if not os.path.exists(outfile_dir):
            os.mkdir(outfile_dir)
        
        start = time.perf_counter()

        
        if parallelisation == 'mpProcess':
            pr={}
        if parallelisation == 'mpAsync':
            pool = Pool()
            res={}
        if parallelisation == 'dask':
            client = LocalCluster().get_client()
            res=[]
        
        n_obs_group=len(group_of_obs)
        
        group_range=range(n_obs_group)
        #group_range=range(i0,i1)
        
        for i_obs_group in group_range:
            
            #print("i_obs_group=",i_obs_group)
            
            if parallelisation == 'no':
                print("\n i_obs_group/n_obs_group = ",i_obs_group+1,"/",n_obs_group, " no parallelisation")
                get_copernicus_data_for_a_group_of_obs(dataset_id,d_dataset_var,i_obs_group,group_of_obs[i_obs_group],df_in_situ_ini,
                                           l_dataset_stf[dataset_id],delta_px,outfile_dir,cache_copernicus_downloaded_data_index,
                                           copernicus_method,indexation_method,record_format,log_file_cop,analysis_date,
                                           location,gp_crit,verbose=verbose)
            if parallelisation == 'mpProcess':
                if (i_obs_group%100 == 0) :print("\n i_obs_group/n_obs_group = ",i_obs_group+1,"/",n_obs_group, " mp parallelisation active - method process")
                pr[i_obs_group]=mp.Process(target=get_copernicus_data_for_a_group_of_obs,
                                       args=(dataset_id,d_dataset_var,i_obs_group,group_of_obs[i_obs_group],df_in_situ_ini,
                                       l_dataset_stf[dataset_id],delta_px,outfile_dir,cache_copernicus_downloaded_data_index,
                                       copernicus_method,indexation_method,record_format,log_file_cop,analysis_date,
                                       location,gp_crit,verbose))

                pr[i_obs_group].start()
                
            if parallelisation == 'mpAsync':
                if (i_obs_group%100 == 0) :print("\n i_obs_group/n_obs_group = ",i_obs_group+1,"/",n_obs_group, " mp parallelisation active - method async")
                res[i_obs_group] = pool.apply_async(get_copernicus_data_for_a_group_of_obs, [dataset_id,d_dataset_var,i_obs_group,group_of_obs[i_obs_group],df_in_situ_ini,
                                   l_dataset_stf[dataset_id],delta_px,outfile_dir,cache_copernicus_downloaded_data_index,
                                   copernicus_method,indexation_method,record_format,log_file_cop,analysis_date,
                                   location,gp_crit,verbose])
                
            if parallelisation == 'dask':
                if (i_obs_group%100 == 0) :print("\n i_obs_group/n_obs_group = ",i_obs_group+1,"/",n_obs_group, " dask parallelisation active")
                res=client.submit(get_copernicus_data_for_a_group_of_obs,dataset_id,d_dataset_var,i_obs_group,group_of_obs[i_obs_group],df_in_situ_ini,
                                   l_dataset_stf[dataset_id],delta_px,outfile_dir,cache_copernicus_downloaded_data_index,
                                   copernicus_method,indexation_method,record_format,log_file_cop,analysis_date,
                                   location,gp_crit,verbose)
        
        for i_obs_group in group_range:
            
            if parallelisation == 'mpProcess':
                pr[i_obs_group].join()

            if parallelisation == 'mpAsync':
                #print(res[i_obs_group])
                ans = res[i_obs_group].get(timeout=600)
        
        if parallelisation == 'dask':
            res = client.gather(res)
            
        finish = time.perf_counter()
        print(f'It took {finish-start:.3f} second(s) to finish')
            
