#!/usr/bin/env python
# coding: utf-8

# COLOCATION TOOL FOR IN-SITU and COPERNICUS GRIDDED PRODUCTS 

# ADAPTED TO ARGO/GLIDER NEEDS FOR CHL WORKFLOW VALIDATION
# Author: D.Dobler (Euro-Argo ERIC)
# Date: 2025-03-10



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
import portalocker
import pickle

# copernicusmarine.login()
# it saved credentials within:
# ~\.copernicusmarine\.copernicusmarine-credentials


def secure_write_log_files_with_parallel_access(log_file_path,line2write):
    """
    Writes a line in an ascii file and handle concurrent access to the file 
    (e.g. in case of parallelization)

    Parameters
    ----------
    log_file_path : str
        Full path name of the ascii file to modify
    line2write : str
        line to append to the file

    """

    
    locked=True
    nb_tries=0
    while (locked == True) and (nb_tries < 100):
        try:
            file = open(log_file_path, 'a')
            portalocker.lock(file, portalocker.LockFlags.EXCLUSIVE)
            file.write(line2write + '\n')
            portalocker.unlock(file)
            file.close()
            locked=False
        except:
            locked=True
        nb_tries=nb_tries+1

def flush_status_in_file(status_file,dictionary,key="",value=""):
    """
    Writes a dictionary content in an ascii file and handle concurrent access to the file 
    (e.g. in case of parallelization)

    Parameters
    ----------
    status_file : str
        Full path name of the ascii file to modify
    dictionary : dictionary
        values must be strings
    key (optional) : string
        key to add or modify to the dictionary
    value (optional): string
        value associated to key
        
    Returns
    -------
    dictionary : dictionary

    """
    #Add the value to key if in arguments
    if key!="":
        dictionary[key]=value
        
    locked=True
    nb_tries=0
    while (locked == True) and (nb_tries < 100):
        try:
            file = open(status_file, 'w')
            portalocker.lock(file, portalocker.LockFlags.EXCLUSIVE)
            for key,value in dictionary.items():
                file.write(f"{key}:{value}\n")
            portalocker.unlock(file)
            file.close()
            locked=False
        except:
            locked=True
        nb_tries=nb_tries+1
        
    return dictionary
    

def get_cms_data(did,var,lonm,lonp,latm,latp,datm,datp,zm,zp,outd,outf):
    """
    Retrieves copernicus data using the subset method. This funcion :
    - forces the download (skip the prompt asking for a user confirmation of the download)
    - overwrites the downloaded file if it already existed locally

    Parameters
    ----------
    did : str
        dataset_id (see copernicus marine service data portfolio)
    var : str array
        variables to extract from the dataset ["var1","var2"]
    lonm : float
        minimum longitude in degrees and in [dataset_lon_min, dataset_lon_max] domain
    lonp : float
        maximum longitude in degrees and in [dataset_lon_min, dataset_lon_max] domain
    latm : float
        minimum latitude in degrees and in [dataset_lat_min, dataset_lat_max] domain
    latp : float
        maximum latitude in degrees and in [dataset_lat_min, dataset_lat_max] domain
    datm : float
        minimum date in "YYYY-MM-DDThh:mm:ss" format and in [dataset_dat_min, dataset_dat_max] domain
    datp : float
        maximum date in "YYYY-MM-DDThh:mm:ss" format and in [dataset_dat_min, dataset_dat_max] domain
    zm : float
        minimum depth in meters and in [dataset_depth_min, dataset_depth_max] domain
    zp : float
        maximum depth in meters and in [dataset_depth_min, dataset_depth_max] domain
    outd : str
        output directory where the downloaded file will be saved
    outf : str
        name of the downloaded file

    """

    if (copernicusmarine.__version__ < "2.0.0"):
        copernicusmarine.subset(
          dataset_id=did, variables=var,
          minimum_longitude=lonm,maximum_longitude=lonp,
          minimum_latitude=latm,maximum_latitude=latp,
          start_datetime=datm,end_datetime=datp,
          minimum_depth=zm,maximum_depth=zp,
          output_filename = outf,output_directory = outd,
          force_download=True, # Important, otherwise a prompt asks for downloading confirmation.
          overwrite_output_data=True # important because if False (default value): when the output file already exists, it adds a (n) at the end. This can prevent from fetching the correct file name
        )
    else:
        copernicusmarine.subset(
          dataset_id=did, variables=var,
          minimum_longitude=lonm,maximum_longitude=lonp,
          minimum_latitude=latm,maximum_latitude=latp,
          start_datetime=datm,end_datetime=datp,
          minimum_depth=zm,maximum_depth=zp,
          output_filename = outf,output_directory = outd,
          overwrite=True, # important because if False (default value): it adds a (n) at the end. This can prevent from fetching the correct file name
          netcdf_compression_level = 5 # set to 0 by default.s
        )
        


def get_resolution(workflow_name,cache_dir,cache_copernicus_resolution_file,clear_cache=False,verbose=False):
    """
    Get the spatio-temporal resolution and boundaries from the copernicus datasets by either:
    - Extracting and locally saving from the copernicus datasets,
    - Reading from previoulsy saved information in the cache file.
    To force new extraction, the previous cache file must be removed, either by deleting it manually
    or by setting the clear_cache parameter to True.

    Parameters
    ----------
    workflow_name : str
        Name of the workflow for which colocation is needed. So far, only "chl" is supported
    cache_dir : str
        as defined in the configuration file: directory where the cache files are saved
    cache_copernicus_resolution_file : str
        as defined in the configuration file: short file name of the cache for resolution and boundaries.
    clear_cache : Boolean
        if clear_cache is set to True: the already existing cache file is removed. 
    verbose : Boolean
        if verbose is set to True: complementary information is printed in the standard output
    
    Returns
    -------
    l_dataset_stf : dictionary
        dictionary with keys : dataset name
                         values : dictionary with keys: spatio-temporal feature (resolution, limits) 
                                            and values: either floats for latitude and longitude or str for dates
        

    """
    # intialisation (list of the datasets spatio-temporal features or stf)
    i_dataset_stf={}
    l_dataset_stf={}

    # workflow datasets and vars:
    l_dataset=cf.l_dataset
    d_dataset_var=cf.d_dataset_var

    # Test is a cache file with value is present
    if not os.path.exists(cache_dir):
        os.mkdir(cache_dir)

    if (os.path.exists(cache_copernicus_resolution_file)) & (clear_cache):
        os.remove(cache_copernicus_resolution_file)
        if verbose:print("the cache file was cleared")
    
    if not os.path.exists(cache_copernicus_resolution_file):

        # Initialise the cache file
        line2write = "dataset_id;reso_lon_deg;reso_lat_deg;spat_lon_min;spat_lon_max;spat_lat_min;spat_lat_max;reso_time_ns;time_min;time_max"
        file = open(cache_copernicus_resolution_file, 'w')
        file.write(line2write + '\n')
        file.close()
        
        for idataset in l_dataset:
            try:
                if verbose:print("\n\n\n ******* Reading spatio-temporal features for " + idataset)
                
                ds=copernicusmarine.open_dataset(dataset_id=idataset)
                
                i_dataset_stf['reso_lon_deg']=ds.attrs['lon_step']
                i_dataset_stf['reso_lat_deg']=ds.attrs['lat_step']
                i_dataset_stf['spat_lon_min']=ds.attrs['geospatial_lon_min']
                i_dataset_stf['spat_lon_max']=ds.attrs['geospatial_lon_max']
                i_dataset_stf['spat_lat_min']=ds.attrs['geospatial_lat_min']
                i_dataset_stf['spat_lat_max']=ds.attrs['geospatial_lat_max']
                # The global attributes for time min and max (time_coverage_start and time_coverage_end) 
                # are wrongly filled: we must compute min and max from the time coordinate.
                # The time_coverage_resolution global attribute is a string, tricky to parse; here again, it is
                # better to recompute from the time coordinate
                i_dataset_stf['time_min']=np.array(ds['time'].min())
                i_dataset_stf['time_max']=np.array(ds['time'].max())
                i_dataset_stf['reso_time_ns']=np.timedelta64(np.mean(np.array(ds['time'][1:])-np.array(ds['time'][:-1])),'ns')
                
                l_dataset_stf[idataset]=i_dataset_stf
                
                ds.close()
                
            except:
                print("ERROR: while downloading data from cmems")
                print(traceback.format_exc())

            try:
                
                reso_time_str=str(np.timedelta64(i_dataset_stf['reso_time_ns'],'ns')/np.timedelta64(1,'ns'))
                line2write = idataset + ";" + str(i_dataset_stf['reso_lon_deg']) + ";" + str(i_dataset_stf['reso_lat_deg']) + ";" + \
                             str(i_dataset_stf['spat_lon_min']) + ";" + str(i_dataset_stf['spat_lon_max']) + ";" + \
                             str(i_dataset_stf['spat_lat_min']) + ";" + str(i_dataset_stf['spat_lat_max']) + ";" + \
                             reso_time_str + ";" + str(i_dataset_stf['time_min']) + ";" + str(i_dataset_stf['time_max'])
                             
                file = open(cache_copernicus_resolution_file, 'a')
                file.write(line2write + '\n')
                file.close()
                
            except:
                print("ERROR: while writing the cache_resolution file: " + cache_copernicus_resolution_file)
                file.close()
                if os.path.exists(cache_copernicus_resolution_file):
                    os.remove(cache_copernicus_resolution_file)
                print(traceback.format_exc())
    else:
        if verbose: print("Reading spatial resolution and boundaries from cache file")
        Reso_index=pd.read_csv(cache_copernicus_resolution_file,sep=";")
        for idataset in l_dataset:
            i_dataset_stf['reso_lon_deg']=Reso_index['reso_lon_deg'][np.where(Reso_index['dataset_id']==idataset)[0][0]]
            i_dataset_stf['reso_lat_deg']=Reso_index['reso_lat_deg'][np.where(Reso_index['dataset_id']==idataset)[0][0]]
            i_dataset_stf['spat_lon_min']=Reso_index['spat_lon_min'][np.where(Reso_index['dataset_id']==idataset)[0][0]]
            i_dataset_stf['spat_lon_max']=Reso_index['spat_lon_max'][np.where(Reso_index['dataset_id']==idataset)[0][0]]
            i_dataset_stf['spat_lat_min']=Reso_index['spat_lat_min'][np.where(Reso_index['dataset_id']==idataset)[0][0]]
            i_dataset_stf['spat_lat_max']=Reso_index['spat_lat_max'][np.where(Reso_index['dataset_id']==idataset)[0][0]]
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


def qc_from_char2int(val):
    """
    This special function was developed specifically to read QC from Argo Parameters as written in Argo NetCDF files.
    It transformed the char array into an integer array

    Parameters
    ----------
    val : char str
        a char string usually containing argo QC values in the range [0 9], it can be used for another purpose, 
        but it must contain characters that can be casted into integers
    
    Returns
    -------
    out_val : int array
        the corresponding int array        

    """
    
    shape_ini=val.shape
    tmp=((np.array(val)).astype('|S1')).tobytes().decode()
    tmp=tmp.replace(' ','6') # beware: this is only for computational issue, QC6 is unused usually
    tmp=list(tmp)
    out_val=np.array(tmp,dtype='int')
    out_val=np.reshape(out_val,shape_ini)
    return out_val


def get_file_from_url(URL,DL_FILE):
    """
    This function get a file accessible from a web query by its URL

    Parameters
    ----------
    URL : string
        the URL of the file to fetch on the web
    DL_FILE : string
        the output file as saved locally, with the complete pathway (path and filename)       

    """
    response = requests.get(URL)
    if response.status_code == 404:
        print("No " + URL + " found in the gdac")
    else:
        open(DL_FILE, "wb").write(response.content)
        print(URL + " found, locally downloaded in " + DL_FILE)


def get_and_read_Argo_meta_index(argo_dir):
    """
    This function locally copies and reads Argo meta index file from the web

    Parameters
    ----------
    argo_dir : string
        repository where the Argo meta file will be saved.
    
    Returns
    -------
    META_index : panda dataframe 
        a panda dataframe with the content of Argo meta index      

    """
    META_index_file=argo_dir + "ar_index_global_meta.csv"
    URL = "https://data-argo.ifremer.fr/ar_index_global_meta.txt"
    get_file_from_url(URL,META_index_file)
    META_index=pd.read_csv(META_index_file,header=9,sep="/",names=['dac','wmo','remaining'],dtype={'dac': 'str', 'wmo': 'str', 'remaining': 'str'})
    return META_index


def get_dac_from_meta_index(argo_dir,wmo):
    """
    This function provides the dac responsible for the decoding of an Argo float by its wmo ID. It uses the argo meta index downloaded from the web.

    Parameters
    ----------
    argo_dir : string
        repository where the Argo meta file will be saved.
    wmo : string
        the wmo id of the float
    
    Returns
    -------
    dac : string
        the corresponding data assembly center (dac)      
        
    """
    META_index=get_and_read_Argo_meta_index(argo_dir)
    dac=META_index['dac'][np.where(META_index['wmo']==wmo)[0][0]]
    return dac


def get_argo_data_from_direct_access(argo_dir,wmo,workflow_name,dl=True):
    """
    This function optionally downloads Argo data multi-profile file from the web and/or
    returns a dataframe and a dataset with all the necessary information regarding
    colocation

    Parameters
    ----------
    argo_dir : string
        the local repository where argo data files are stored    
    wmo : string
        the wmo id of the float
    workflow_name : string
        the name of the workflow (only 'chl' is supported so far)
    dl (optional): boolean (set to True by default)
        if set to True, Argo data are collected on the GDAC at https://data-argo.ifremer.fr/
        if set to False, a local copy should already exists.
    
    Returns
    -------
    df : panda dataframe
        this dataframe only contains the necessary information for colocation: cycle, date, lat, lon, date_qc and position_qc
    ds : xarray dataset
        this dataset contains df information and pressure and chlorophyll-a arrays.
        
    """
    SPROF_FILE = argo_dir + "/" + wmo + "_Sprof.nc"
    
    if dl:
        
        dac = get_dac_from_meta_index(argo_dir,wmo)
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
    DIRECTION=ds.variables['DIRECTION'][:]
    if workflow_name == 'chl':
        CHLA=ds.variables['CHLA'][:]
        CHLA_QC=ds.variables['CHLA_QC'][:]
        
    #print(ds)
    #ds.close()
    
    ref_date=np.datetime64("1950-01-01T00:00:00")
    dates=JULD*86400*np.timedelta64(1, 's')+ref_date
    
    # output for colocation computation (1-D)
    df=pd.DataFrame({'CYCLE':cycles,'DATE': dates, 'LAT': latitudes, 'LON': longitudes, 'DATE_QC': dates_qc, 'POSITION_QC': position_qc})

    # output for colocation display
    n_prof,n_levels=PRES.shape
    prof=np.arange(n_prof)
    levels=np.arange(n_levels)
    
    
    if workflow_name == 'chl':
        ds = xr.Dataset(
            data_vars=dict(
                DATE=(["prof"], dates),
                LAT=(["prof"], latitudes,{'units':ds.variables['LATITUDE'].units}),
                LON=(["prof"], longitudes,{'units':ds.variables['LONGITUDE'].units}),
                CYCLE=(["prof"], cycles),
                DIRECTION=(["prof"], DIRECTION),
                DATE_QC=(["prof"], dates_qc),
                POSITION_QC=(["prof"], position_qc),
                PRES=(["prof", "levels"], PRES),
                CHLA=(["prof", "levels"], CHLA,{'units':ds.variables['CHLA'].units}),
                CHLA_QC=(["prof", "levels"], CHLA_QC),
                
            ),
            coords=dict(
                prof=prof,
                n_levels=levels,
            ),
            attrs=dict(description="Observation related data"),
        )
    else:
        ds = xr.Dataset(
            data_vars=dict(
                DATE=(["prof"], dates),
                LAT=(["prof"], latitudes,{'units':ds.variables['LATITUDE'].units}),
                LON=(["prof"], longitudes,{'units':ds.variables['LONGITUDE'].units}),
                CYCLE=(["prof"], cycles),
                DATE_QC=(["prof"], dates_qc),
                POSITION_QC=(["prof"], position_qc),
                PRES=(["prof", "levels"], PRES),
            ),
            coords=dict(
                prof=prof,
                n_levels=levels,
            ),
            attrs=dict(description="Observation related data"),
        )

    return df,ds


def get_argo_data_from_cerbere_access(cerbere_dir,wmo,workflow_name,dl=False):
    """
    This function optionally copies cerbere data files from xxx (not yet plugged) and/or
    returns a dataframe and a dataset with all the necessary information regarding
    colocation

    Parameters
    ----------
    cerbere_dir : string
        the local repository where cerbere data files are stored    
    wmo : string
        the wmo id of the float
    workflow_name : string
        the name of the workflow (only 'chl' is supported so far)
    dl (optional) : boolean (set to False by default)
        if set to True, cerbere data are collected on xxx (not yet plugged)
        if set to False, a local copy should already exists.
    
    Returns
    -------
    df : panda dataframe
        this dataframe only contains the necessary information for colocation: cycle, date, lat, lon, date_qc and position_qc
    ds : xarray dataset
        this dataset contains df information and pressure and chlorophyll-a arrays.
        
    """
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
    
    if workflow_name == 'chl':
        ds = xr.Dataset(
            data_vars=dict(
                DATE=(["prof"], dates),
                LAT=(["prof"], latitudes,{'units':ds.variables['LATITUDE'].units}),
                LON=(["prof"], longitudes,{'units':ds.variables['LONGITUDE'].units}),
                CYCLE=(["prof"], cycles),
                DATE_QC=(["prof"], dates_qc),
                POSITION_QC=(["prof"], position_qc),
                PRES=(["prof", "levels"], PRES),
                CHLA=(["prof", "levels"], CHLA,{'units':ds.variables['CHLA'].units}),
                CHLA_QC=(["prof", "levels"], CHLA_QC),
                
            ),
            coords=dict(
                prof=prof,
                n_levels=levels,
            ),
            attrs=dict(description="Observation related data"),
        )
    else:
        ds = xr.Dataset(
            data_vars=dict(
                DATE=(["prof"], dates),
                LAT=(["prof"], latitudes,{'units':ds.variables['LATITUDE'].units}),
                LON=(["prof"], longitudes,{'units':ds.variables['LONGITUDE'].units}),
                CYCLE=(["prof"], cycles),
                DATE_QC=(["prof"], dates_qc),
                POSITION_QC=(["prof"], position_qc),
                PRES=(["prof", "levels"], PRES),
            ),
            coords=dict(
                prof=prof,
                n_levels=levels,
            ),
            attrs=dict(description="Observation related data"),
        )

    return df,ds


# ## II.d - get all observations for one workflow


def get_argo_data_from_index(argo_dir,workflow_name,dl=True,verbose=False):
    """
    This function optionally downloads Argo data bio-profile index file from the web and/or
    returns a dataframe with all the necessary information regarding
    colocation

    Parameters
    ----------
    argo_dir : string
        the local repository where argo data files are stored
    workflow_name : string
        the name of the workflow (only 'chl' is supported so far)
    dl (optional) : boolean (set to True by default)
        if set to True, Argo bio-profile index file is collected the GDAC at https://data-argo.ifremer.fr/
        if set to False, a local copy should already exists.
    
    Returns
    -------
    df : panda dataframe
        this dataframe only contains the necessary information for colocation: cycle, date, lat, lon, date_qc and position_qc
        p.s. There is no xarray dataset outputs as the pressure and chlorophyll-A information are not within the index file.
        
    """
    BIO_Index_file=argo_dir + "/argo_bio-profile_index.csv"
    
    if dl:
        URL = "https://data-argo.ifremer.fr/argo_bio-profile_index.txt"
        get_file_from_url(URL,BIO_Index_file)
        print("Downloading in-situ data from ", URL)
        
    
    if verbose: print("Reading in-situ data from ", BIO_Index_file)
    BIO_Index=pd.read_csv(BIO_Index_file,header=8,sep=",")

    # Removing lines with incomplete coordinates:
    #print("intial size : ",BIO_Index.shape)
    BIO_Index.drop(BIO_Index.index[BIO_Index['date'].isnull()],inplace=True)
    #print("after removing null dates : ",BIO_Index.shape)
    BIO_Index.drop(BIO_Index.index[BIO_Index['latitude'].isnull()],inplace=True)
    #print("after removing null latitudes : ",BIO_Index.shape)
    BIO_Index.drop(BIO_Index.index[BIO_Index['longitude'].isnull()],inplace=True)
    #print("after removing null longitudes : ",BIO_Index.shape)

    if workflow_name == 'chl':
        # Keeping lines including the worklow parameter:
        BIO_Index=BIO_Index[BIO_Index['parameters'].str.contains("CHLA")]
        #print("after selecting chla lines : ",BIO_Index.shape)
  
   
    latitudes=np.array(BIO_Index['latitude'])
    longitudes=np.array(BIO_Index['longitude'])
    dates=BIO_Index['date'].astype(str)
    dates_str_iso_8601=dates.str[:4]+"-"+dates.str[4:6]+"-"+dates.str[6:8]+"T"+dates.str[8:10]+":"+dates.str[10:12]+":"+dates.str[12:14]
    dates_dt64=np.array(dates_str_iso_8601,dtype='datetime64')
    #print(dates_dt64)
    
    # as the index file already filters on good qc, values are set to 1:
    dates_qc=np.ones(dates_dt64.shape)
    position_qc=np.ones(dates_dt64.shape)
    
    df=pd.DataFrame({'DATE': dates_dt64, 'LAT': latitudes, 'LON': longitudes, 'DATE_QC': dates_qc, 'POSITION_QC': position_qc})

    return df


# ### II.e - Distance computation function

def compute_earth_radius_elliptical(lat_deg):
    """
    This function returns the earth radius at a given latitude, assuming an
    elliptical earth model.

    Parameters
    ----------
    lat_deg : float
        the latitude in degrees
    
    Returns
    -------
    earth_radius_m : float
        the corresponding earth radius in meters
        
    """
    # 
    
    if type(lat_deg)==np.ndarray:
        lat_deg=lat_deg.astype('float64')
        
    a=6378137 # equatorial radius
    b=6356750 # polar radius
    e=np.sqrt(1-(b/a)**2)
    lat_rad=lat_deg*np.pi/180
    earth_radius_m=a*np.sqrt(1-e**2*(np.sin(lat_rad))**2)
    
    return earth_radius_m
    
def compute_distance(lonA_deg=0,latA_deg=0,lonB_deg=1,latB_deg=0,verbose=False):
    """
    This function computes the distance between two Earth positions or set of positions A and B given known geographical coordinates.
    if both A and B are sets of positions (arrays), then A and B should have the same dimension and the distance is point by point distance. 
    if A (resp. B) is an array and B (resp. A) is a float, then the distance is computed between all A posisions (resp. B) and B points
    if A and B are floats: the distance between the 2 positions is computed.

    Parameters
    ----------
    lonA_deg : float or float array
        the longitude of point A or points A array in degrees
    latA_deg : float or float array
        the latitude of point A or points A array in degrees
    lonB_deg : float or float array
        the longitude of point B or points B array in degrees
        
    latB_deg : float or float array
        the latitude of point B or points B array in degrees
    verbose (optional): boolean (set to False by default)
        if set to True: additionnal printing are made in the standard output
        else: no print
    
    Returns
    -------
    distance_m : float of float array
        the corresponding distance in meters
        
    """
    
    # force float64 for input data to deal with default
    # ndarray dtype which is float32
    # and in this case, the computation is done in float32 
    # which can lead to up to 8% relative
    # error a distance of 1/12 deg (8/9 km).

    lonA_deg=np.array(lonA_deg).astype('float64')
    latA_deg=np.array(latA_deg).astype('float64')
    lonB_deg=np.array(lonB_deg).astype('float64')    
    latB_deg=np.array(latB_deg).astype('float64')
    
    
    #then compute earth_radius_m median
    #Earth_radius_m=6376*10**3 # in [m]
    lat_med_deg=0.5*(latA_deg+latB_deg)
    Earth_radius_m=compute_earth_radius_elliptical(lat_med_deg)
    #print(lat_med_deg,Earth_radius_m)
    
    # first, put them in radians
    lonA_rad=lonA_deg*np.pi/180
    latA_rad=latA_deg*np.pi/180
    lonB_rad=lonB_deg*np.pi/180
    latB_rad=latB_deg*np.pi/180
    #print(lonA_rad,latA_rad,lonB_rad,latB_rad)
    
    if ((len(lonA_deg.shape)!=0) & (len(lonB_deg.shape) !=0)):
        if (len(lonA_deg) > len(lonB_deg)): distance_m=np.zeros(lonA_deg.shape)
        else: distance_m=np.zeros(lonB_deg.shape)
    if ((len(lonA_deg.shape)==0) & (len(lonB_deg.shape) !=0)):
        distance_m=np.zeros(lonB_deg.shape)
    if ((len(lonA_deg.shape)!=0) & (len(lonB_deg.shape) ==0)):
        distance_m=np.zeros(lonA_deg.shape)
    if ((len(lonA_deg.shape)==0) & (len(lonB_deg.shape) ==0)): 
        distance_m=0.0
    
    eps=1e-13
    is_A_an_array=False
    is_B_an_array=False
    try: 
        if np.size(lonA_deg) > 1:is_A_an_array=True
    except: 
        eps=eps
    try: 
        if np.size(lonB_deg) > 1:is_B_an_array=True
    except: 
        eps=eps
    
    if verbose:
        print("is_A_an_array,is_B_an_array:")
        print(is_A_an_array,is_B_an_array)
    
    
    if ((is_A_an_array==True) & (is_B_an_array==True)):
        #check where there is equality:
        i_neq=np.where((abs(lonA_deg-lonB_deg)>eps) | (abs(latA_deg-latB_deg) > eps))

        # then compute distance_m in [m]
        distance_m[i_neq]=Earth_radius_m[i_neq]*np.arccos(np.sin(latA_rad[i_neq])*np.sin(latB_rad[i_neq]) + \
                                 np.cos(latA_rad[i_neq])*np.cos(latB_rad[i_neq])*np.cos(lonB_rad[i_neq]-lonA_rad[i_neq]))
    
        
    if ( (is_A_an_array==False) & (is_B_an_array==True)):
        #check where there is equality:
        i_neq=np.where((abs(lonA_deg-lonB_deg)>eps) | (abs(latA_deg-latB_deg) > eps))
        # then compute distance_m in [m]
        AA=np.sin(latA_rad)*np.sin(latB_rad[i_neq])
        BB=np.cos(latA_rad)*np.cos(latB_rad[i_neq])*np.cos(lonB_rad[i_neq]-lonA_rad)
        cos_val=AA+BB
        distance_m[i_neq]=Earth_radius_m[i_neq]*np.arccos(cos_val)
    
    
    if ((is_A_an_array==True) & (is_B_an_array==False)):
        #check where there is equality:
        i_neq=np.where((abs(lonA_deg-lonB_deg)>eps) | (abs(latA_deg-latB_deg) > eps))

        # then compute distance_m in [m]
        distance_m[i_neq]=Earth_radius_m[i_neq]*np.arccos(np.sin(latA_rad[i_neq])*np.sin(latB_rad) + \
                                 np.cos(latA_rad[i_neq])*np.cos(latB_rad)*np.cos(lonB_rad-lonA_rad[i_neq]))
    
    
    if ((is_A_an_array==False) & (is_B_an_array==False)):
        if (abs(lonA_deg-lonB_deg)>eps) | (abs(latA_deg-latB_deg) > eps):
            distance_m=Earth_radius_m*np.arccos(np.sin(latA_rad)*np.sin(latB_rad)+ 
                                 np.cos(latA_rad)*np.cos(latB_rad)*np.cos(lonB_rad-lonA_rad))
    
            
    
    
    return distance_m


def get_bbox_from_df(df,i_dataset_stf,delta_px,verbose=False):
    
    """
    For a given copernicus dataset, this function computes the spatio-temporal bounding box around df coordinates, accounting for 
    dataset resolution (i_dataset_stf) and additionnal number of pixels (on the east/west/north/east/before/after: delta_px)
    It handles when df only contains 1 point.
    It handles 180° crossing and limit cases such as when it is crossed adding delta_px but not sufficiently to reach the first pixel on the other side.

    Parameters
    ----------
    df : panda dataframe
        output of either get_argo_data_from_index, get_argo_data_from_cerbere_access or get_argo_data_from_direct_access functions
        it shall contain the necessary information for colocation: cycle, date, lat, lon, date_qc and position_qc
    i_dataset_stf : dictionary 
        spatio-temporal resolution and limits associated to the dataset. It means the function is called for a given copernicus dataset
        with keys : spatio-temporal feature (resolution, limits) 
        and values : either floats for latitude and longitude or str for dates
        keys and values should be as created by the get_resolution function.
    delta_px : dictionary
        as defined in the configuration file. keys are 'x', 'y' and 't' strings and values are integers
    verbose (optional) : boolean (set to False by default)
        if set to True: additionnal printing are made in the standard output
        else: no print
    
    Returns
    -------
    bbox : dictionary
        with keys: spatio-temporal feature (bbox_dates_min, bbox_dates_max, bbox_lat_min, bbox_lat_max, bbox_lon_west, bbox_lon_east, 
                                            cross_180, spatial_extension_square_deg, temporal_extension_days) 
        and values: - floats for latitudes, longitudes and spatio-temporal extensions 
                    - string for dates (ISO 8601) 
                    - int for cross_180
        cross_180 is set to 1 if the bounding box is across the 180° longitude line, with bbox_lon_west positive and bbox_lon_east negative. It is set to 0 otherwise.
        spatial_extension_square_deg and temporal_extension_days are not really used afterwards, here for performance assessment and cross-checking.
        
    """
    
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
    # Assumption: the medium cube will never reach 180° large
    # To prevent headaches with min and max, the bbox boundaries for longitude are expressed
    # as western boundary (bbox_lon_west) and eastern boundary (bbox_lon_east), making it clear
    # of what limits we are talking about. 
    if lonp-lonm < 180:
        
        cross_180=0
        bbox_lon_west=lonm-delta_lon
        bbox_lon_east=lonp+delta_lon
        
        #print("1",bbox_lon_west,bbox_lon_east)
        
        # Here we handle boundaries conditions if copernicus longitude limits are 
        # crossed while applying delta_lon
        if (bbox_lon_west < -180 - i_dataset_stf['reso_lon_deg']):
            cross_180 = 1
            bbox_lon_west = bbox_lon_west + 360
            
        #print("2",bbox_lon_west,bbox_lon_east)

        if (bbox_lon_east > 180 + i_dataset_stf['reso_lon_deg']):
            cross_180 = 1
            bbox_lon_east = bbox_lon_east - 360
        #print("3",bbox_lon_west,bbox_lon_east)
        
        # Here is handled when limits are overcrossed but not enough to reach the first pixel on the other side.
        if ((bbox_lon_west >= -180 - i_dataset_stf['reso_lon_deg']) & (bbox_lon_west < i_dataset_stf['spat_lon_min'])):
            bbox_lon_west = i_dataset_stf['spat_lon_min']
        if ((bbox_lon_east <= 180  + i_dataset_stf['reso_lon_deg']) & (bbox_lon_east > i_dataset_stf['spat_lon_max'])):
            bbox_lon_east = i_dataset_stf['spat_lon_max']
        
    else:
        cross_180=1
        lonp=np.min(lon_ok[lon_ok>0])
        lonm=np.max(lon_ok[lon_ok<0])
        
        # When 180 is crossed and longitude domain is ]-180,180], 
        # the bbox_lon_west is positive and bbox_lon_east is negative
        bbox_lon_west=lonp-delta_lon
        bbox_lon_east=lonm+delta_lon
        
        # Here are handled boundary conditions if the 180 crossing is too small
        # with respect to copernicus longitude limits (i.e. the small part in the 
        # other side of 180° line is smaller than the resolution of Copernicus dataset).
        if bbox_lon_west > i_dataset_stf['spat_lon_max']:
            cross_180=0
            bbox_lon_west=i_dataset_stf['spat_lon_min']
            
        if bbox_lon_east < i_dataset_stf['spat_lon_min']:
            cross_180=0
            bbox_lon_east=i_dataset_stf['spat_lon_max']
        
        #print("4",bbox_lon_west,bbox_lon_east)
    
    
    
    dlon=bbox_lon_east-bbox_lon_west
        
    # compute time boundaries
    bbox_dates_min=str(np.datetime_as_string(datm-np.timedelta64(delta_px['t'],'D'),'s'))
    bbox_dates_max=str(np.datetime_as_string(datp+np.timedelta64(delta_px['t'],'D'),'s'))
    # account for copernicus dataset time boundaries
    bbox_dates_min=np.datetime_as_string(max(np.datetime64(bbox_dates_min),np.datetime64(i_dataset_stf['time_min'])),unit='s')
    bbox_dates_max=np.datetime_as_string(min(np.datetime64(bbox_dates_max),np.datetime64(i_dataset_stf['time_max'])),unit='s')

    # compute spatio-temporal extensions:
    spatial_extension_square_deg=dlon*dlat
    temporal_extension_days=np.timedelta64(datp-datm,'D') / np.timedelta64(1, 'D')
    spatial_extension_square_deg_max=max(spatial_extension_square_deg,spatial_extension_square_deg_max)
    temporal_extension_days_max=max(temporal_extension_days_max,temporal_extension_days)

    # some debug printing:
    if verbose:
        print("obs_lat_min     = {0:.2f}\t\t, obs_lat_max     = {1:.2f}".format(latm,latp))
        print("bbox_lat_min    = {0:.2f}\t\t, bbox_lat_max    = {1:.2f}\t dlat = {2:.2f}".format(bbox_lat_min,bbox_lat_max,dlat))
        print("obs_lon_east     = {0:.2f}\t\t, obs_lon_west     = {1:.2f}".format(lonm,lonp))
        print("bbox_lon_west    = {0:.2f}\t\t, bbox_lon_east    = {1:.2f}\t dlon = {2:.2f}".format(bbox_lon_west,bbox_lon_east,dlon))
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
    bbox['bbox_lon_west']=bbox_lon_west
    bbox['bbox_lon_east']=bbox_lon_east
    bbox['cross_180']=cross_180
    bbox['spatial_extension_square_deg']=spatial_extension_square_deg
    bbox['temporal_extension_days']=temporal_extension_days

    return bbox

def get_data_to_colocate(df,dataset_id,i_dataset_stf,delta_px,cache_copernicus_downloaded_data_index,
                         verbose=False,log4debug=False,log_file_col_1="",log_file_col_2=""):
    
    """
    This function computes the remaining data to colocate given a copernicus dataset, data coordinates, a cache for already downloaded colocated data

    Parameters
    ----------
    df : panda dataframe
        output of either get_argo_data_from_index, get_argo_data_from_cerbere_access or get_argo_data_from_direct_access functions
        it shall contain the necessary information for colocation: cycle, date, lat, lon, date_qc and position_qc
    dataset_id : string
        as defined in the l_dataset output of get_workflow_dataset_and_var function
    i_dataset_stf : dictionary 
        spatio-temporal resolution and limits associated to the dataset. It means the function is called for a given copernicus dataset
        with keys: spatio-temporal feature (resolution, limits) 
        and values: either floats for latitude and longitude or str for dates
        keys and values should be as created by the get_resolution function.
    delta_px : dictionary
        as defined in the configuration file. keys are 'x', 'y' and 't' strings and values are integers
    cache_copernicus_downloaded_data_index : string
        complete filename with path where the cache file is stored. The file should include the following fields:
        'dataset_id': 'str', 'date_min': 'str', 'date_max': 'str', 'lat_min' : 'float','lat_max':'float',
        'lon_west' : 'float','lon_east':'float', 'cross_180' : 'int','file_name':'str','i_group':'int'
    verbose (optional) : boolean (set to False by default)
        if set to True: additionnal printing are made in the standard output
        else: no print
    log4debug (optional) : boolean (set to False by default)
        if set to True, create additionnal log files, used to test the robustness of the function.
    log_file_col_1 (optional) : string (set to "" by default)
        complete filename with path where the input df data will be stored if log4debug is set to True
    log_file_col_2 (optional) : string (set to "" by default)
        complete filename with path where the results of search in cache will be stored if log4debug is set to True
    
    Returns
    -------
    df_to_colocate : panda dataframe 
        subset of df, containing subset of df data that are not in cache 
        it contains the necessary information for colocation: cycle, date, lat, lon, date_qc and position_qc
    colocated_files : dictionary
        it contains, for each observation with colocated data in local cache repository, the full name of the associated colocated file.
    """
    
    #initialise dictionary matching colocated file:
    colocated_files={}
    
    # Test whether the information are not yet downloaded:
    cache_index=pd.read_csv(cache_copernicus_downloaded_data_index,sep=";",dtype={'dataset_id': 'str', 'date_min': 'str', 'date_max': 'str',
                                                                'lat_min' : 'float','lat_max':'float','lon_west' : 'float','lon_east':'float',
                                                                'cross_180' : 'int','file_name':'str','i_group' : 'int'})

    
    # extract the observations from df:
    dates_qc=np.array(df['DATE_QC'])
    position_qc=np.array(df['POSITION_QC'])
    latitudes=np.array(df['LAT'])
    longitudes=np.array(df['LON'])
    dates=np.array(df['DATE'],dtype='datetime64')
    
    n_obs=np.size(dates)
    obs_ids=np.arange(n_obs)
    
    # if log4debug == True:
        # file=open(log_file_col_1,'w')
        # file.write("i_obs;date;lat;lon;date_qc;pos_qc\n")
        # file.close()

        # for i_obs in obs_ids:
            # df_iobs=df.iloc[i_obs]
            # line2write_fmt="{0:d};{1:s};{2:.6f};{3:.6f};{4:.0f};{5:.0f}"
            # line2write=line2write_fmt.format(i_obs,str(df_iobs['DATE']),df_iobs['LAT'],df_iobs['LON'],df_iobs['DATE_QC'],df_iobs['POSITION_QC'])
            # if (i_obs%100 == 0) :print(line2write)
            # file=open(log_file_col_1,'a')
            # file.write(line2write + '\n')
            # file.close()
    
    # keep only correctly located observations, and within the dataset time domain.
    i_good_obs=np.where(   ((position_qc==1) | (position_qc==2) | (position_qc==5) | (position_qc==8)) & \
                           ((dates_qc==1) | (dates_qc==2) | (dates_qc==5) | (dates_qc==8))  )[0]
    
    df_ok=df.iloc[i_good_obs]
    obs_ids_ok=obs_ids[i_good_obs]
    
    
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
        for i_obs in obs_ids_ok:
            
            df_iobs=df_ok.iloc[i_obs]
            
            bbox=get_bbox_from_df(df_iobs,i_dataset_stf,delta_px)
            
            
            if (df_iobs['DATE'] <= np.datetime64(i_dataset_stf['time_min'])) | (df_iobs['DATE'] >= np.datetime64(i_dataset_stf['time_max'])):
                in_cache="no"
                cache_file=""
                bbox['bbox_dates_max']=df_iobs['DATE']
            else:
     
                # dataset_id;date_min;date_max;lat_min;lat_max;lon_west;lon_east;cross_180;file_name
                
                # test_presence needs to deal with cyclic longitude as well, which can be tedious and errors-prone
                # There are four cases to handle:
                # case 1: 180 is not crossed, both obs bbox x-limits are inside [lon_west (min), lon_east (max)] in the cache index.
                # case 2: 180 is crossed, both obs bbox x-limits are inside [-180,lon_east (max)] in the cache index.
                # case 3: 180 is crossed, both obs bbox x-limits are inside [lon_west (min), 180] in the cache index.
                # case 4: 180 is crossed, obs bbox lon_west are inside [lon_west (min), 180] and obs bbox lon_east is inside [-180,lon_east (max)]
                # N.B. case 1 and 4 could be merged in fact, but I let them like this for undersantding purpose.
                
                cache_west_lim=cache_index['lon_west']-eps
                cache_east_lim=cache_index['lon_east']+eps
                
                test_presence=np.where((cache_index['dataset_id'] == dataset_id) &
                                   (cache_index['date_min']   <= bbox['bbox_dates_min']) &
                                   (cache_index['date_max']   >= bbox['bbox_dates_max']) &
                                   (cache_index['lat_min']-eps    <= bbox['bbox_lat_min']) &
                                   (cache_index['lat_max']+eps    >= bbox['bbox_lat_max']) &
                                   ( ((cache_west_lim    <= bbox['bbox_lon_west']) & # case 1
                                      (cache_east_lim    >= bbox['bbox_lon_west']) & # case 1
                                      (cache_west_lim    <= bbox['bbox_lon_east']) & # case 1
                                      (cache_east_lim    >= bbox['bbox_lon_east']) & # case 1
                                      (cache_index['cross_180']  == 0)) |            # case 1
                                     ((cache_east_lim    >= bbox['bbox_lon_west']) & # case 2
                                      (cache_east_lim    >= bbox['bbox_lon_east']) & # case 2
                                      (cache_index['cross_180']  == 1)) |            # case 2
                                     ((cache_west_lim    <= bbox['bbox_lon_west']) & # case 3
                                      (cache_west_lim    <= bbox['bbox_lon_east']) & # case 3
                                      (cache_index['cross_180']  == 1)) |            # case 3
                                     ((cache_west_lim    <= bbox['bbox_lon_west']) & # case 4
                                      (cache_east_lim    >= bbox['bbox_lon_east']) & # case 4
                                      (cache_index['cross_180']  == 1)) )            # case 4
                                   )[0]

                
                #print(str(bbox['bbox_dates_min']))
                if test_presence.size == 0:
                    in_cache="no"
                    cache_file=""
                    if len(i_obs_to_colocate)==0:
                        i_obs_to_colocate=[i_obs]
                    else:
                        i_obs_to_colocate.append(i_obs)
                else:
                    in_cache=str(cache_index['i_group'][test_presence[0]])
                    cache_file=cache_index['file_name'][test_presence[0]]
                    colocated_files[i_obs]=cache_file

            line2write_fmt="{0:d};{1:s};{2:.6f};{3:.6f};{4:s};{5:s};{6:.6f};{7:.6f};{8:.6f};{9:.6f};{10:d};{11:s};{12:s}\n"
            line2write=line2write_fmt.format(i_obs,str(df_iobs['DATE']),df_iobs['LAT'],df_iobs['LON'],str(bbox['bbox_dates_min']),
                                             str(bbox['bbox_dates_max']),bbox['bbox_lat_min'],bbox['bbox_lat_max'],bbox['bbox_lon_west'],
                                             bbox['bbox_lon_east'],bbox['cross_180'],in_cache,cache_file)
            
            if log4debug:
                file=open(log_file_col_2,'a')
                file.write(line2write)
                file.close()
                if (i_obs%100 == 0): print(line2write)
                
            
            if verbose & (not log4debug):print(line2write)
                                   
        
        print("There are {0:d}/{1:d} points already colocated.".format(len(colocated_files.keys()),n_obs))
        print("There remains {0:d}/{1:d} points to colocate.".format(len(i_obs_to_colocate),n_obs))
        print("The remaining points, if existant, are too recent to colocate with copernicus dataset")
        df_to_colocate=df_ok.iloc[i_obs_to_colocate]
        
    else:
        print("There is nothing in the " + cache_copernicus_downloaded_data_index + "file, all data need to be colocated")
        df_to_colocate=df_ok
        

    return df_to_colocate,colocated_files



def create_obs_groups(df_in_situ_ini,gp_crit,i_dataset_stf,verbose=False,log4debug=False,log_file_grp=""):
    """
    This function creates groups of close-by observations

    Parameters
    ----------
    df_in_situ_ini : panda dataframe
        output of either get_argo_data_from_index, get_argo_data_from_cerbere_access or get_argo_data_from_direct_access functions
        it shall contain the necessary information for colocation: cycle, date, lat, lon, date_qc and position_qc
    gp_crit : dictionary
        grouping criteria in number of close-by pixels, as defined in the configuration file. 
        keys are 'gp_max_x_n', 'gp_max_y_n' and 'gp_max_t_n' strings and values are integers.
    i_dataset_stf : dictionary 
        spatio-temporal resolution and limits associated to the dataset. It means the function is called for a given copernicus dataset
        with keys: spatio-temporal feature (resolution, limits) 
        and values: either floats for latitude and longitude or str for dates
        keys and values should be as created by the get_resolution function.
    verbose (optional) : boolean (set to False by default)
        if set to True: additionnal printing are made in the standard output
        else: no print
    log4debug (optional) : boolean (set to False by default)
        if set to True, create additionnal log files, used to test the robustness of the function.
    log_file_grp (optional) : string (set to "" by default)
        complete filename with path where the results of grouping will be stored if log4debug is set to True
    
    Returns
    -------
    group_of_obs : dictionary
        keys are integer, corresponding to the group id.
        values are integer arrays containing the list of df_in_situ_ini indexes (and not keys) corresponding to the group id
    group_of_obs_too_old : integer array
        df_in_situ_ini indexes for which the date if before the first date of the copernicus dataset
    group_of_obs_too_recent : integer array
        df_in_situ_ini indexes for which the date if after the last date of the copernicus dataset
    
    """
    # outputs initialisation
    # create a dictionary with the indexes of the various groups
    group_of_obs={}
    group_of_obs_too_old={}
    group_of_obs_too_recent={}
    

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
    n_elt_ini=len(list_obs_id)
    #list_obs_id = df_in_situ_ini.keys()
    print("Initial number of observation:", n_elt_ini)
    
    #initialise log file for investigations
    if log4debug:
        file=open(log_file_grp,'w')
        file.write("i_obs;date;lat;lon;date_qc;pos_qc;id_group;comment;nb_elt_group\n")
        file.close()
    
    i_group=0
    


    # First discard observations that are too old for the copernicus dataset:
    print("Discarding observations that are too old")
    i_too_old=np.where( (df_in_situ_ini['DATE'] < i_dataset_stf['time_min']) )[0]
    i_obs_group_n=list_obs_id[i_too_old]
    group_of_obs_too_old=i_obs_group_n
    list_obs_id=np.delete(list_obs_id,i_too_old)
    print("Number of too old obs = {0:d}; Remaining obs to group = {1:d}".format(len(i_obs_group_n),len(list_obs_id)))
    if log4debug:
        for i_obs_too_old in i_obs_group_n:
            line2write_fmt="{0:d};{1:s};{2:.4f};{3:.4f};{4:.0f};{5:.0f};no_group;too_old\n"
            line2write=line2write_fmt.format(i_obs_too_old,str(df_in_situ_ini['DATE'][i_obs_too_old]),df_in_situ_ini['LAT'][i_obs_too_old],
                                             df_in_situ_ini['LON'][i_obs_too_old],df_in_situ_ini['DATE_QC'][i_obs_too_old],
                                             df_in_situ_ini['POSITION_QC'][i_obs_too_old])
            file=open(log_file_grp,'a')
            file.write(line2write)
            file.close()
        

    # Second discard observations that are too recent for the copernicus dataset:
    print("Discarding observations that are too recent")
    i_too_recent=np.where( (df_in_situ_ini['DATE'] > i_dataset_stf['time_max']) )[0]
    i_obs_group_n=list_obs_id[i_too_recent]
    group_of_obs_too_recent=i_obs_group_n
    list_obs_id=np.delete(list_obs_id,i_too_recent)
    print("Number of too recent obs = {0:d}; Remaining obs to group = {1:d}".format(len(i_obs_group_n),len(list_obs_id)))
    if log4debug:
        for i_obs_too_recent in i_obs_group_n:
            line2write_fmt="{0:d};{1:s};{2:.4f};{3:.4f};{4:.0f};{5:.0f};no_group;too_recent\n"
            line2write=line2write_fmt.format(i_obs_too_recent,str(df_in_situ_ini['DATE'][i_obs_too_recent]),df_in_situ_ini['LAT'][i_obs_too_recent],
                                             df_in_situ_ini['LON'][i_obs_too_recent],df_in_situ_ini['DATE_QC'][i_obs_too_recent],
                                             df_in_situ_ini['POSITION_QC'][i_obs_too_recent])
            file=open(log_file_grp,'a')
            file.write(line2write)
            file.close()

    earth_radius_eq=compute_earth_radius_elliptical(0)
    
    while (len(list_obs_id) > 0) & (i_group<=df_in_situ_ini.shape[0]) :
    #while (len(list_obs_id) > 0) & (i_group<=600) :

        df_0=df_in_situ_ini.iloc[list_obs_id[0]]
        lon = df_0['LON']
        lat = df_0['LAT']
        dat = df_0['DATE']
        # lon = df_in_situ_ini['LON'][list_obs_id[0]]
        # lat = df_in_situ_ini['LAT'][list_obs_id[0]]
        # dat = df_in_situ_ini['DATE'][list_obs_id[0]]
    
        df_remain=df_in_situ_ini.iloc[list_obs_id]
        lon_obs_left=df_remain['LON']
        lat_obs_left=df_remain['LAT']
        dat_obs_left=df_remain['DATE']
        # lon_obs_left=df_in_situ_ini['LON'][list_obs_id]
        # lat_obs_left=df_in_situ_ini['LAT'][list_obs_id]
        # dat_obs_left=df_in_situ_ini['DATE'][list_obs_id]
    
        #compute equatorial equivalent distances
        dist_m_2_deg_at_equat=(180/(np.pi*earth_radius_eq))
        dist_lon_deg=compute_distance(lon,0,lon_obs_left,np.zeros(lon_obs_left.shape)) * dist_m_2_deg_at_equat
        dist_lat_deg=compute_distance(0,lat,np.zeros(lat_obs_left.shape),lat_obs_left) * dist_m_2_deg_at_equat
        dist_time_dt64=abs(dat-dat_obs_left)
    
        if verbose:
            print("First observation position: {0:}, {1:.3f}°N {2:.3f}°E".format(dat,lat,lon))
            print(dist_lon_deg[:5])
            print((lon-df_remain['LON'])[:5])
            print(dist_lat_deg[:5])
            print((lat-df_remain['LAT'])[:5])
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
                df_i_obs=df_in_situ_ini.iloc[i_obs]
                line2write_fmt="{0:d};{1:s};{2:.4f};{3:.4f};{4:.0f};{5:.0f};{6:d};;{7:d}\n"
                line2write=line2write_fmt.format(i_obs,str(df_i_obs['DATE']),df_i_obs['LAT'], 
                                                 df_i_obs['LON'],df_i_obs['DATE_QC'],
                                                 df_i_obs['POSITION_QC'],i_group,len(i_close_by))
                file=open(log_file_grp,'a')
                file.write(line2write)
                file.close()
                
        completion_rate=100*(1-len(list_obs_id)/n_elt_ini)
        if(i_group%100==0):print("i_group={0:d};nb_elt_group={1:d};n_elt_left_to_group={2:d};completion={3:.1f}%".format(i_group,len(i_close_by),len(list_obs_id),completion_rate))
            
        i_group = i_group + 1  
    if verbose: print(group_of_obs)
        
    return group_of_obs,group_of_obs_too_old,group_of_obs_too_recent



def get_copernicus_data_for_a_group_of_obs(dataset_id,d_dataset_var,i_obs_group,group_of_obs,df_in_situ_ini,i_dataset_stf,delta_px,outfile_dir,
                                           cache_copernicus_downloaded_data_index,copernicus_method,indexation_method,record_format,log_file_cop,
                                           analysis_date,location,gp_crit,spa_lim,tim_lim,verbose=True):
    """
    This function uses copernicusmarine library to download colocated copernicus data for a given group of observations. To download the data, several options 
    are possible: "subset" method or "lazy" loading. When lazy loading is selected, the index selection method can be tuned ("sel", "isel", or "index") 
    and the output format as well ("values", "NetCDF" or a "computation"). N.B.: the lazy loading method was implemented and tested for performance assessment. However, it did not 
    prove to be more performant, thus the default choice should be the subset method.

    Parameters
    ----------
    dataset_id : string
        as defined in the l_dataset output of get_workflow_dataset_and_var function
    d_dataset_var : 
        dictionary with keys: string from l_dataset
                         values: string array with the variables to extract from the dataset.
        as defined in the d_dataset_var output of get_workflow_dataset_and_var function
    i_obs_group : integer
        specific key of group_of_obs dictionary to handle
    group_of_obs : dictionary
        keys are integer, corresponding to the group id.
        values are integer arrays containing the list of df_in_situ_ini indexes corresponding to the group id (please note these
        are indexes of the df_in_situ_ini dataframe, not keys).
    df_in_situ_ini : panda dataframe
        output of either get_argo_data_from_index, get_argo_data_from_cerbere_access or get_argo_data_from_direct_access functions
        it shall contain the necessary information for colocation: cycle, date, lat, lon, date_qc and position_qc
    i_dataset_stf : dictionary 
        spatio-temporal resolution and limits associated to the dataset. It means the function is called for a given copernicus dataset
        with keys: spatio-temporal feature (resolution, limits) 
        and values: either floats for latitude and longitude or str for dates
        keys and values should be as created by the get_resolution function.
    delta_px : dictionary
        as defined in the configuration file. keys are 'x', 'y' and 't' strings and values are integers
    outfile_dir : string
        as defined in the configuration file. Repository where copernicus colocated files will be stored
    cache_copernicus_downloaded_data_index : string
        as defined in the configuration file. Complete filename. The cache file indexes copernicus files.  
    copernicus_method : string
        as defined in the configuration file: either 'lazy' or 'subset' (n.b. performance tests proved that lazy loading was not more performant)
    indexation_method : string
        as defined in the configuration file: either 'sel' or 'isel' or 'index' (only used in case of copernicus_method = 'lazy')
    record_format : string
        as defined in the configuration file: either 'values' or 'NetCDF' or 'computation' (only used in case of copernicus_method = 'lazy')
    log_file_cop : string
        complete filename with path where performance information will be logged (n.b. in case of parallel processing, this has to be corrected).
    analysis_date : string
        date of analysis for log purpose
    location : string
        location or condition when the test was run, for log purpose
    gp_crit : dictionary
        grouping criteria in number of close-by pixels, as defined in the configuration file. 
        keys are 'gp_max_x_n', 'gp_max_y_n' and 'gp_max_t_n' strings and values are integers.
        only use in this function for log purpose
    spa_lim : float
        limit in square degrees for the spatial area covered by the group of observation above which the download is not performed (secure if issue in groups)
    tim_lim : float
        limit in days for the time span covered by the group of observation above which the download is not performed (secure if issue in groups)
    verbose (optional) : boolean (set to False by default)
        if set to True: additionnal printing are made in the standard output
        else: no print
    
    Returns
    -------
    none
    
    """
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
    bbox_lon_west=bbox['bbox_lon_west']
    bbox_lon_east=bbox['bbox_lon_east']
    cross_180=bbox['cross_180']
    spatial_extension_square_deg=bbox['spatial_extension_square_deg']
    temporal_extension_days=bbox['temporal_extension_days']
    
    if (spatial_extension_square_deg > spa_lim) or (temporal_extension_days > tim_lim):
        print("WARNING: the size of box containing the group of observations is above the configured limit; \n " + \
              "rerun Step 5: create groups of in-situ data to make sure your groups are clean")
        line2write_fmt="{0:s};{1:s};{2:d}_{3:d}_{4:d};{5:s};{6:s};{7:s};{8:.0f};{9:.5f};{10:.2f};{11:.2f};{12:s}"
        line2write=line2write_fmt.format(analysis_date,location,gp_crit['gp_max_x_n'],gp_crit['gp_max_y_n'],gp_crit['gp_max_t_n'],dataset_id,"",
                                         record_format,i_obs_group,time.time()-stime,spatial_extension_square_deg,temporal_extension_days,
                                         "WARNING: the size of box containing the group of observations is above the configured limit; copernicus data not downloaded.")
        print(line2write)
        secure_write_log_files_with_parallel_access(log_file_cop,line2write)
        return

    if copernicus_method == 'lazy':
        print("Subsetting data with the 'lazy' load")
        outfile_name="{0:s}_{1:s}_{2:s}_{3:.1f}_{4:.1f}_{5:.1f}_{6:.1f}.nc".format(dataset_id,bbox_dates_min[:10],bbox_dates_max[:10],bbox_lat_min,
                                                                        bbox_lat_max,bbox_lon_west,bbox_lon_east)
        print("outfile_name=",outfile_name)
        
        
        ii_dat=np.where((dat_cop > np.datetime64(bbox_dates_min)) & (dat_cop < np.datetime64(bbox_dates_max)))
        ii_lat=np.where((lat_cop > bbox_lat_min) & (lat_cop < bbox_lat_max))
        if cross_180 == 0:
            ii_lon=np.where((lon_cop > bbox_lon_west) & (lon_cop < bbox_lon_east))
        if cross_180 == 1:
            ii_lon=np.where((lon_cop < bbox_lon_east) | (lon_cop > bbox_lon_west))

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
        
        outfile_name="{0:s}_{1:s}_{2:s}_{3:.1f}_{4:.1f}_{5:.1f}_{6:.1f}.nc".format(dataset_id,tmin[:10],tmax[:10],ymin,ymax,bbox_lon_west,bbox_lon_east)
        print("outfile_name=",outfile_name)
        if cross_180 == 1 :
            # There is a need to split the request in two:
            
            outfile_name_1=outfile_name + "_1.nc"
            xmin=i_dataset_stf['spat_lon_min']
            xmax=min(bbox_lon_east,i_dataset_stf['spat_lon_max'])
            #print("cross_180=1, part 1",bbox_lon_west,bbox_lon_east,i_dataset_stf['spat_lon_min'],i_dataset_stf['spat_lon_max'],xmin,xmax)
            if verbose:
                print("180 was crossed, splitting the subset request")
                print("First request arguments:")
                print(dataset_id,d_dataset_var[dataset_id],
                     xmin,xmax,ymin,ymax,tmin,tmax,0,0,
                     outfile_dir,outfile_name_1)
            get_cms_data(dataset_id,d_dataset_var[dataset_id],xmin,xmax,ymin,ymax,tmin,tmax,0,0,outfile_dir,outfile_name_1)
            
            outfile_name_2=outfile_name + "_2.nc"
            xmin=max(bbox_lon_west,i_dataset_stf['spat_lon_min'])
            xmax=i_dataset_stf['spat_lon_max']
            #print("cross_180=1, part 2",xmin,xmax)
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
            
            # delete temporary files
            os.remove(outfile_dir + "/" + outfile_name_1)
            os.remove(outfile_dir + "/" + outfile_name_2)
        else:
            xmin=max(bbox_lon_west,i_dataset_stf['spat_lon_min'])
            xmax=min(bbox_lon_east,i_dataset_stf['spat_lon_max'])
            #print("cross_180=0",xmin,xmax)
            get_cms_data(dataset_id,d_dataset_var[dataset_id],xmin,xmax,ymin,ymax,tmin,tmax,0,0,outfile_dir,outfile_name)

    
    # Saving in the cache file 
    #"dataset_id;date_min;date_max;lat_min;lat_max;lon_west;lon_east;cross_180;file_name"
    line2write_fmt="{0:s};{1:s};{2:s};{3:.6f};{4:.6f};{5:.6f};{6:.6f};{7:d};{8:s};{9:d};{10:s}"
    line2write=line2write_fmt.format(dataset_id,tmin,tmax,ymin,ymax,bbox_lon_west,bbox_lon_east,cross_180,outfile_name,i_obs_group,analysis_date)
    secure_write_log_files_with_parallel_access(cache_copernicus_downloaded_data_index,line2write)

    
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
    secure_write_log_files_with_parallel_access(log_file_cop,line2write)


def get_copernicus_mini_cubes(workflow_name,df_obs,ds_obs,dataset_id,i_dataset_stf,delta_px,outfile_dir,cache_copernicus_downloaded_data_index,d_dataset_var,verbose=False):
    """
    This function compute mini-cubes around observations using medium cubes already locally downloaded in cache.

    Parameters
    ----------
    workflow_name : string
        the name of the workflow (only 'chl' is supported so far)
    df_obs : panda dataframe
        output of either get_argo_data_from_index, get_argo_data_from_cerbere_access or get_argo_data_from_direct_access functions
        it shall contain the necessary information for colocation: cycle, date, lat, lon, date_qc and position_qc
    ds_obs : xarray dataset
        this dataset contains df information and pressure and chlorophyll-a arrays.
    dataset_id : string
        as defined in the l_dataset output of get_workflow_dataset_and_var function
    i_dataset_stf : dictionary 
        spatio-temporal resolution and limits associated to the dataset. It means the function is called for a given copernicus dataset
        with keys: spatio-temporal feature (resolution, limits) 
        and values: either floats for latitude and longitude or str for dates
        keys and values should be as created by the get_resolution function.
    delta_px : dictionary
        as defined in the configuration file. keys are 'x', 'y' and 't' strings and values are integers
    outfile_dir : string
        as defined in the configuration file. Repository where copernicus colocated files are stored
    cache_copernicus_downloaded_data_index : string
        as defined in the configuration file. Complete filename. The cache file indexes copernicus files. 
    d_dataset_var : 
        dictionary with keys: string from l_dataset
                         values: string array with the variables to extract from the dataset.
        as defined in the d_dataset_var output of get_workflow_dataset_and_var function
    verbose (optional) : boolean (set to False by default)
        if set to True: additionnal printing are made in the standard output
        else: no print
    
    Returns
    -------
    colocated_data : dictionary of dictionaries of xarray datasets
        keys are 'copernicus' and 'insitu'
        values are dictionnaries with 
            keys are observation ids (prof index in ds_obs)
            values are xarray datasets containing colocated mini-cube from copernicus dataset or in-situ observation dataset (coordinates and workflow parameter value)
        
    
    """
    colocated_data={}
    colocated_copernicus_data={}
    colocated_insitu_data={}
    
    # For each observation in df_obs, find the corresponding line in the cache index and the corresponding copernicus NetCDF local copy.
    df_to_colocate,colocated_files = get_data_to_colocate(df_obs,dataset_id,i_dataset_stf,delta_px,cache_copernicus_downloaded_data_index)
    
    delta_lon=i_dataset_stf['reso_lon_deg']*delta_px['x']
    delta_lat=i_dataset_stf['reso_lon_deg']*delta_px['y']
    
    #print(colocated_files)
    
    for i_obs in colocated_files.keys():
        
        if (i_obs%100 == 0):print("i_obs=",i_obs)
        
        #read copernicus medium cube
        ds_cop=xr.open_dataset(outfile_dir+"/"+colocated_files[i_obs])
        # Followinglines discarded, were needed when NetCDF4 library was used
        # ds_cop.set_auto_mask(False)
        # cop_dat=np.array(nc.num2date(ds_cop.variables['time'][:], ds_cop.variables['time'].units, only_use_cftime_datetimes=False)).astype('datetime64[ms]') 
        cop_lat=ds_cop.variables['latitude'][:]
        cop_lon=ds_cop.variables['longitude'][:]
        cop_dat=ds_cop.variables['time'][:]

        ds_iobs=ds_obs.isel(prof=i_obs)
        #read obs coordinates:
        obs_dat_qc=np.array(ds_iobs['DATE_QC'])
        obs_pos_qc=np.array(ds_iobs['POSITION_QC'])
        obs_lat=np.array(ds_iobs['LAT'])
        obs_lon=np.array(ds_iobs['LON'])
        obs_dat=np.array(ds_iobs['DATE'],dtype='datetime64')
        
        if workflow_name == 'chl':
            obs_chl=np.array(ds_iobs['CHLA'])
            obs_chl_qc=np.array(ds_iobs['CHLA_QC'])

        #print(obs_lon,cop_lon)
        dist_lon=compute_distance(obs_lon,0,cop_lon,np.zeros(cop_lon.shape))
        dist_lat=compute_distance(0,obs_lat,np.zeros(cop_lat.shape),cop_lat)
        dist_dat=abs(cop_dat-obs_dat)

        dist_lon_thresh=compute_distance(delta_lon,0,0,0)
        dist_lat_thresh=compute_distance(0,delta_lat,0,0)
        dist_dat_thresh=i_dataset_stf['reso_time_ns']*delta_px['t']

        i_dat_close_by=np.where((dist_dat<dist_dat_thresh))[0]
        i_lat_close_by=np.where((dist_lat<dist_lat_thresh))[0]
        i_lon_close_by=np.where((dist_lon<dist_lon_thresh))[0]
        
        # print("\n ds_cop")
        # print(ds_cop)
        ds_cop_colocated=ds_cop.isel(time=i_dat_close_by,latitude=i_lat_close_by,longitude=i_lon_close_by)
        # print("\n ds_cop_colocated")
        # print(ds_cop_colocated)
        
        # print("\n",i_obs,obs_dat,obs_lat,obs_lon)
        # print(colocated_files[i_obs])
        # print(i_dat_close_by)
        # print(cop_dat[i_dat_close_by])
        # print(i_lat_close_by)
        # print(cop_lat[i_lat_close_by])
        # print(i_lon_close_by)
        # print(cop_lon[i_lon_close_by])
        
        colocated_copernicus_data[i_obs]=ds_cop_colocated
        colocated_insitu_data[i_obs]=ds_iobs
    
    
    colocated_data['copernicus']=colocated_copernicus_data
    colocated_data['insitu']=colocated_insitu_data
        
    # print("\n HEHE : On y est enfin, après un peu de sueur ...")
    # print("colocated_data")
    # print(colocated_data)

    return colocated_data

if __name__ == '__main__':
    
    """
    This library extract data from copernicus marine service, colocated with in-situ observations.
    There are 9 steps:
    Step 1: retrieve parameters from the configuration file Colocation_cfg.py
    Step 2: get in-situ data
    Step 3: get copernicus dataset spatio-temporal resolution and limits
    Step 4: get remaining data to colocate
    Step 5: create groups of in-situ data (medium-cubes definition)
    Step 6: download copernicus data
    Step 7: extract mini-cube from downloaded data
    Step 8: a few plots - TO BE CODED
    Steps can be run at demand
    
    There are 3 temporary arguments for debug purposes, which will be discarded in the end: 
     - the parallelisation method either 'no' or 'mpProcess' or 'mpAsync' or 'dask'
     - the first group to treat
     - the last group to treat
     
    Parameters
    ----------
    none
    
    Returns
    -------
    none
    
    """
    
    

    print("\n#STEP 1: READING CONFIGURATION") 
    import Colocation_cfg as cf
    
    # input data selection
    access_type=cf.access_type
    cerbere_dir=cf.cerbere_dir
    argo_dir=cf.argo_dir
    wmo=cf.wmo
    
    # colocation parameterization
    workflow_name=cf.workflow_name
    l_dataset=cf.l_dataset
    d_dataset_var=cf.d_dataset_var
    gp_crit=cf.gp_crit
    delta_px=cf.delta_px
    spa_lim=cf.spatial_extension_square_deg_limit
    tim_lim=cf.temporal_extension_days_limit
    copernicus_method=cf.copernicus_method
    indexation_method=cf.indexation_method
    record_format=cf.record_format
    outdir_cop=cf.outdir_cop
    outdir_col_plots=cf.outdir_col_plots
    outfig_dir=outdir_col_plots + cf.wmo + "/"
    if not os.path.exists(outdir_cop):os.mkdir(outdir_cop)
    if not os.path.exists(outdir_col_plots):os.mkdir(outdir_col_plots)
    if not os.path.exists(outfig_dir):os.mkdir(outfig_dir)
    
    # applying constraints
    if copernicus_method == 'subset':
        record_format="NetCDF"
        indexation_method=""
        
    # paralellisation option
    parallelisation=cf.parallelisation
    i0=cf.grp_deb
    grp_end=cf.grp_end
    igrp_2_colocate=cf.igrp_2_colocate
    # debug parameterization:
    #i0,i1=int(sys.argv[2]),int(sys.argv[3])
    #parallelisation = sys.argv[1]
    
    # Steps_to_run
    steps_2_run=cf.steps_2_run
    
    # cache definition
    cache_dir=cf.cache_dir
    cache_copernicus_resolution_file=cf.cache_copernicus_resolution_file
    cache_copernicus_downloaded_data_index=cf.cache_copernicus_downloaded_data_index
    cache_group_of_obs_prefix=cf.cache_group_of_obs_prefix
    clear_cache_copernicus_resolution=cf.clear_cache_copernicus_resolution
    clear_cache_copernicus_downloaded_data_index=cf.clear_cache_copernicus_downloaded_data_index
    clear_cache_group_of_obs=cf.clear_cache_group_of_obs
    
    # additionnal logs 
    log4debug=cf.log4debug
    log_file_cop=cf.log_file_cop
    log_file_col_1_prefix=cf.log_file_col_1_prefix
    log_file_col_2_prefix=cf.log_file_col_2_prefix
    log_file_grp_prefix=cf.log_file_grp_prefix
    location=cf.location
    
    # standard output log
    verbose=cf.verbose
    
    
    # Initialise the status codes
    status_file=cf.status_file
    if os.path.exists(status_file):
        os.remove(status_file)
    ST_notstarted='0'
    ST_started='1'
    ST_completed='2' 
    ST_skipped='6'
    ST_error='404'
    
    status={}
    status['stat_step_1']=ST_completed
    status['stat_step_2']=ST_notstarted
    status['stat_step_3']=ST_notstarted
    status['stat_step_4']=ST_notstarted
    status['stat_step_5']=ST_notstarted
    status['stat_step_6']=ST_notstarted
    status['stat_step_7']=ST_notstarted
    status['stat_step_8']=ST_notstarted
    status=flush_status_in_file(status_file,status)


    if verbose:
        print("Estimate of the number of copernicus points to fetch: {:d}".format(gp_crit['gp_max_x_n']*gp_crit['gp_max_y_n']*gp_crit['gp_max_t_n']))

    

    # ### III.b - IN-SITU data selection
    print("\n#STEP 2: GET IN SITU DATA FROM:",access_type,"...")
    if 2 in steps_2_run:
        dl=True
        status['stat_step_2']=ST_started
    else:
        dl=False
        status['stat_step_2']=ST_skipped
    flush_status_in_file(status_file,status)
        
    if access_type == 'ARGO_DIRECT':
        df_in_situ_ini,ds_in_situ=get_argo_data_from_direct_access(argo_dir,wmo,workflow_name,dl=dl)
    if access_type == 'ARGO_CERBERE':
        df_in_situ_ini,ds_in_situ=get_argo_data_from_cerbere_access(cerbere_dir,wmo,workflow_name,dl=dl)
    if access_type == 'ARGO_INDEX':
        # When access_type is ARGO_INDEX, no ds is output as it is not intended to create and draw for the whole index.
        # Additionnaly, only the index is used, which means no data is read. If reading was required, this would mean opening
        # a fair number of files.
        df_in_situ_ini=get_argo_data_from_index(argo_dir,workflow_name,dl=dl)
        
    if 2 in steps_2_run:
        print("...completed")
        status['stat_step_2']=ST_completed
    else:
        print("...skipped, data read from local repository")
    flush_status_in_file(status_file,status)


    # ### III.c - spatial resolution and boundaries of the copernicus datasets
    print("\n#STEP 3: GET COPERNICUS DATASETS SPATIO-TEMPORAL RESOLUTION ...")
    status['stat_step_3']=ST_started
    flush_status_in_file(status_file,status)
    l_dataset_stf=get_resolution(workflow_name,cache_dir,cache_copernicus_resolution_file,clear_cache=clear_cache_copernicus_resolution,verbose=verbose)
    if not clear_cache_copernicus_resolution:
        print("...completed, resolution downloaded from copernicus")
    else:
        print("...completed, resolution read from cache file")
    status['stat_step_3']=ST_completed
    flush_status_in_file(status_file,status)

    
    # Initialise the performance log file header
    line2write="date;location;group_crit;dataset_id;copernicus_method;record_format;cycle_step;" +\
               "execution_time[s];spatial_extension[square_degrees];temporal_extension[days];cache file size[B]"
    file = open(log_file_cop, 'a')
    file.write(line2write + '\n')
    file.close()
    analysis_date=str(np.datetime64('today','D'))

    # If asked, clear and reinitialise the cache file
    if not os.path.exists(cache_dir):
        os.mkdir(cache_dir)

    if (os.path.exists(cache_copernicus_downloaded_data_index)) & (clear_cache_copernicus_downloaded_data_index):
        os.remove(cache_copernicus_downloaded_data_index)
        if verbose:print("the copernicus_downloaded_data_index_cache file was cleared")

    if not os.path.exists(cache_copernicus_downloaded_data_index):
        line2write = "dataset_id;date_min;date_max;lat_min;lat_max;lon_west;lon_east;cross_180;file_name;i_group"
        file = open(cache_copernicus_downloaded_data_index, 'w')
        file.write(line2write + '\n')
        file.close()

    # ### III.d - group extraction by geographical criterion
    status_downloading_colocated_data_percent={}
    for dataset_id in l_dataset:
    #for dataset_id in [l_dataset[0]]:
        i_dataset_stf=l_dataset_stf[dataset_id]
        
        outfile_dir=outdir_cop + workflow_name + "_" + cf.l_dataset_short_name[dataset_id] + "/"
        if not os.path.exists(outfile_dir):os.mkdir(outfile_dir)
        
        print("\n#STEP 4: GET REMAINING IN-SITU DATA TO COLOCATE FROM CACHE INDEX OF ALREADY LOCALLY DOWNLOADED COPERNICUS DATA ...")
        # Test the existence of an index-cache file and if it exists, assess the existence of already downloaded data
        if 4 in steps_2_run:
            status=flush_status_in_file(status_file,status,'stat_step_4',ST_started)
            status=flush_status_in_file(status_file,status,'stat_step_4_'+ dataset_id,ST_started)
            df_to_colocate,colocated_files =get_data_to_colocate(df_in_situ_ini,dataset_id,i_dataset_stf,
                                                delta_px,cache_copernicus_downloaded_data_index,verbose=verbose,log4debug=log4debug,
                                                log_file_col_1=log_file_col_1_prefix + dataset_id + ".csv",log_file_col_2=log_file_col_2_prefix + dataset_id + ".csv")
            status=flush_status_in_file(status_file,status,'stat_step_4_'+ dataset_id,ST_completed)
            print("...completed")
        else:
            
            df_to_colocate=df_in_situ_ini
            status=flush_status_in_file(status_file,status,'stat_step_4_'+ dataset_id,ST_skipped)
            print("...skipped")
        
        # group observation to colocate in spatio-temporal medium cubes
        print("\n#STEP 5: CREATE GROUPS OF IN-SITU OBSERVATIONS USING CLOSE-BY IN SPACE AND TIME CRITERIA ...")
        if (2 in steps_2_run) or (4 in steps_2_run) or (5 in steps_2_run) or (not os.path.exists(cache_group_of_obs_prefix + dataset_id + ".pkl")):
            stime=time.time()
            status=flush_status_in_file(status_file,status,'stat_step_5',ST_started)
            status=flush_status_in_file(status_file,status,'stat_step_5_'+ dataset_id,ST_started)
            
            group_of_obs,group_of_obs_too_old,group_of_obs_too_recent=create_obs_groups(df_to_colocate,gp_crit,i_dataset_stf,
                                                                                        verbose=verbose,log4debug=log4debug,
                                                                                        log_file_grp=log_file_grp_prefix + dataset_id + ".csv")
            print('Execution time: {0:.1f} s'.format(time.time()-stime))
            # For debug purpose: save variable 
            with open(cache_group_of_obs_prefix + dataset_id + ".pkl", 'wb') as file:
                pickle.dump(group_of_obs, file)
            
            status=flush_status_in_file(status_file,status,'stat_step_5_'+ dataset_id,ST_completed)
            print("...completed")
        else:
            
            with open(cache_group_of_obs_prefix  + dataset_id + ".pkl", 'rb') as file:
                group_of_obs = pickle.load(file)
            
            status=flush_status_in_file(status_file,status,'stat_step_5_'+ dataset_id,ST_skipped)
            print("...skipped, read from saved variable ...")
            
        
        
        print("\n#STEP 6: DOWNLOAD COPERNICUS DATA USING ", parallelisation, " parallelisation method.")
        status=flush_status_in_file(status_file,status,'stat_step_6',ST_started)
        status=flush_status_in_file(status_file,status,'stat_step_6_'+ dataset_id,ST_started)
        
        print("\n Workflow {0:s}; dataset {1:s} ".format(workflow_name,dataset_id))
        print("Variables to extract: ",d_dataset_var[dataset_id])

        # If method is lazy, pre-load index for the dataset
        if copernicus_method=='lazy':
            ds_cop=copernicusmarine.open_dataset(dataset_id=dataset_id)
            lat_cop=ds_cop['latitude']
            lon_cop=ds_cop['longitude']
            dat_cop=ds_cop['time']

        
        start = time.perf_counter()

        
        if parallelisation == 'mpProcess':
            pr={}
        if parallelisation == 'mpAsync':
            # define the number of core to use.
            # Here, we keep 2 cores for other works
            nb_processes=max(1,mp.cpu_count()-2)
            pool = Pool(processes=nb_processes)
            res={}
        # if parallelisation == 'dask':
            # client = LocalCluster().get_client()
            # res=[]
        
        n_obs_group=len(group_of_obs)
        
        if np.size(igrp_2_colocate)>0:
            group_range=igrp_2_colocate
        else:
            if grp_end==-1: 
                i1=n_obs_group
            else:
                i1=min(grp_end,n_obs_group)
                #print(i1,grp_end,n_obs_group)
            group_range=range(i0,i1)
        
        #print("group_range=",group_range)
        n_group_range=len(group_range)
        
        i=-1
        for i_obs_group in group_range:
            i=i+1
            
            if parallelisation == 'no':
                print("\n i_obs_group/n_obs_group = ",i_obs_group+1,"/",n_obs_group, " no parallelisation")
                get_copernicus_data_for_a_group_of_obs(dataset_id,d_dataset_var,i_obs_group,group_of_obs[i_obs_group],df_to_colocate,
                                           i_dataset_stf,delta_px,outfile_dir,cache_copernicus_downloaded_data_index,
                                           copernicus_method,indexation_method,record_format,log_file_cop,analysis_date,
                                           location,gp_crit,spa_lim,tim_lim,verbose=verbose)
            if parallelisation == 'mpProcess':
                if (i%100 == 0) :print("\n i_obs_group/n_obs_group = ",i_obs_group+1,"/",n_obs_group, " mp parallelisation active - method process")
                pr[i_obs_group]=mp.Process(target=get_copernicus_data_for_a_group_of_obs,
                                       args=(dataset_id,d_dataset_var,i_obs_group,group_of_obs[i_obs_group],df_to_colocate,
                                       i_dataset_stf,delta_px,outfile_dir,cache_copernicus_downloaded_data_index,
                                       copernicus_method,indexation_method,record_format,log_file_cop,analysis_date,
                                       location,gp_crit,spa_lim,tim_lim,verbose))

                pr[i_obs_group].start()
                
            if parallelisation == 'mpAsync':
                if (i%100 == 0) :print("\n i_obs_group/n_obs_group = ",i_obs_group+1,"/",n_obs_group, " mp parallelisation active - method async")
                res[i_obs_group] = pool.apply_async(get_copernicus_data_for_a_group_of_obs, [dataset_id,d_dataset_var,i_obs_group,group_of_obs[i_obs_group],
                                                                                             df_to_colocate,i_dataset_stf,delta_px,outfile_dir,
                                                                                             cache_copernicus_downloaded_data_index,copernicus_method,
                                                                                             indexation_method,record_format,log_file_cop,analysis_date,
                                                                                             location,gp_crit,spa_lim,tim_lim,verbose])
        
        i=0
        for i_obs_group in group_range:
            
            if parallelisation == 'mpProcess':
                pr[i_obs_group].join()

            if parallelisation == 'mpAsync':
                #print(res[i_obs_group])
                ans = res[i_obs_group].get(timeout=600)

            if (i%10==0):
                status_downloading_colocated_data_percent[dataset_id]=100*i/n_group_range
                
                completion_rate="{:.1f} %".format(status_downloading_colocated_data_percent[dataset_id])
                status=flush_status_in_file(status_file,status,'stat_step_6_' + dataset_id + "_percent",completion_rate)
            i=i+1
        
        # if parallelisation == 'dask':
            # res = client.gather(res)
        
        status=flush_status_in_file(status_file,status,'stat_step_6_' + dataset_id + "_percent","100 %")
        status=flush_status_in_file(status_file,status,'stat_step_6_' + dataset_id,ST_completed)
        finish = time.perf_counter()
        print(f'It took {finish-start:.3f} second(s) to finish')

    if 4 in steps_2_run:
        status=flush_status_in_file(status_file,status,'stat_step_4',ST_completed)
    status=flush_status_in_file(status_file,status,'stat_step_5',ST_completed)
    status=flush_status_in_file(status_file,status,'stat_step_6',ST_completed)

    
    for dataset_id in l_dataset:
        print("\n#STEP 7: EXTRACT MINI-CUBES around observations...")
        print(dataset_id)
        
        outfile_dir=outdir_cop + workflow_name + "_" + cf.l_dataset_short_name[dataset_id] + "/"
        if not os.path.exists(outfile_dir):os.mkdir(outfile_dir)
        
        if (7 in steps_2_run) and (access_type != 'ARGO_INDEX'):
            status=flush_status_in_file(status_file,status,'stat_step_7',ST_started)
            status=flush_status_in_file(status_file,status,'stat_step_7_' + dataset_id,ST_started)
            colocated_data=get_copernicus_mini_cubes(workflow_name,df_in_situ_ini,ds_in_situ,dataset_id,i_dataset_stf,delta_px,outfile_dir,cache_copernicus_downloaded_data_index,d_dataset_var,verbose=verbose)
            status=flush_status_in_file(status_file,status,'stat_step_7_' + dataset_id,ST_completed)

        else:
            if (7 in steps_2_run) and (access_type == 'ARGO_INDEX'):
                print("ERROR: DO NOT EXTRACT MINI-CUBES FOR THE WHOLE INDEX, please")
                status=flush_status_in_file(status_file,status,'stat_step_7_' + dataset_id,ST_error)
                status=flush_status_in_file(status_file,status,'stat_step_7',ST_error)
            else:
                print(" ... skipped")
                status=flush_status_in_file(status_file,status,'stat_step_7_' + dataset_id,ST_skipped)
                
                
        print("\n#STEP 8: display observations...")
        
        if (8 in steps_2_run) & (7 in steps_2_run) & (access_type != 'ARGO_INDEX'):
            status=flush_status_in_file(status_file,status,'stat_step_8',ST_started)
            status=flush_status_in_file(status_file,status,'stat_step_8_' + dataset_id,ST_started)
            
            n_obs=len(colocated_data['copernicus'].keys())
            
            
            # value to extract:
            cop_var_name=d_dataset_var[dataset_id][0]
            cop_var_MIN,cop_var_MAX=99999,-99999
            for i_obs in range(n_obs):
                var_val=(colocated_data['copernicus'][i_obs])[cop_var_name].values
                cop_var_MIN=np.nanmin([np.nanmin(var_val),cop_var_MIN])
                cop_var_MAX=np.nanmax([np.nanmax(var_val),cop_var_MAX])
            cop_var_all_values=(colocated_data['copernicus'][0])[cop_var_name].values
            cop_var_UNIT=(colocated_data['copernicus'][0])[cop_var_name].units
            
            print((colocated_data['insitu'][0])["CHLA"])
            obs_var_UNIT=(colocated_data['insitu'][0])["CHLA"].units
            
            for i_obs in range(n_obs):
                
                #obs point:
                #ds_iobs=ds_in_situ.isel(prof=i_obs)
                ds_iobs=colocated_data['insitu'][i_obs]
                #read obs coordinates:
                obs_dati=np.array(ds_iobs['DATE'],dtype='datetime64')
                obs_cyci=np.array(ds_iobs['CYCLE'])
                obs_diri=np.array(ds_iobs['DIRECTION'],dtype='|S1').tobytes().decode()
                obs_lati=np.array(ds_iobs['LAT'])
                obs_loni=np.array(ds_iobs['LON'])
                obs_chli=np.array(ds_iobs['CHLA'])
                obs_prei=np.array(ds_iobs['PRES'])
                obs_dat=[obs_dati, obs_dati]
                obs_lat=[obs_lati, obs_lati]
                obs_lon=[obs_loni, obs_loni]
                
                i_notnan=np.where(obs_chli < 99999)
                
                print("len(i_notnan[0])=",len(i_notnan[0]))
                if len(i_notnan[0])==0:
                    obs_prei_notnan=[np.nan]
                    obs_chli_notnan=[np.nan]
                else:
                    obs_prei_notnan=obs_prei[i_notnan]
                    obs_chli_notnan=obs_chli[i_notnan]
                
                obs_chl=[obs_chli_notnan[0], obs_chli_notnan[0]]

                #copernicus points:
                
                print("Extracting colocated copernicus variable:",cop_var_name)
                cop_ds=colocated_data['copernicus'][i_obs]
                cop_dat=cop_ds["time"].values
                cop_lat=cop_ds["latitude"].values
                cop_lon=cop_ds["longitude"].values
                cop_var=cop_ds[cop_var_name].values
                
                # Adapt longitude domain when 180 is crossed (test this on 6903024 / 145A)
                if (np.max(cop_lon)-np.min(cop_lon)>180):
                    cop_lon[np.where(cop_lon<0)]=cop_lon[np.where(cop_lon<0)]+360
                    if obs_loni < 0:obs_loni=obs_loni+360
                    obs_lon=[obs_loni, obs_loni]
                
                from time import mktime
                from datetime import datetime
                obs_t = [mktime(datetime.strptime(np.datetime_as_string(i,unit='s'), "%Y-%m-%dT%H:%M:%S").timetuple())/86400 for i in obs_dat]
                t0 = np.min(obs_t)
                obs_t=obs_t-t0
                cop_t=[mktime(datetime.strptime(np.datetime_as_string(i,unit='s'), "%Y-%m-%dT%H:%M:%S").timetuple())/86400 for i in cop_dat]
                cop_t=cop_t-t0
                
                COP_t,COP_lat,COP_lon=np.meshgrid(cop_t,cop_lat,cop_lon)
                nt,ny,nx=COP_t.shape
                COP_t=np.squeeze(np.reshape(COP_t,(nt*ny*nx,1)))
                COP_lat=np.squeeze(np.reshape(COP_lat,(nt*ny*nx,1)))
                COP_lon=np.squeeze(np.reshape(COP_lon,(nt*ny*nx,1)))
                COP_var=np.squeeze(np.reshape(cop_var,(nt*ny*nx,1)))

                
                if not os.path.exists(outfig_dir):
                    os.mkdir(outfig_dir)
                # TO BE ADAPTED WHEN NOT ARGO
                fig_name_prefix=outfig_dir + dataset_id + "_" + wmo + "_" + "{:03d}".format(obs_cyci) + obs_diri
                
                fig = plt.figure()
                ax = fig.add_subplot(projection='3d')
                if (len(i_notnan[0])==0) | (cop_var_name !='CHL'):
                    scat1=ax.scatter(obs_lon, obs_lat, obs_t,c='black',s=30)
                else:
                    scat1=ax.scatter(obs_lon, obs_lat, obs_t, c=obs_chl, cmap='jet', vmin=np.nanmin(COP_var), vmax=np.nanmax(COP_var),s=30)
                #scat2=ax.scatter(COP_lon, COP_lat, COP_t, c=COP_var, cmap='jet', vmin=cop_var_MIN, vmax=cop_var_MAX,s=4)
                scat2=ax.scatter(COP_lon, COP_lat, COP_t, c=COP_var, cmap='jet', vmin=np.nanmin(COP_var), vmax=np.nanmax(COP_var),s=4)
                plt.colorbar(scat2,pad=0.15,label="copernicus " + cop_var_name + "[" + cop_var_UNIT + "]")
                plt.title("Argo float wmo id " + wmo + " cycle " + "{:d}".format(obs_cyci) + obs_diri + \
                          "\n{0:s} , {1:.3f}°N, {2:.3f}°E, at {3:.1f}dbar \n In-situ CHLA = {4:.2f} [{5:s}]".format(np.datetime_as_string(obs_dati,unit='s'),
                          obs_lati,obs_loni,obs_prei_notnan[0],obs_chli_notnan[0],obs_var_UNIT))
                
                
                ax.set_xlabel('lon')
                ax.set_xlim(np.min(cop_lon),np.max(cop_lon))
                ax.set_xticks(np.linspace(np.min(cop_lon),np.max(cop_lon),10))
                ax.set_ylabel('lat')
                
                ax.set_ylim(np.min(cop_lat),np.max(cop_lat))
                ax.set_yticks(np.linspace(np.min(cop_lat),np.max(cop_lat),10))
                
                ax.set_zlabel('date [delta days from obs]',rotation=90)
                ax.set_zlim(np.min(cop_t),np.max(cop_t))
                ax.set_zticks(np.linspace(np.min(cop_t),np.max(cop_t),10))
                
                ax.tick_params(axis='x',labelsize=6)
                ax.tick_params(axis='y',labelsize=6)
                ax.tick_params(axis='z',labelsize=6)

                print("The plot will be saved in " + fig_name_prefix + ".png")
                # plt.savefig
                plt.savefig(fig_name_prefix + ".png", dpi=200)
                pickle.dump(fig, open(fig_name_prefix + ".pkl", 'wb'))
                
                #plt.show()
                plt.close()
                
                # To see the figures in interactive mode:
                # figx = pickle.load(open(fig_name_prefix + "pkl", 'rb'))
                # figx.show() # Show the figure, edit it, etc.!
                
            status=flush_status_in_file(status_file,status,'stat_step_8_' + dataset_id,ST_completed)
        else:
            if (8 in steps_2_run) and (access_type == 'ARGO_INDEX'):
                print("ERROR: DO NOT PLOT FOR THE WHOLE INDEX, please")
                status=flush_status_in_file(status_file,status,'stat_step_8_' + dataset_id,ST_error)
                status=flush_status_in_file(status_file,status,'stat_step_8',ST_error)
            else:
                print(" ... skipped")
                status=flush_status_in_file(status_file,status,'stat_step_8',ST_skipped)
                status=flush_status_in_file(status_file,status,'stat_step_8_' + dataset_id,ST_skipped)
                
    if (7 in steps_2_run) & (access_type != 'ARGO_INDEX'):
        status=flush_status_in_file(status_file,status,'stat_step_7',ST_completed)
    if (8 in steps_2_run) & (access_type != 'ARGO_INDEX'):
        status=flush_status_in_file(status_file,status,'stat_step_8',ST_completed)
        
       
        
