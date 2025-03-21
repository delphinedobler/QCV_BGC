# Colocation parameterisation file

# Tool repository:
tool_dir="C:/Users/ddobler/Documents/08_DD_scripts/09_FAIR-EASE/Colocation/"

# I - input data selection
##########################
access_type='ARGO_DIRECT' # 'ARGO_DIRECT' or 'ARGO_CERBERE' or 'ARGO_INDEX' for the moment. this parameter will be used afterwards for plugging cerberized data
insitu_data_dir=tool_dir+"/insitu-data/"
argo_dir=insitu_data_dir + "argo/"
cerbere_dir=insitu_data_dir + "cerbere/"
wmo='6901578' # long journey float
#wmo='6903024' # crosses 180 line (cycles 139 to 145 are on the West side of the line, the others on the East side)


# II - colocation parameterization
##################################
workflow_name='chl'

if workflow_name == "chl":
    dataset_rrs='cmems_obs-oc_glo_bgc-reflectance_my_l3-multi-4km_P1D'
    rrs_var=['RRS412','RRS443','RRS490','RRS555','RRS670']
    dataset_chl='cmems_obs-oc_glo_bgc-plankton_my_l3-multi-4km_P1D'
    chl_var=['CHL']
    dataset_Kd='cmems_obs-oc_glo_bgc-transp_my_l3-multi-4km_P1D'
    Kd_var=['KD490']
    dataset_rrs_short_name='cms_reflectance'
    dataset_chl_short_name='cms_plankton'
    dataset_Kd_short_name='cms_transp'

    l_dataset=[dataset_chl,dataset_rrs,dataset_Kd]
    #l_dataset=[dataset_Kd]
    l_dataset_short_name={}
    l_dataset_short_name[dataset_chl]=dataset_chl_short_name
    l_dataset_short_name[dataset_rrs]=dataset_rrs_short_name
    l_dataset_short_name[dataset_Kd]=dataset_Kd_short_name
    d_dataset_var={}
    d_dataset_var[dataset_chl]=chl_var
    d_dataset_var[dataset_rrs]=rrs_var
    d_dataset_var[dataset_Kd]=Kd_var


# Depending on your capacity, tune the grouping options
gp_crit={}
gp_crit['gp_max_x_n']=50#15#25#50#100 # i.e. within gp_max_x_n*reso_lon_deg, e.g. 200*0.04 = 8 deg
gp_crit['gp_max_y_n']=50#15#25#50#100
gp_crit['gp_max_t_n']=100#30#50#100#200
# two limits to prevent downloading when grouping is in error (e.g. step 2 is rerun but not step 5)
spatial_extension_square_deg_limit=30 # in square_degrees
temporal_extension_days_limit=215 # in days

# defining the box around the observation
delta_px={}
delta_px['x']=5
delta_px['y']=5
delta_px['t']=5
# delta_px['x']=6 # add one pixel when dealing with the entire index: this helps is 1e-4 digit is changed in location.
# delta_px['y']=6
# delta_px['t']=6

copernicus_method='subset' # 'lazy' or 'subset' : I kept both, can be tuned
indexation_method='sel' # 'sel' or 'isel' or 'index' (in case of lazy access)
record_format='NetCDF' # 'values' or 'NetCDF' or 'computation': either data are get (.values) or locally saved in a NetCDF file. Used for performance assessments

outdir_cop=tool_dir + "copernicus-data/"
                                                                                               
outdir_col_plots = tool_dir + "colocated_plots/"
                                            


# III - paralellisation option
##############################
parallelisation = 'mpAsync' # 'no' or 'mpProcess' or 'mpAsync'
grp_deb=0  # for debug or testing: use grp_deb and grp_end to specify which groups you want to download.
grp_end=-1 # use grp_deb=0 and grp_end=-1 to perform the complete set.
igrp_2_colocate=[]# set igrp_2_colocate=[] to use grp_deb and grp_end feature
# igrp_2_colocate=[2187,2191,2224,2230,2232,2238,2239,2259,2262,2337,2561,2614,2617,2681,2697,2709,2712,2736]

# IV - Steps_to_run
##################################
steps_2_run=[2,4,5,7,8]
#Steps that can be run or not at demand. It is advised to rerun all steps for consistency.
#Step selection should be triggered only for specific reasons (debugs, techniques, etc.)
#Step 2: get in-situ data (read existing in cache if not in steps_2_run) # if you run step 2, step 5 will be run as well
#Step 4: get remaining data to colocate (else all data from step 1 are re-colocated with copernicus data) # if you run step 4, make sure you run step 5 as well.
#Step 5: create groups of in-situ data (medium-cubes definition) (else {cache_group_of_obs_prefix}_{dataset_id}.pkl is used)
#Step 7: extract mini-cube from downloaded data
#Step 8: a few plots
#Steps that are run anyway:
#Step 1: read configuration file
#Step 3: get copernicus spatio-temporal resolution and limits (to download again, delete cache (clear_cache_copernicus_resolution=True), otherwise, read from cache
#Step 6: download copernicus data for missing observations

# V - cache definition
######################

cache_dir = tool_dir + "cache_files/"
cache_copernicus_resolution_file = cache_dir + "/cache_datasets_for_"+workflow_name+"_workflow_spatial_resolution.csv"
cache_copernicus_downloaded_data_index = cache_dir + "/cache_downloaded_data_index.csv"
#cache_files prefixes (will be suffixed by + dataset_id + ".pkl" in the code)
cache_group_of_obs_prefix = cache_dir + "/cache_group_of_obs_"
clear_cache_copernicus_resolution=False
clear_cache_copernicus_downloaded_data_index=False
clear_cache_group_of_obs=False

# VI - additionnal logs 
#######################

log4debug=True
log_dir = tool_dir + "logs/"
log_file_cop = log_dir + "log_copernicus_dl_performance.csv"
#log_files prefixes (will be suffixed by + dataset_id + ".csv" in the code)
log_file_col_1_prefix = log_dir + "log_get_data_to_colocate_obs_list_"
log_file_col_2_prefix = log_dir + "log_get_data_to_colocate_goodobs_with_bbox_list_"
log_file_grp_prefix = log_dir + "log_create_obs_group_function_"
status_file = log_dir + "status.txt"
location="office"

# VII - standard output additionnal prints
##########################################
verbose=False 
# N.B.the copernicus library can not yet be turned into quiet mode (but this works for informative prints)
      
