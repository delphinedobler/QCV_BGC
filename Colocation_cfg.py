# Colocation parameterisation file


# I - input data selection
##########################
access_type='ARGO_DIRECT' # 'ARGO_DIRECT' or 'ARGO_CERBERE' or 'ARGO_INDEX' for the moment. this parameter will be used afterwards for plugging cerberized data
insitu_data_dir="C:/Users/ddobler/Documents/08_DD_scripts/09_FAIR-EASE/insitu-data/"
argo_dir=insitu_data_dir + "argo/"
cerbere_dir=insitu_data_dir + "cerbere/"
#wmo='6901578' # long journey float
wmo='6903024' # crosses 180 line (cycles 139 to 145 are on the West side of the line, the others on the East side)


# II - colocation parameterization
##################################
workflow_name='chl'

# Depending on your capacity, tune the grouping options
gp_crit={}
gp_crit['gp_max_x_n']=50#15#25#50#100 # i.e. within gp_max_x_n*reso_lon_deg, e.g. 200*0.04 = 8 deg
gp_crit['gp_max_y_n']=50#15#25#50#100
gp_crit['gp_max_t_n']=100#30#50#100#200

# defining the box around the observation
delta_px={}
delta_px['x']=5
delta_px['y']=5
delta_px['t']=5

copernicus_method='subset' # 'lazy' or 'subset' : I kept both, can be tuned
indexation_method='sel' # 'sel' or 'isel' or 'index' (in case of lazy access)
record_format='NetCDF' # 'values' or 'NetCDF' or 'computation': either data are get (.values) or locally saved in a NetCDF file. Used for performance assessments

outfile_dir="copernicus-data/worflow_{0:s}_xn_{1:03d}_yn_{2:03d}_tn_{3:03d}".format(workflow_name,gp_crit['gp_max_x_n'],
                                                                                                gp_crit['gp_max_x_n'],
                                                                                                gp_crit['gp_max_t_n'])
outfig_dir="colocated_plots/" + wmo + "/"
                                            


# III - paralellisation option
##############################
parallelisation = 'no' # 'no' or 'mpProcess' or 'mpAsync' or 'dask'


# IV - Steps_to_run
##################################
steps_2_run=[7,8]
#Steps that can be run at demand:
#Step 1: get in-situ data (read existing in cache if not in steps_2_run)
#Step 4: get remaining data to colocate (else all data from step 1 are re-colocated with copernicus data)
#Step 5: create groups of in-situ data (medium-cubes definition) (else {cache_group_of_obs_prefix}_{dataset_id}.pkl is used)
#Step 7: extract mini-cube from downloaded data
#Step 8: a few plots - TO BE CODED

# V - cache definition
######################

cache_dir = "cache_files"
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
log_file_cop='logs/log_copernicus_dl_performance.csv'
#log_files prefixes (will be suffixed by + dataset_id + ".csv" in the code)
log_file_col_1_prefix="logs/log_get_data_to_colocate_obs_list_"
log_file_col_2_prefix="logs/log_get_data_to_colocate_goodobs_with_bbox_list_"
log_file_grp_prefix="logs/log_create_obs_group_function_"
location="office"

# VII - standard output additionnal prints
##########################################
verbose=True 
# N.B.the copernicus library can not yet be turned into quiet mode (but this works for informative prints)
      
