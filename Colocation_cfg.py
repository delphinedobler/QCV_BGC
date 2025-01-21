# Colocation parameterisation file

# TODO: use a config file for the following parameters
# choose the input depending on your needs (the output can be tuned
access_type='ARGO_INDEX' # 'ARGO_DIRECT' or 'ARGO_CERBERE' or 'ARGO_INDEX' for the moment. this parameter will be used afterwards for plugging cerberized data
insitu_data_dir="C:/Users/ddobler/Documents/08_DD_scripts/09_FAIR-EASE/insitu-data/"
argo_dir=insitu_data_dir + "argo/"
cerbere_dir=insitu_data_dir + "cerbere/"

#wmo='6901578' # long journey float
wmo='6903024' # crosses 180 line (cycles 139 to 145 are on the West side of the line, the others on the East side)
workflow_name='chl'


copernicus_method='subset' # 'lazy' or 'subset' : I kept both, can be tuned
indexation_method='sel' # 'sel' or 'isel' or 'index' (in case of lazy access)
record_format='NetCDF' # 'values' or 'NetCDF' or 'computation': either data are get (.values) or locally saved in a NetCDF file. Used for performance assessments
verbose=False # the copernicus library can not yet be turned into quiet mode (but this works for informative prints)


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

clear_copernicus_resolution_cache=False
clear_copernicus_downloaded_data_index_cache=True
cache_dir="cache_files"
cache_copernicus_downloaded_data_index = cache_dir + "/cache_downloaded_data_index.csv"


log_file_cop='logs/log_copernicus_dl_performance.csv'
log_file_col_1="logs/log_get_data_to_colocate_obs_list.csv"
log_file_col_2="logs/log_get_data_to_colocate_goodobs_with_bbox_list.csv"
log_file_grp="logs/log_create_obs_group_function.csv"
location="office"
outfile_dir="copernicus-data/worflow_{0:s}_xn_{1:03d}_yn_{2:03d}_tn_{3:03d}".format(workflow_name,gp_crit['gp_max_x_n'],
                                                                                                gp_crit['gp_max_x_n'],
                                                                                                gp_crit['gp_max_t_n'])