# Colocation parameterisation file

# Tool repository:
tool_dir="C:/Users/ddobler/Documents/08_DD_scripts/09_FAIR-EASE/Colocation/"

# I - input data selection
##########################
access_type='ARGO_INDEX' # 'ARGO_DIRECT' or 'ARGO_CERBERE' or 'ARGO_INDEX' for the moment. this parameter will be used afterwards for plugging cerberized data
insitu_data_dir=tool_dir+"/insitu-data/"
argo_dir=insitu_data_dir + "argo/"
cerbere_dir=insitu_data_dir + "cerbere/"
#wmo='6901578' # long journey float
wmo='6903024' # crosses 180 line (cycles 139 to 145 are on the West side of the line, the others on the East side)


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

    #l_dataset=[dataset_chl,dataset_rrs,dataset_Kd]
    l_dataset=[dataset_Kd]
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
# delta_px['x']=5
# delta_px['y']=5
# delta_px['t']=5
delta_px['x']=6 # add one pixel when dealing with the entire index: this helps is 1e-4 digit is changed in location.
delta_px['y']=6
delta_px['t']=6

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
# igrp_2_colocate=[2187,2191,2224,2230,2232,2238,2239,2259,2262,2337,2561,2614,2617,2681,2697,2709,2712,2736,2783,
# 2818,2819,2859,2860,2861,2862,2863,2864,2865,2866,2867,2868,2883,2886,2890,2891,2909,2910,2914,2935,
# 2936,2941,2942,3039,3089,3103,3104,3105,3116,3125,3129,3135,3187,3252,3257,3292,3319,3336,3435,3436,
# 3441,3443,3448,3454,3458,3493,3500,3519,3520,3544,3730,3827,3829,3841,3917,3985,4023,4031,4093,4101,
# 4107,4110,4113,4145,4170,4171,4182,4185,4188,4189,4203,4221,4222,4223,4233,4292,4293,4296,4299,4353,
# 4403,4445,4446,4447,4528,4529,4534,4554,4620,4625,4626,4640,4682,4812,4813,4814,4815,4816,4817,4818,
# 4829,4866,4879,4885,4917,4919,4920,4930,4937,4962,4966,4974,4987,5003,5011,5014,5020,5021,5035,5061,
# 5062,5108,5194,5199,5203,5250,5268,5279,5286,5402,5406,5408,5409,5411,5413,5464,5465,5466,5476,5477,
# 5562,5570,5588,5604,5605,5649,5671,5695,5757,5777,5780,5782,5820,5825,5828,5829,5830,5836,5863,5905,
# 5907,5937,5938,5939,5940,5942,5949,5982,5983,5984,6000,6003,6004,6039,6096,6097,6098,6099,6117,6140,
# 6145,6214,6257,6279,6301,6318,6319,6320,6326,6328,6332,6334,6348,6350,6360,6362,6365,6406,6407,6466,
# 6530,6531,6532,6533,6540,6616,6618,6623,6633,6723,6802,6803,6955,7016,7056,7077,7079,7080,7081,7082,
# 7110,7112,7115,7151,7154,7155,7161,7180,7204,7206,7259,7336,7339,7340,7342,7346,7350,7367,7371,7372,
# 7373,7375,7376,7379,7380,7382,7383,7384,7386,7399,7420,7423,7435,7442,7465,7478,7519,7525,7554,7584,
# 7586,7594,7595,7601,7634,7637,7649,7650,7653,7654,7656,7677,7701,7702,7704,7705,7706,7708,7738,7739,
# 7742,7755,7805,7810,7811,7818,7821,7825,7826,7829,7842,7846,7851,7852,7855,7865,7875,7921,7934,7937,
# 7938,7939,7941,7972,7976,7978,7979,7980,7981,7982,7983,7984,7985,7989,7990,7992,7993,7994,7995,8003,
# 8004,8008,8010,8011,8012,8014,8015,8016,8035,8037,8067,8072,8090,8091,8096,8104,8105,8109,8110,8141,
# 8147,8171,8172,8173,8174,8175,8181,8204,8205,8210,8212,8215,8219,8265,8276,8282,8290,8291,8292,8297,
# 8332,8347,8448,8491,8492,8534,8535,8536,8537,8538,8549,8557,8562,8567,8570,8604,8606,8608,8609,8610,
# 8611,8614,8621,8622,8623,8631,8632,8635,8636,8637,8641,8642,8644,8645,8646,8656,8657,8660,8663,8679,
# 8682,8693,8709,8711,8716,8722,8724,8726,8727,8741,8803,8804,8806,8807,8809,8815,8817,8818,8827,8828,
# 8832,8834,8839,8841,8844,8939,8940,8941,8949,8950,8951,8956,9032,9033,9081,9082,9084,9087,9096,9147,
# 9148,9171,9175,9176,9177,9178,9179,9184,9188,9189,9270,9300,9385,9387,9388,9393,9394,9395,9420,9421,
# 9435,9441,9444,9465,9466,9494,9495,9496,9497,9498,9554,9555,9559,9616,9638,9639,9641,9671,9710,9745,
# 9746,9748,9783,9784,9785,9786,9811,9938,9947,9962,9992,10000,10016,10018,10032,10108,10109,10110,10113,10115,
# 10130,10133,10273,10355,10367,10368,10376,10379,10401,10414,10417,10418,10420,10612,10613,10614,10615,10617,10618,10647,
# 10648,10649,10650,10651,10655,10656,10660,10666,10668,10673,10674,10675,10676,10677,10678,10679,10681,10682,10685,10687,
# 10689,10690,10691,10693,10734,10737,10738,10739,10744,10745,10746,10749,10753,10754,10759,10760,10761,10762,10763,10771,
# 10772,10783,10852,10855,10856,10861,10867,10868,10875,10901,10903,10904,10964,10965,10966,10968,10969,10971,10973,10987,
# 10989,10990,10991,11005,11010,11011,11012,11014,11016,11025,11026,11027,11028,11034,11076,11078,11079,11080,11081,11082,
# 11083,11086,11093,11098,11101,11103,11104,11105,11106,11107,11109,11125,11126,11128,11132,11133,11134,11144,11158,11159,
# 11160,11161,11162,11173,11174,11175,11178,11179,11181,11182,11190,11192,11203,11206,11232,11234,11276,11278,11279,11303,
# 11311,11317,11407,11445,11467,11532,11542,11567,11573,11574,11594,11638]

# IV - Steps_to_run
##################################
steps_2_run=[2,5]
#Steps that can be run or not at demand:
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
location="office"

# VII - standard output additionnal prints
##########################################
verbose=False 
# N.B.the copernicus library can not yet be turned into quiet mode (but this works for informative prints)
      
