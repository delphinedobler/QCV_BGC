# Colocation tool

## Purpose
Colocation tool allows to:
- pre-download into cache files copernicus data colocated with a full batch of in-situ observations, grouped by spatio-temporal criteria (a.k.a. medium cubes)
- download remaining missing colocated copernicus data for a specific set of observations
- output colocated data as a {"Format still to be defined - either dataset or dataframe or netCDF"}, for each observation, using colocation limits criteria (a.k.a mini-cubes)
- display colocated data

## How to
To use this tool:
1) make sure you installed the required libraries: pip install -r requirements.txt
2) locally update you copernicus credentials (you only need it once):<br>
> import copernicusmarine <br> copernicusmarine.login()

3) update the configuration file: Colocation_cfg.py
4) run the tool:
python Colocation.py

## Configuration parameters
To be written

## inputs
To be written

## outputs
To be written
