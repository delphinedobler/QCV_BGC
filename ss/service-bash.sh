#!/bin/bash

SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
cd $SCRIPT_DIR

# - - functions
source ../.env
source ./tool_functions.sh

# get docker-compose file name
dcfile=$(get_dc_file)

# execute command
docker-compose -f ../$dcfile exec -u serviceuser pok_app /bin/bash
