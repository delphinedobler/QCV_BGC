#!/bin/bash
SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
cd $SCRIPT_DIR

# - - functions
source ../.env
source ./tool_functions.sh

# Ensure a default config deployment file exists
SRC_CONFIG_FILE=../env-templates/config-deployment.yaml
TARGET_CONFIG_FILE=../dev-mount-config/config-deployment.yaml

if [ -f $TARGET_CONFIG_FILE ]; then
   echo "File $TARGET_CONFIG_FILE exists."
   echo ".. if not wanted ,  remove it yourself."
else
   echo "Copy config deployment-file in dev-mount-config..."
   cp $SRC_CONFIG_FILE $TARGET_CONFIG_FILE
fi
echo "Don't forget to verify that values defined in /runtime/config/config-deployment.yaml are correct."


# get docker-compose file name
dcfile=$(get_dc_file)

# execute command
docker-compose -f ../$dcfile up -d
