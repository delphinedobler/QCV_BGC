#!/bin/bash
SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
cd $SCRIPT_DIR

############################################################
# clean project                                            #
############################################################


read -p "Warning : you will loose your .env config and your secrets. Confirm (y/n) " response
if [[ $response == "y" ]]; then
    rm ../.env
    rm -rf ../dev-mount*
    rm -rf ../dev-secrets*
    echo ../cleaning done.
else
    echo aborting
fi




