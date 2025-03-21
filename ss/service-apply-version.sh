#!/bin/bash
SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
cd $SCRIPT_DIR

API_VERSION=`cat ../APIVersion.tag | grep API_VERSION | awk -F"=" '{print $2}'`


YEAR=`echo ${API_VERSION} | awk -F"-" '{print $1}'`
MONTH=`echo ${API_VERSION} | awk -F"-" '{print $2}'`
DATAMODEL=`echo ${API_VERSION} | awk -F"-" '{print $3}'`
ENVIRONMENT=`echo ${API_VERSION} | awk -F"-" '{print $4}'`


echo YEAR        = $YEAR
echo MONTH       = $MONTH
echo DATAMODEL   = $DATAMODEL
echo ENVIRONMENT = $ENVIRONMENT


sed -i "s|DOCKER_TAG=.*|DOCKER_TAG=$API_VERSION|g" ../.env
sed -i "s|APP_PACKAGE_VERSION=.*|APP_PACKAGE_VERSION=$API_VERSION|g" ../.env
