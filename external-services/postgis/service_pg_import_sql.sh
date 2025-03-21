#!/bin/bash
# Pokapok PG Import Dockerstacks - 2022-01


echo ==========================================================================
echo Pokapok PG Import Dockerstack $0
echo --------------------------------------------------------------------------
echo  usage : $0  targetBDDName  dumpFile.sql
echo --------------------------------------------------------------------------

# -------------------------------------------------------
# ---- default variables

SOURCESTACK=${PWD##*/} 
PGHOST=localhost
USERNAME=`awk -F"=" '$1 ~ /^POSTGRES_USER/ {print $2}' .env`
USERPWD=`awk -F"=" '$1 ~ /^POSTGRES_PASSWORD/ {print $2}' .env`
IMPORT_DIR=`awk -F"=" '$1 ~ /^PROJECT_SQL_IMPORTS/ {print $2}' .env`
PGPORT=`awk -F"=" '$1 ~ /^POSTGRES_PORT/ {print $2}' .env`
LOCAL_IMPORTDIR=/imports

# Human Readable (HR) start time
HRTIME=$(date --utc +%Y-%m-%d_%H-%M-%S)

echo --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
echo IMPORT DIR content :
ls -lah $IMPORT_DIR
echo --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 



if [ -z $1 ];
then
    echo ... no parameters given
    exit 1
fi

if [ -z $2 ];
then
    echo ... 2 parameters are required
    exit 1
fi



DUMPFILE=$IMPORT_DIR/$1
DATABASE=$2

echo IMPORTING FROM : $IMPORT_DIR 
echo IN DATABASE    : $DATABASE

# echo ---
# todo : option create database
#echo "docker-compose exec -u postgres postGIS psql -U postgres -c \"CREATE DATABASE $targetBDD\""
#docker-compose exec -u postgres postGIS psql -U postgres -c "CREATE DATABASE \"$targetBDD\""

echo ---
docker-compose exec -u postgres postGIS pg_restore --dbname $DATABASE --no-privileges --role=$USERNAME  $LOCAL_IMPORTDIR/$1
echo - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - end
