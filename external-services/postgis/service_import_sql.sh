#!/bin/bash
# Pokapok PG Import Dockerstacks - 2022-01


echo ==========================================================================
echo Pokapok PG Import Dockerstack $0

# -------------------------------------------------------
# ---- default variables

SOURCESTACK=${PWD##*/} 
PGHOST=localhost
USERNAME=`awk -F"=" '$1 ~ /^POSTGRES_USER/ {print $2}' .env`
USERPWD=`awk -F"=" '$1 ~ /^POSTGRES_PASSWORD/ {print $2}' .env`
PGPORT=`awk -F"=" '$1 ~ /^POSTGRES_PORT/ {print $2}' .env`
LOCAL_IMPORTDIR=/imports

# Human Readable (HR) start time
HRTIME=$(date --utc +%Y-%m-%d_%H-%M-%S)

# Source file
SOURCEFILE=./exports/${SOURCESTACK}/${SOURCESTACK}_@${DATABASE}@_${HRTIME}


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


DATABASE=$2
SQLDUMP=$1



echo IMPORTING DATABASE = $DATABASE
echo from $SQLDUMP

# echo ---
# echo "docker-compose exec -u postgres postGIS psql -U postgres -c \"CREATE DATABASE $targetBDD\""
# docker-compose exec -u postgres postGIS psql -U postgres -c "CREATE DATABASE \"$targetBDD\""

echo ---
docker-compose exec -u postgres postGIS psql $DATABASE -f $LOCAL_IMPORTDIR/$SQLDUMP
echo - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - end
