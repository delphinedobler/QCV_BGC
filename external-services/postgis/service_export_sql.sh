#!/bin/bash
# Pokapok PGDUMP Dockerstacks - 2022-01


echo ==========================================================================
echo Pokapok PGDUMP Dockerstack $0

# -------------------------------------------------------
# ---- default variables

SOURCESTACK=${PWD##*/} 
PGHOST=localhost
USERNAME=`awk -F"=" '$1 ~ /^POSTGRES_USER/ {print $2}' .env`
USERPWD=`awk -F"=" '$1 ~ /^POSTGRES_PASSWORD/ {print $2}' .env`
EXPORT_DIR=`awk -F"=" '$1 ~ /^PROJECT_SQL_EXPORTS/ {print $2}' .env`

DATABASE=$1
PGPORT=`awk -F"=" '$1 ~ /^POSTGRES_PORT/ {print $2}' .env`
LOCAL_EXPORTDIR=/exports

# Human Readable (HR) start time
HRTIME=$(date --utc +%Y-%m-%d_%H-%M-%S)

# target file
TARGETFILE=${EXPORT_DIR}/${HOSTNAME}_${SOURCESTACK}_${DATABASE}_${HRTIME}.sql

echo --- job infos
CMD="pg_dump --no-owner -U ${USERNAME} ${DATABASE} > ${TARGETFILE}"
echo COMMAND = ${CMD}
echo note : add option -s   for schema only export
echo ---
echo --- --- --- 
echo EXPORTING DATABASE = ${DATABASE}
echo from ${SOURCESTACK}
echo by user : ${USERNAME}
echo --- --- ---
# docker-compose exec -u postgres postGIS pg_dump --create --dbname=postgresql://${USERNAME}:${USERPWD}@127.0.0.1:${PGPORT}/${DATABASE} > ${TARGETFILE}
docker-compose exec -u postgres postGIS pg_dump --no-owner --dbname=postgresql://${USERNAME}:${USERPWD}@127.0.0.1/${DATABASE} > ${TARGETFILE}
echo - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 