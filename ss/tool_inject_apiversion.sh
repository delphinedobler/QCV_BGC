#!/bin/bash

SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
cd $SCRIPT_DIR

# - - functions
source ../.env
source ./tool_functions.sh

# ..
NOWADAY=$(date --utc +%Y-%m-%d)

# - - Hardware Architecture
echo_indent 4 "INFO: Docker platform = $DOCKER_PLATFORM"
echo_indent 2 ""

# ================================================

echo_indent 0 ""
echo_indent 0 "-> injecting current API and GIT versions into flag files :"
echo_indent 0 ""




# API

#   source files
API_DATAMODEL_VERSION_FILE=../API_VERSION_DATAMODEL.txt
API_ENGINE_VERSION_FILE=../API_VERSION_ENGINE.txt

echo_indent 2 "API informations, source files :"
echo_indent 4 "Datamodel version : $API_DATAMODEL_VERSION_FILE"
echo_indent 4 "Engine version    : $API_ENGINE_VERSION_FILE"

#   retrieve current API version informations (should also be GIT branches : see branch naming pratices)
API_DATAMODEL_VERSION=$(cat ${API_DATAMODEL_VERSION_FILE} | grep API_DATAMODEL_VERSION | awk -F"=" '{print $2}')
API_ENGINE_VERSION=$(cat ${API_ENGINE_VERSION_FILE} | grep API_ENGINE_VERSION | awk -F"=" '{print $2}')

#   forge API_VERSION.tag
API_CURRENT_VERSION_TAG=${API_DATAMODEL_VERSION}-${API_ENGINE_VERSION}_${NOWADAY}

#   target file
API_VERSION_FILE=../CURRENT_VERSION_API.tag

#   update git engine version tag
universal_sed API_CURRENT_VERSION_TAG "$API_CURRENT_VERSION_TAG" $API_VERSION_FILE
echo_indent 2 "API informations, target file : ${API_VERSION_FILE}"
echo_indent 4 "-> API_CURRENT_VERSION_TAG = $API_CURRENT_VERSION_TAG"
echo_indent 0 ""


# BRANCH
#   retrieve current commit informations
GITBRANCH=`git rev-parse --abbrev-ref HEAD`
GITCOMMIT=`git rev-parse --short HEAD`
GIT_COMMIT_DATE=$(git show -s --format=%cd --date=short)
echo_indent 2 "GIT informations from current branch/commit"

#   forge git commit tag
GIT_CURRENT_VERSION_TAG=${GITBRANCH}-${GITCOMMIT}_${GIT_COMMIT_DATE}

#   target file
GIT_VERSION_FILE=../CURRENT_VERSION_GIT.tag

#   update git commit tag
universal_sed GIT_CURRENT_VERSION_TAG "$GIT_CURRENT_VERSION_TAG" $GIT_VERSION_FILE
echo_indent 2 "GIT informations, target file : ${GIT_VERSION_FILE}"
echo_indent 4 "-> GIT_CURRENT_VERSION_TAG = $GIT_CURRENT_VERSION_TAG"
echo_indent 0 ""










