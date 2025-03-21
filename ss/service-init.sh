#!/bin/bash
SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
cd $SCRIPT_DIR

echo "======================================="
echo "= service initialization               "
echo "= -------------------------------------"


NOW=$(date --utc +%Y%m%d%H%M%S)



# - - functions
source ./tool_functions.sh


# - - start options
if [ $# -gt 0 ]; then
    if [ $1 == "-f" ]; then
        OPTION="-f"
        echo_indent 2 "Force removing previous ../.env"
        mv ../.env ../.env.old.$NOW
        echo_indent 2 ""
    fi
else
    echo ..
fi

# Generate services .env
services=("postgis" "cache-storage" "mongodb" "proxy")
for service in "${services[@]}"; do
  echo "Generating .env for $service..."
  if ! (cd ../external-services/$service && ./service_init.sh $OPTION); then
    echo "Error: Failed to generate env for $service."
    exit 1
  fi
done



# - - Hardware Architecture

echo_indent 4 "INFO: Docker platform = $DOCKER_PLATFORM"
echo_indent 2 ""



# - - Regenerate template
./service-gen-env-template.sh


# ---------------------------------------------------------------------------------------------------

# - - - current user informations
MY_UID=$(id -u)
MY_GID=$(id -g)
POK_DEFAULT_DEVELOPER_SSH_KEY="$HOME/.ssh/pok-ed_$USER"

# - - - create ../.env file if not exists
echo_indent 2  "creating ../.env if not exists - - - -"

if [ -f "../.env" ]
then
    echo_indent 4  "../.env already exists"
    rm -f ../.env.example
else
    echo_indent 4  "new ../.env from ../.env.example"

    cp ../.env.example ../.env
    rm -f ../.env.example

    # -- POKAPOK USER ENV DEFAULTS
    # --- SSH KEY
    if [ -f "$POK_DEFAULT_DEVELOPER_SSH_KEY" ]
    then
        echo_indent 4  "found Pokapok standard user ssh key : $POK_DEFAULT_DEVELOPER_SSH_KEY"
        universal_sed DEVELOPER_SSH_KEY "$POK_DEFAULT_DEVELOPER_SSH_KEY" ../.env
    else
        echo_indent 4  "no default public key found."
    fi

    # --- USER NAME
    if [[ -z "${POK_MY_FULL_NAME}" ]]
    then
        echo_indent 4  "no default user name found."
    else
        echo_indent 4  "found Pokapok standard user name : $POK_MY_FULL_NAME"
        SAFE_USER_NAME=`quote_spaces "$POK_MY_FULL_NAME"`
        universal_sed DEVELOPER_GIT_USER_NAME "$SAFE_USER_NAME" ../.env
    fi

    # --- USER EMAIL
    if [[ -z "${POK_MY_EMAIL}" ]]
    then
        echo_indent 4  "no default user email found."
    else
        echo_indent 4  "found Pokapok standard user email : $POK_MY_EMAIL"
        universal_sed DEVELOPER_GIT_USER_EMAIL "$POK_MY_EMAIL" ../.env
    fi


fi

echo_indent 2 ""


# - - - DOCKER_PLATFORM ,  LOCAL USER & GROUP

universal_sed DOCKER_PLATFORM "$DOCKER_PLATFORM" ../.env
universal_sed FORCE_USER_ID "$MY_UID" ../.env
universal_sed FORCE_GROUP_ID "$MY_GID" ../.env


# - - - Create SSH dev-secrets folder and copy developer key inside
echo_indent 2  "creating ./dev-secrets-ssh folder - - - -"

# - - - get the DEVELOPER_SSH_KEY from env and expands to full path
DEVELOPER_SSH_KEY=$(eval echo $(cat ../.env | grep DEVELOPER_SSH_KEY | awk -F"=" '{print $2}') )

# - - - copy key in a local folder for developper volume mount
mkdir -p ../dev-secrets-ssh
cp ../system-configs/SSH/config ../dev-secrets-ssh

echo_indent 2 ""


echo_indent 2  "creating ./dev-secrets folder - - - -"

mkdir -p ../dev-secrets
touch ../dev-secrets/dummy

# - - - - copy private key
if  [ -f "${DEVELOPER_SSH_KEY}" ]
then
    cp ${DEVELOPER_SSH_KEY} ../dev-secrets-ssh/git-ssh-key
else
    echo_indent 4  "private key not found,  it should have been : ${DEVELOPER_SSH_KEY}"
    echo_indent 6  "-> fill in the .env  DEVELOPER_SSH_KEY  information and execute ./service-init.sh again"
fi

# - - - - copy public key
if  [ -f "${DEVELOPER_SSH_KEY}.pub" ]
then
    cp ${DEVELOPER_SSH_KEY}.pub ../dev-secrets-ssh/git-ssh-key.pub
else
    echo_indent 4  "public key not found , it should have been : ${DEVELOPER_SSH_KEY}.pub"
fi


echo_indent 2  "creating ./dev-mount folder for developement volumes mounting - - - -"
echo_indent 2  ""


# - - - - folder for publishing embedded resources on git repository
mkdir -p ../embedded-resources

# - - - - folder for publishing development documentation on git repository
#            !! Note : this is NOT sphinx or roxygen folder : only human generated files
mkdir -p ../app-dev-deploy-doc


# - - - - other useful folders for developer's debug
mkdir -p ../dev-mount-config
mkdir -p ../dev-mount-data
mkdir -p ../dev-mount-log
mkdir -p ../dev-mount-run
mkdir -p ../dev-external-libs



# - - - - secure local folder permissions
chmod 700 ../dev-secrets-ssh
chmod 600 ../dev-secrets-ssh/*
chmod 700 ../dev-secrets
chmod 600 ../dev-secrets/*



echo_indent 3  "!! Warning : Please copy all other credential files in ./dev-secrets"
echo_indent 4  "note : if needed, keycloak config file must have the following name : ServiceClient.json"
echo_indent 4  "note : Ensure that all secret files are private : chmod 600 !"

# - - end of task
echo "="
echo "= service init done."
echo "======================================="
echo
