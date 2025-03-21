#!/bin/bash

SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
cd $SCRIPT_DIR

# - - functions
source ../.env
source ./tool_functions.sh

# get docker-compose file name
dcfile=$(get_dc_file)


dockercmd_demo() {
    # launch the demo
    mkdir -p ./dev-demo-data/log
    mkdir -p ./dev-demo-data/out

    cmd="docker run \
        -v ./dev-demo-data/log:/runtime/log \
        -v ./dev-demo-data/out:/runtime/data-out \
        -e CSTM_NAME=exampleRun \
        -e APP_START_OPTION=DEMO \
        -e FORCE_USER_ID=$(id -u) \
        -e FORCE_GROUP_ID=$(id -g) \
        $DOCKER_IMAGENAME_V"

    echo executing this command :
    echo $cmd
    echo

    $cmd
}



# redirect to chosen start option
case $1 in
"PYTHON_SHINY")
    echo docker command not yet defined
    ;;
"R_SHINY")
    echo docker command not yet defined
    ;;
"PYTHON_FASTAPI")
    echo docker command not yet defined
    ;;
"PYTHON_APP")
    echo docker command not yet defined
    ;;
"DEMO")
    dockercmd_demo
    ;;
"BASH")
    echo docker command not yet defined
    ;;
*)
    echo docker command not yet defined
    ;;
esac


