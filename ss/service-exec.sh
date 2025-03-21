#!/bin/bash

source ../.env

# - - start options
if [ $# -gt 0 ]; then
    thecommand=$1
else
    thecommand=/bin/bash
fi

docker run -it \
    --user root \
    DOCKER_IMAGENAME_V $thecommand