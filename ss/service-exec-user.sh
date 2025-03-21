#!/bin/bash

# - - start options
if [ $# -gt 0 ]; then
    thecommand=$1
else
    thecommand=/bin/bash
fi

docker run -it \
    --user root \
    -e FORCE_USER_ID=$(id -u) \
    -e FORCE_GROUP_ID=$(id -g) \
    DOCKER_IMAGENAME_V $thecommand