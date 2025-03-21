#!/bin/bash

# This is the docker entrypoint that passes execution to a non root user
# -----------------------------------------------------------------------------------------------------------------


# -- if not launched by root.. exit
if [ "$(id -u)" != 0  ]; then
    echo "non production mode ..."
    bash
fi

# -- functions

# pretty indent for echoing in this script
echo_log() {
    local indent_level=$1
    local content=$2
    local indent=""

    # Create indentation string based on indent level
    for ((i=0; i<indent_level; i++)); do
        indent+="    "
    done

    # Print the content with the indent
    echo "= ${indent}${content}"
    echo "= ${indent}${content}" >> $RUN_LOG

}

env_defaults() {
    USER_ID=${FORCE_USER_ID:-1000}
    GROUP_ID=${FORCE_GROUP_ID:-1000}
    DOCKER_ENTRYPOINT_BYPASS=${DOCKER_ENTRYPOINT_BYPASS:-0}
    FORCE_ENTRYPOINT_CHOWN=${FORCE_ENTRYPOINT_CHOWN:-0}
}

entrypoint_actions() {

    echo_log 0 "----------------------------------"
    echo_log 0 ""
    echo_log 0 "Service account is  serviceuser : UID=$USER_ID GID=$GROUP_ID"
    echo_log 0 ""


    USER_NOT_DEFAULT=0

    # change user id
    if [ $GROUP_ID != 1000 ]; then
        groupmod -g ${GROUP_ID} -o servicegroup >> $RUN_LOG
        echo_log 2 "GID custom : $GROUP_ID"
        USER_NOT_DEFAULT=1
    else
        echo_log 2 "GID default : $GROUP_ID"
    fi

    # change user group
    if [ $USER_ID != 1000 ]; then
        usermod -u ${USER_ID} -o serviceuser >> $RUN_LOG
        echo_log 2 "UID custom : $USER_ID"
        USER_NOT_DEFAULT=1
    else
        echo_log 2 "UID default : $GROUP_ID"
    fi

    # change main directories owners / permissions
    if [ "$USER_NOT_DEFAULT" -eq 1 ]; then
    
        ls -la /home/serviceuser
        du -hs /home/serviceuser
        if [ "$FORCE_ENTRYPOINT_CHOWN" -eq 1 ]; then
            echo_log 2 "will change owner of : /home/serviceuser"
            chown -R serviceuser:servicegroup /home/serviceuser
            echo_log 2 "directory owner updated : /home/serviceuser"
        fi

        ls -la /runtime
        du -hs /runtime
        if [ "$FORCE_ENTRYPOINT_CHOWN" -eq 1 ]; then
            echo_log 2 "will change owner of : /home/serviceuser"
            chown -R serviceuser:servicegroup /runtime
            echo_log 2 "directory owner updated : /runtime"
        fi

    fi

    echo_log 2 "actions done."
}

# -- command line arguments

# # capture entrypoint context
# if [ $# -gt 0 ]; then
#     bypass_entrypoint=$1
# else
#     bypass_entrypoint=$DOCKER_ENTRYPOINT_BYPASS
# fi

# -----------------------------------------------------------------------------------------------------------------
# -- entrypoint code

# human readable time
START_TIME=$(date --utc +%Y%m%d%H%M%S)

# Load default Env values
env_defaults

# runtime log
RUN_LOGDIR=/runtime/run
RUN_LOG=$RUN_LOGDIR/pok_appstart.log


if [ "$DOCKER_ENTRYPOINT_BYPASS" -eq 0 ]; then
    echo_log 0 "==================================="
    echo_log 0 "= = = = = = = = = = = = = = = = = ="
    echo_log 0 "entrypoint account       : $(whoami)"
    # echo_log 0 "service user account     : serviceuser ( UID=$USER_ID  GID=$GROUP_ID )"
    echo_log 0 "lifecycle begins   : $START_TIME"
    echo_log 0 "entrypoint bypass  : $DOCKER_ENTRYPOINT_BYPASS" 
    echo_log 0 "entrypoint chown   : $FORCE_ENTRYPOINT_CHOWN" 
    echo_log 0 "start option       : $APP_START_OPTION"
    echo_log 0 "exec mode          : $APP_EXEC_MODE" 
    echo_log 0 ""
    echo_log 0 ".. executing actions"

    entrypoint_actions
    exit_code=$?
    echo_log 0 "actions done : exit code = $exit_code"
    
    echo_log 0 ""

    # echo_log 2 "passing command : [$@]"
    # exec /usr/sbin/gosu serviceuser "$@"
   
    echo_log 0 "= = = = = = = = = = = = = = = = = ="
    echo_log 0 "starting application with : gosu serviceuser /app/launcher/start-app.sh"

    exec /usr/sbin/gosu serviceuser /app/launchers/start-app.sh
    echo_log 0 ""

    END_TIME=$(date --utc +%Y%m%d%H%M%S)
    echo_log 0 "lifecycle end : $END_TIME"
fi



