#!/bin/bash

# -------------------------------------- Pokapok Boiletplate Component ------------------------------------------
# --- start-app.sh ---
# ---
# --- This start script is mandatory in this project boilerplate.
# --- Its name is static and it is launched by the docker image default command : entrypoint.sh
# --- It is aimed at selecting the application starting mode and options
# ---------------------------------------------------------------------------------------------------------------

# useful functions
echo_log() {
    local indent_level=$1
    local content=$2
    local indent=""
    local logfile=${APP_LOGFILE:-/runtime/log/default.log}

    # Create indentation string based on indent level
    for ((i=0; i<indent_level; i++)); do
        indent+="  "
    done

    # Print the content with the indent
    echo ". ${indent}${content}"
    echo ". ${indent}${content}" >> $logfile
}

# be sure to have bashrc context
source /home/serviceuser/.bashrc
# Application ENV vars custom




# Application ENV vars defaults
export NETWORK_BACK_INTERNAL_PORT=${NETWORK_BACK_INTERNAL_PORT:-8841}
export NETWORK_BACK_EXPOSED_PORT=${NETWORK_BACK_EXPOSED_PORT:-8841}
export NETWORK_FRONT_INTERNAL_PORT=${NETWORK_FRONT_INTERNAL_PORT:-8842}
export NETWORK_FRONT_EXPOSED_PORT=${NETWORK_FRONT_EXPOSED_PORT:-8842}
export APP_EXEC_MODE=${APP_EXEC_MODE:-PROD}
export APP_START_OPTION=${APP_START_OPTION:-PYTHON_APP}
export APP_PACKAGE_NAME=${APP_PACKAGE_NAME:-cerberecontrib}

# main folders
export APP_DATA_ROOT=/runtime
export APP_CODE_ROOT=/app
export SERVICE_ACCOUNT_HOME=/home/serviceuser

# execution log
export APP_LOGDIR=$APP_DATA_ROOT/log
export APP_LOG=$APP_LOGDIR/$APP_PACKAGE_NAME

# current date
export START_DATE=$(date --utc +%Y%m%d%H%M%S)


# capture interactive options
if [ $# -gt 0 ]; then
    start_option=$1
else
    start_option=$APP_START_OPTION
    echo_log 0 "  Chosen start option is : $start_option "
fi

if [ $# -gt 1 ]; then
    exec_mode=$2
else
    exec_mode=$APP_EXEC_MODE
    echo_log 0 "  Chosen exec mode is : $exec_mode "
fi

# apply DEV / PROD modes
if [[ $exec_mode=="DEV" ]]; then
    echo_log 0 "  DEV Mode actions : none defined"
else
    echo_log 0 "  PROD Mode actions : none defined"
fi
echo_log 0 " "


echo_log 0 "- - - - - - - - - - - - - - - - - - "
echo_log 0 "- -  start-app.sh - - - - - - - - - $(whoami)"
echo_log 0 "- - - - - - - - - - - - - - - - - - $START_DATE"
echo_log 0 " "
echo_log 0 " service account home is  : $SERVICE_ACCOUNT_HOME "
echo_log 0 " "
echo_log 0 "  starting option is      : $start_option "
echo_log 0 "  exec mode is            : $exec_mode "
echo_log 0 "  app source code path is : $APP_CODE_ROOT "
echo_log 0 "  app source code path is : $APP_CODE_ROOT "
echo_log 0 "  Read Write data vol is  : $APP_DATA_ROOT/data "
echo_log 0 "  Read only data vol is   : $APP_DATA_ROOT/resources "
echo_log 0 "  App Log vol is          : $APP_DATA_ROOT/logs "
echo_log 0 "  Run log vol is          : $APP_DATA_ROOT/run "
echo_log 0 " "
echo_log 0 "  front internal port is : $NETWORK_FRONT_INTERNAL_PORT "
echo_log 0 "  front external port is : $NETWORK_FRONT_EXPOSED_PORT "
echo_log 0 "  back internal port is : $NETWORK_BACK_INTERNAL_PORT "
echo_log 0 "  back external port is : $NETWORK_BACK_EXPOSED_PORT "
echo_log 0 " "echo_log 0 " "



# generic commands and paths
R_CMD="/opt/miniforge3/envs/env-develop/bin/Rscript"
PYTHON_CMD="/opt/miniforge3/envs/env-develop/bin/python"
RSHINY_CMD="/opt/miniforge3/envs/env-develop/bin/Rscript"
PYSHINY_CMD="/opt/miniforge3/envs/env-develop/bin/shiny run"
PYSHINY_APP_PATH="/home/serviceuser/service/app/src_python/exampleproject/main.py"
# FASTAPI_CMD="cd /home/serviceuser/service/app && gunicorn -c /home/serviceuser/service/app/src_python/example_fastapi/gunicorn_conf.py -w $(($(nproc) * 2)) -k uvicorn.workers.UvicornWorker"
FASTAPI_CMD="gunicorn -c /home/serviceuser/service/app/src_python/example_fastapi/gunicorn_conf.py \
    -w $(($(nproc) * 2)) -k uvicorn.workers.UvicornWorker \
    src_python.example_fastapi.app:app --bind 0.0.0.0:$NETWORK_BACK_INTERNAL_PORT --timeout 5000 \
    --log-level debug --access-logfile - --error-logfile -"

# redirect to chosen start option
case ${exec_mode} in
"DEV")
    # pyshiny
    PYSHINY_OPTIONS="--reload --host 0.0.0.0 --port $NETWORK_FRONT_INTERNAL_PORT"

    # rshiny
    RSHINY_AUTORELOAD=TRUE

    # fastapi
    FASTAPI_OPTIONS=""

    # python app

    # log
    APP_LOGFILE="$APP_LOG-$exec_mode"
    ;;
"PROD")
    # pyshiny
    PYSHINY_OPTIONS="--host 0.0.0.0 --port $NETWORK_FRONT_INTERNAL_PORT"

    # rshiny
    RSHINY_AUTORELOAD=FALSE

    # fastapi
    FASTAPI_OPTIONS=""

    # python app

    # log
    APP_LOGFILE="$APP_LOG-$exec_mode"

    ;;
*)
    echo_log 0 "  no valid exec mode passed ( APP_EXEC_MODE=DEV or APP_EXEC_MODE=PROD)... fatal"
    exit 1
    ;;
esac

# redirect to chosen start option
case ${start_option} in
"PYTHON_SHINY")
    APP_LOGFILE=$APP_LOGFILE.py_shiny.log
    $PYSHINY_CMD $PYSHINY_OPTIONS $PYSHINY_APP_PATH > $APP_LOGFILE 2>&1 &
    SHINY_PID=$!
    echo "Shiny started with PID: $SHINY_PID"
    echo Shiny internal port is $NETWORK_FRONT_INTERNAL_PORT 
    
    # Wait only in DEV mode
    if [ "$exec_mode" = "DEV" ]; then
        wait
    fi
    ;;
"R_SHINY")
    APP_LOGFILE=$APP_LOGFILE.r_shiny.log
    $RSHINY_CMD ./src_R/R_shiny/Shiny_example.R > $APP_LOGFILE
    ;;
"PYTHON_FASTAPI")

    # Signal handling
    terminate() {
        echo "Terminating processes..."
        # Stop Gunicorn
        if [ -n "$GUNICORN_PID" ]; then
            echo "Stopping Gunicorn (PID: $GUNICORN_PID)..."
            kill -SIGTERM "$GUNICORN_PID" 2>/dev/null
            wait "$GUNICORN_PID"
        fi

        # # Stop Python routine
        # if [ -n "$ROUTINE_PID" ]; then
        #     echo "Stopping Python routine (PID: $ROUTINE_PID)..."
        #     kill -SIGTERM "$ROUTINE_PID" 2>/dev/null
        #     wait "$ROUTINE_PID"
        # fi

        echo "All processes terminated."
        exit 0
    }

    # Trap signals
    trap terminate SIGINT SIGTERM

    APP_FASTAPI_LOGFILE=$APP_LOGFILE.py_fastapi.log

    # $FASTAPI_CMD src_python.example_fastapi.app:app --bind 0.0.0.0:$NETWORK_FRONT_INTERNAL_PORT --timeout 5000 --log-level debug --access-logfile - --error-logfile > $APP_LOGFILE &

    cd /home/serviceuser/service/app
    echo "Starting Gunicorn..."
    $FASTAPI_CMD > "$APP_FASTAPI_LOGFILE" 2>&1 &
    GUNICORN_PID=$!
    echo "Gunicorn started with PID: $GUNICORN_PID"
    echo fast api internal port is $NETWORK_BACK_INTERNAL_PORT 

    # cd /home/serviceuser/service/app && uvicorn src_python.boilerplate_fs_api.cache.routines:app > $BACKGROUND_ROUTINE_LOG &

    # Start Python routine in background (sync cache with database)
    # ROUTINE_CMD="/opt/miniforge3/envs/env-develop/bin/python /home/serviceuser/service/app/src_python/example_fastapi/routines.py"
    # echo "Starting Python sync routine..."
    # $ROUTINE_CMD > "$APP_LOGFILE.py_routine.log" 2>&1 &
    # ROUTINE_PID=$!
    # echo "Python sync routine started with PID: $ROUTINE_PID"

    # Wait only in DEV mode
    if [ "$exec_mode" = "DEV" ]; then
        wait
    fi
    ;;
"PYTHON_APP")
    APP_LOGFILE=$APP_LOGFILE.py_app.log
    echo_log 0 "  ..python application mode"
    echo_log 0 ""
    $PYTHON_CMD /app/src_python/qcvcolocapi/Colocation.py > $APP_LOGFILE
    ;;
"DEMO")
    APP_LOGFILE=$APP_LOGFILE.py_app.log
    echo_log 0 "python application mode : DEMO"
    echo_log 0 "ending demo."
    ;;
"BASH")
    echo_log 0 "  ..bash mode"
    echo_log 0 ""
    /bin/bash
    ;;
*)
    echo_log 0 "  no option passed ..... defaulting to bash "
    /bin/bash
    ;;
esac

echo_log 0 ""
echo_log 0 "  launcher exited..."
