#!/bin/bash
START_DATE=$(date --utc +%Y%m%d%H%M%S)
LAUNCHER_LOG=/runtime/log/nohup_back_$START_DATE.log
LAUNCHER_RUN=/runtime/run/nohup_back_$START_DATE.run


# Execute the command with nohup and redirect output to log file
nohup start-app.sh PYTHON_FASTAPI > $LAUNCHER_LOG 2>&1 &

# Get the PID of the last background process
PID=$!

# Print the PID to the log file
echo $PID > $LAUNCHER_RUN