#!/bin/bash

# lastrun=$(ls -t /runtime/run/nohup_front_*.log | head -n1)
# PID=$(cat $lastrun)

# echo Trying to kill : $PID

# kill -9 $PID

killall shiny
killall python