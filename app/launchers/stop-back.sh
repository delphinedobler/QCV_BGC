#!/bin/bash
lastrun=$(ls -t /runtime/log/nohup_back_*.log | head -n1)
PID=$(cat $lastrun)

echo Trying to kill : $PID
kill -9 $PID
echo
echo ----
ps -ux | grep gunicorn
echo ----
echo
echo Trying to kill all gunicorn processes

killall gunicorn

killall python