#!/bin/bash

# stop backend processes
echo "Stopping backend processes..."
if ! /home/serviceuser/service/app/launchers/stop-back.sh; then
  echo "Error: Failed to stop backend processes."
  exit 1
fi

# stop frontend processes
echo "Stopping frontend processes..."
if ! /home/serviceuser/service/app/launchers/stop-front.sh; then
  echo "Error: Failed to stop frontend processes."
  exit 1
fi


echo "|=== Fullstack Application stopped Successfully. ===|"