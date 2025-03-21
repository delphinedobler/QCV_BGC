#!/bin/bash


# Start backend services
echo "Starting backend services..."
if ! /home/serviceuser/service/app/launchers/start-back.sh; then
  echo "Error: Failed to start backend services."
  exit 1
fi

# Start frontend services
echo "Starting frontend services..."
if ! /home/serviceuser/service/app/launchers/start-front.sh; then
  echo "Error: Failed to start frontend services."
  exit 1
fi


echo "|=== Fullstack Application Started Successfully. ===|"