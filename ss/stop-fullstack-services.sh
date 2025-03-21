#!/bin/bash

# Start app main service
if ! docker-compose down; then
    echo "Error: Failed to stop main service."
    exit 1
  fi

# Start Docker services
services=("postgis" "cache-storage" "mongodb")
for service in "${services[@]}"; do
  echo "Starting $service..."
  if ! (cd ./external-services/$service && docker-compose down); then
    echo "Error: Failed to stop $service."
    exit 1
  fi
done


# Start the proxy service last
echo "Stopping proxy..."
if ! (cd ./external-services/proxy && docker-compose down); then
  echo "Error: Failed to stop proxy service."
  exit 1
fi

echo "All services stopped successfully."