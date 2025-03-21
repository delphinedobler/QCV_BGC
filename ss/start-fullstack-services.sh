#!/bin/bash

NETWORK_HOSTNAME=$(cat .env | grep NETWORK_HOSTNAME | awk -F"=" '{print $2}')


# Ensure a default config deployment file exists
SRC_CONFIG_FILE=./env-templates/config-deployment.yaml
TARGET_CONFIG_FILE=./dev-mount-config/config-deployment.yaml

if [ -f $TARGET_CONFIG_FILE ]; then
   echo "File $TARGET_CONFIG_FILE exists."
   echo ".. if not wanted ,  remove it yourself."
else
   echo "Copy config deployment-file in dev-mount-config..."
   cp $SRC_CONFIG_FILE $TARGET_CONFIG_FILE
fi


# Start app main service
if ! docker-compose up -d; then
    echo "Error: Failed to start main service."
    exit 1
  fi

# Start Docker services
services=("postgis" "cache-storage" "mongodb")
for service in "${services[@]}"; do
  echo "Starting $service..."
  if ! (cd ./external-services/$service && docker-compose up -d); then
    echo "Error: Failed to start $service."
    exit 1
  fi
done


# Start the proxy service last
echo "Creating conf file for proxy..."
cp ./external-services/proxy/proxy/conf/default.conf.template ./external-services/proxy/proxy/conf/default.conf
echo "Replacing placeholders in proxy configuration file..."
if ! ./ss/sed-proxy-conf.sh; then
  echo "Error: Failed to replace placeholders in proxy configuration file."
  exit 1
fi
echo "Starting proxy service..."
echo "Ensure ./external-services/proxy/proxy/conf/default.conf is configured correctly."
echo "Example: check proxy pass, ports and ensure server name can be resolved."
MYIP=`ip route get 1.1.1.1 | grep -oP 'src \K[0-9.]+'`
echo "NOTE : If you are deploying locally, you may need to map a server name to the external web ip. In this case, add the following line in /etc/hosts:¨
echo ¨ $MYIP $NETWORK_HOSTNAME"
if ! (cd ./external-services/proxy && docker-compose up -d); then
  echo "Error: Failed to start proxy service."
  exit 1
fi

echo "All services started successfully."