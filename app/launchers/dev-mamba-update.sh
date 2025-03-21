#!/bin/bash

echo "must sudo.. enter your passwd"

echo "Getting serviceuser ssh"
chmod 777 -R /home/serviceuser/.ssh
sudo cp -r /home/serviceuser/.ssh /root/

chmod 700 /home/serviceuser/.ssh
chmod 600 /home/serviceuser/.ssh/*

sudo chmod 700 /root/.ssh
sudo ls -la /root/.ssh
sudo chmod 600 -R /root/.ssh/

if [ "$(id -u)" != 0  ]; then
    sudo /opt/miniforge3/bin/mamba env update -n env-develop -f /app/environment.yml
fi

echo "Removing serviceuser ssh from root"
sudo rm -rf /root/.ssh

echo "Please launch dev-py-dl-pckgs.sh if you want to install some packages in development mode"