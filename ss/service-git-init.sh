#!/bin/bash

DEVELOPER_GIT_REPO_NAME=$(cat ../.env | grep DEVELOPER_GIT_REPO_NAME | awk -F"=" '{print $2}')

read -p "Are you sure? " -n 1 -r
echo    # (optional) move to a new line
if [[ $REPLY =~ ^[Yy]$ ]]
then
    rm -rf ../.git
    echo previous git tree erased.
    git init --initial-branch=main
    echo new git repo ready.
fi

git remote add origin ${DEVELOPER_GIT_REPO_NAME}
echo new git remote add
