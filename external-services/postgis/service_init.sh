#!/bin/bash

NOW=$(date --utc +%Y%m%d%H%M%S)

# - - functions
source ../../ss/tool_functions.sh

# - - start options
if [ $# -gt 0 ]; then
    if [ $1 == "-f" ]; then
        echo_indent 2 "Force removing previous ../.env"
        mv .env .env.old.$NOW
        echo_indent 2 ""
    fi
else
    echo ..
fi

if [ -f ".env" ]
then
    echo env exists
else
    echo new env from example
    cp .env.example .env
fi






