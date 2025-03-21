#!/bin/bash

SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
cd $SCRIPT_DIR


# this script generates the .env.example template from fragments

# get all env example fragments
example_tpl=$(ls .././env-templates/.env.example.* | sort)

if [ -z "$example_tpl" ]; then
    echo " fatal ... No example skeleton files found in the current directory."
    exit 1
fi

# Concatenate the contents of the fragments files
cat $example_tpl > ../.env.example