#!/bin/bash

# adds a new https credential to the git credentials cache

# param 1 : git repository FQDN
# param 2 : LOGIN
# param # : PASSWD

# activate cache
git config --global credential.helper store

# approve cache entry
echo "protocol=https
host="$1"
username="$2"
password="$3"
" | git credential approve
