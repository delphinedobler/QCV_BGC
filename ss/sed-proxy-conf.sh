
#!/bin/bash

SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
cd $SCRIPT_DIR

# - - functions
source ./tool_functions.sh

# Load environment variables from .env file if the line does not start with '#'
export $(grep -v '^#' ../.env | xargs)




universal_sed_placeholder "NETWORK_HOSTNAME" "$NETWORK_HOSTNAME" ../external-services/proxy/proxy/conf/default.conf
universal_sed_placeholder "NETWORK_FRONT_INTERNAL_PORT" "$NETWORK_FRONT_INTERNAL_PORT" ../external-services/proxy/proxy/conf/default.conf
universal_sed_placeholder "NETWORK_BACK_INTERNAL_PORT" "$NETWORK_BACK_INTERNAL_PORT" ../external-services/proxy/proxy/conf/default.conf