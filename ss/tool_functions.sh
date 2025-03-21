#!/bin/bash

# reusable tools for deployment scripts

quote_spaces() {
  local value="$1"

  # Check if the value contains spaces
  if [[ "$value" =~ [[:space:]] ]]; then
    echo "\"$value\""
  else
    echo "$value"
  fi
}


get_envvar() {
  local var_name="$1"
  local file="../.env"

  if [[ -f "$file" ]]; then
    grep -E "^${var_name}=" "$file" | cut -d '=' -f2
#   else
#     echo "File not found: $file"
  fi
}

get_dc_file(){
  if [ -z "$APP_DC_ALTERNATIVE" ]; then
    echo docker-compose.yml
  else
    echo docker-compose-$APP_DC_ALTERNATIVE.yml
  fi
}

# pretty indent for echoing in this script
echo_indent() {
    local indent_level=$1
    local content=$2
    local indent=""

    # Create indentation string based on indent level
    for ((i=0; i<indent_level; i++)); do
        indent+="    "
    done

    # Print the content with the indent
    echo "= ${indent}${content}"

}

ARCH_TYPE=$(uname)
case $ARCH_TYPE in
    Darwin)
        #echo_indent 2 "INFO: Arch = MacOS"
        DOCKER_PLATFORM=linux/arm64/v8
        ;;
    Linux)
        #echo_indent 2  "INFO: Arch = Linux"
        DOCKER_PLATFORM=linux/amd64
        ;;
    *)
        #echo_indent 2  "ERROR: Arch is not defined! Are you using Linux or MacOS? Exiting!"
        exit 1
        ;;         
esac

# sed working on Linux and MacOS
universal_sed() {
    local env_key="$1"
    local env_value="$2"
    local env_file=$3

    #echo "debug env_file = $env_file"


    case $ARCH_TYPE in
        Darwin)
            sed -i '' -e "s|$env_key=.*|$env_key=$env_value|g" $env_file
            ;;
        Linux)
            sed -i "s|$env_key=.*|$env_key=$env_value|g" $env_file
            ;;      
    esac

}

# Replace placeholders in default.conf with environment variables
universal_sed_placeholder() {
    local placeholder="$1"
    local value="$2"
    local file="$3"

    case "$(uname)" in
        Darwin)
            sed -i '' -e "s|$placeholder|$value|g" "$file"
            ;;
        Linux)
            sed -i "s|$placeholder|$value|g" "$file"
            ;;
    esac
}

