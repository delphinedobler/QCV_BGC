# Ce script initialise un projet local
# Ce qu'il fait : 
    # lire les options parameters (--change-opt ; --change-pyversion)
    # si .env existe, source .env et replace les default values ou current values
    # demande les informations suivantes à l'utilisateur, affiche la current value, et demande si on la conserve
        # - NETWORK_FRONT_INTERNAL_PORT
        # - NETWORK_FRONT_EXPOSED_PORT 
        # - DEVELOPER_GIT_REPO_NAME ; default = ''
        # - DOCKER_NAME
        # - DOCKER_NAMESPACE
        # - nom de la clef ssh (sans le .pub); valeur par defaut = id_rsa

    # crée le .env à partir du .env.example si inexistant
    # Si $INSTALL_OPT renseigné dans --change-opt ==> sed dans le .env la version de INSTALL_OPT
    # Si $PYTHON_VERSION renseigné dans --change-pyversion ==> sed dans le .env la version de PYTHON
    # fait les remplacements automatiques l'UID/GID
    # fait les remplacements prompt depuis les demandes
    # copie les clefs ssh et securise le dossier secret
    # build le docker (avec option --no-cache directement dans les options de lancement)
    # change la variable d'environnement INSTALL_OPT=LOCK

#!/bin/bash

#       where the script is
script_dir="$(cd "$(dirname "$0")" && pwd)"

#       Is the terminal is a zsh ? --> if the command exist
# Check if the command is available
if command -v $zsh_try_command &> /dev/null; then
    # The command exists, so you can execute it
    echo "Command '$zsh_try_command' exists, executing..."
    ZSH_VERSION=$(zsh --version| awk -F" " '{print $2}')
    echo $ZSH_VERSION
else
    echo "blblblblbl"
    ZSH_VERSION=
fi



# ██████╗ ███████╗███████╗ █████╗ ██╗   ██╗██╗  ████████╗    ██╗   ██╗ █████╗ ██╗     ██╗   ██╗███████╗███████╗
# ██╔══██╗██╔════╝██╔════╝██╔══██╗██║   ██║██║  ╚══██╔══╝    ██║   ██║██╔══██╗██║     ██║   ██║██╔════╝██╔════╝
# ██║  ██║█████╗  █████╗  ███████║██║   ██║██║     ██║       ██║   ██║███████║██║     ██║   ██║█████╗  ███████╗
# ██║  ██║██╔══╝  ██╔══╝  ██╔══██║██║   ██║██║     ██║       ╚██╗ ██╔╝██╔══██║██║     ██║   ██║██╔══╝  ╚════██║
# ██████╔╝███████╗██║     ██║  ██║╚██████╔╝███████╗██║        ╚████╔╝ ██║  ██║███████╗╚██████╔╝███████╗███████║
# ╚═════╝ ╚══════╝╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚══════╝╚═╝         ╚═══╝  ╚═╝  ╚═╝╚══════╝ ╚═════╝ ╚══════╝╚══════╝

change_opt_default=
change_pyversion_default=

START_OPTION_DEFAULT="BASH"
INTERNAL_PORT_DEFAULT=5585
EXPOSED_PORT_DEFAULT=5585
DEVELOPER_GIT_REPO_NAME_DEFAULT=git@gitlab.com:clement.weber/testproj_boilerplate.git
DOCKER_NAME_DEFAULT="testproj_docker"
DOCKER_NAMESPACE_DEFAULT="clementweber"
VOL_DATA_DEFAULT=/datafiles2/EXT/ARGO/content/202307-ArgoData
DEVELOPER_SSH_KEY_PRIVATE_DEFAULT=/nfs-home/pokapok/cweber/.ssh/pokapok_cweber




#    █████╗    █████╗    █████╗    █████╗    █████╗    █████╗    █████╗    █████╗    █████╗    █████╗    █████╗    █████╗    █████╗    █████╗    █████╗
#    ╚════╝    ╚════╝    ╚════╝    ╚════╝    ╚════╝    ╚════╝    ╚════╝    ╚════╝    ╚════╝    ╚════╝    ╚════╝    ╚════╝    ╚════╝    ╚════╝    ╚════╝



# ███████╗██╗   ██╗███╗   ██╗ ██████╗████████╗██╗ ██████╗ ███╗   ██╗███████╗
# ██╔════╝██║   ██║████╗  ██║██╔════╝╚══██╔══╝██║██╔═══██╗████╗  ██║██╔════╝
# █████╗  ██║   ██║██╔██╗ ██║██║        ██║   ██║██║   ██║██╔██╗ ██║███████╗
# ██╔══╝  ██║   ██║██║╚██╗██║██║        ██║   ██║██║   ██║██║╚██╗██║╚════██║
# ██║     ╚██████╔╝██║ ╚████║╚██████╗   ██║   ██║╚██████╔╝██║ ╚████║███████║
# ╚═╝      ╚═════╝ ╚═╝  ╚═══╝ ╚═════╝   ╚═╝   ╚═╝ ╚═════╝ ╚═╝  ╚═══╝╚══════╝

#       Function to prompt the user for input with options
function prompt_user {
    local variable_name="$1"
    local default_value="$2"
    local current_value="$3"

    # TODO: si on demande la DEVELOPER_SSH_KEY il faut absolument dire qu'il faut le path absolu
    if [[ -z "$current_value" ]]; then
        read -p "Enter new value for $variable_name (default/current/overwrite): (default=$default_value)" new_value
    else
        read -p "Enter new value for $variable_name (default/current/overwrite): (current=$current_value)" new_value
    fi

    case "$new_value" in
        "default")
            out_value=$default_value
            ;;
        "current")
            out_value=$current_value
            ;;
        "overwrite")
            read -p "Enter the new value for $variable_name: " new_variable_value
            # eval "$variable_name=\"$new_variable_value\""
            out_value=$new_variable_value
            ;;
        *)
            if [[ -z "$new_value" ]]; then
                if [[ -z "$current_value" ]]; then
                    out_value=$default_value
                else
                    out_value=$current_value
                fi
            else
                out_value=$new_value
            fi
            ;;
    esac
    # OUTPUT OF THE FUNCTION
    echo "$out_value"
}

#       Function to replace the value in the .env file
function replace_value_in_env {
    local variable_name="$1"
    local new_value="$2"

    echo "------"
    pwd
    sed_original="sed -i '' '/^'$variable_name'=/s|^'$variable_name'=.*|'$variable_name'='$new_value'|' .env"
    test_zsh_4_sed "$sed_original"
    echo "------"
}

#       function to prompt and replace in .env
function prompt_and_replace {
    local variable_name="$1"
    local default_value="$2"
    local current_value="$3"

    # Prompt the user for input and get the new value
    local new_value
    new_value=$(prompt_user "$variable_name" "$default_value" "$current_value")

    # Replace the value in the .env file
    replace_value_in_env "$variable_name" "$new_value"

    # OUTPUT OF FUNCTION
    echo $new_value
}

function test_zsh_4_sed() {
    local command="$1"

    echo "zsh testing ..."
    if [ -n "$ZSH_VERSION" ]; then
        eval $command

    else
        eval "${command//"''"/}"
    fi
}

#    █████╗    █████╗    █████╗    █████╗    █████╗    █████╗    █████╗    █████╗    █████╗    █████╗    █████╗    █████╗    █████╗    █████╗    █████╗
#    ╚════╝    ╚════╝    ╚════╝    ╚════╝    ╚════╝    ╚════╝    ╚════╝    ╚════╝    ╚════╝    ╚════╝    ╚════╝    ╚════╝    ╚════╝    ╚════╝    ╚════╝


# ██╗███╗   ██╗██╗████████╗
# ██║████╗  ██║██║╚══██╔══╝
# ██║██╔██╗ ██║██║   ██║   
# ██║██║╚██╗██║██║   ██║   
# ██║██║ ╚████║██║   ██║   
# ╚═╝╚═╝  ╚═══╝╚═╝   ╚═╝   

#       Check if .env file exists and source it (if it does) to load environment variables
if [[ -f .env ]]; then
    echo " ------> .env already exists => load info"
    source .env
else
    echo ------> new env from example
    cp .env.example .env
    # source .env
fi
echo "-------"
echo $INSTALL_OPT


#  ██████╗ ██████╗ ████████╗███████╗       ██╗        █████╗ ██████╗  ██████╗ ███████╗
# ██╔═══██╗██╔══██╗╚══██╔══╝██╔════╝       ██║       ██╔══██╗██╔══██╗██╔════╝ ██╔════╝
# ██║   ██║██████╔╝   ██║   ███████╗    ████████╗    ███████║██████╔╝██║  ███╗███████╗
# ██║   ██║██╔═══╝    ██║   ╚════██║    ██╔═██╔═╝    ██╔══██║██╔══██╗██║   ██║╚════██║
# ╚██████╔╝██║        ██║   ███████║    ██████║      ██║  ██║██║  ██║╚██████╔╝███████║
#  ╚═════╝ ╚═╝        ╚═╝   ╚══════╝    ╚═════╝      ╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝ ╚══════╝

#       bloc : passing options and arguments 
while [[ $# -gt 0 ]]; do
    case "$1" in
        --change-opt)
            shift
            change_opt="$1"
            ;;
        --change-pyversion)
            shift
            change_pyversion="$1"
            ;;
        --use-cache)
            shift
            use_cache="$1"
            ;;
        *)
            # Handle other arguments or options as needed
            ;;
    esac
    shift
done


#       If --change-opt option was not provided, use the value from .env or set a default
if [[ -z "$change_opt" ]]; then
    if [[ -n "$INSTALL_OPT" ]]; then
        change_opt="$INSTALL_OPT"
    else
        change_opt="$change_opt_default"
    fi
else
    replace_value_in_env "INSTALL_OPT" "$change_opt"
fi

#       If --change-pyversion option was not provided, use the value from .env or set a default
if [[ -z "$change_pyversion" ]]; then
    if [[ -n "$PYTHON_VERSION" ]]; then
        change_pyversion="$PYTHON_VERSION"
    else
        change_pyversion="$change_opt_default"
    fi
else
    replace_value_in_env "PYTHON_VERSION" "$change_pyversion"
fi




# ██████╗ ██████╗  ██████╗ ███╗   ███╗██████╗ ████████╗     █████╗ ███╗   ██╗██████╗     ██████╗ ███████╗██████╗ ██╗      █████╗  ██████╗███████╗
# ██╔══██╗██╔══██╗██╔═══██╗████╗ ████║██╔══██╗╚══██╔══╝    ██╔══██╗████╗  ██║██╔══██╗    ██╔══██╗██╔════╝██╔══██╗██║     ██╔══██╗██╔════╝██╔════╝
# ██████╔╝██████╔╝██║   ██║██╔████╔██║██████╔╝   ██║       ███████║██╔██╗ ██║██║  ██║    ██████╔╝█████╗  ██████╔╝██║     ███████║██║     █████╗  
# ██╔═══╝ ██╔══██╗██║   ██║██║╚██╔╝██║██╔═══╝    ██║       ██╔══██║██║╚██╗██║██║  ██║    ██╔══██╗██╔══╝  ██╔═══╝ ██║     ██╔══██║██║     ██╔══╝  
# ██║     ██║  ██║╚██████╔╝██║ ╚═╝ ██║██║        ██║       ██║  ██║██║ ╚████║██████╔╝    ██║  ██║███████╗██║     ███████╗██║  ██║╚██████╗███████╗
# ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚═╝╚═╝        ╚═╝       ╚═╝  ╚═╝╚═╝  ╚═══╝╚═════╝     ╚═╝  ╚═╝╚══════╝╚═╝     ╚══════╝╚═╝  ╚═╝ ╚═════╝╚══════╝


prompt_and_replace "APP_START_OPTION" "$START_OPTION_DEFAULT" "$APP_START_OPTION"
prompt_and_replace "NETWORK_BACK_INTERNAL_PORT" "$INTERNAL_PORT_DEFAULT" "$NETWORK_BACK_INTERNAL_PORT"
prompt_and_replace "NETWORK_BACK_EXPOSED_PORT" "$EXPOSED_PORT_DEFAULT" "$NETWORK_BACK_EXPOSED_PORT"
prompt_and_replace "NETWORK_FRONT_INTERNAL_PORT" "$INTERNAL_PORT_DEFAULT" "$NETWORK_FRONT_INTERNAL_PORT"
prompt_and_replace "NETWORK_FRONT_EXPOSED_PORT" "$EXPOSED_PORT_DEFAULT" "$NETWORK_FRONT_EXPOSED_PORT"
prompt_and_replace "DEVELOPER_GIT_REPO_NAME" "$DEVELOPER_GIT_REPO_NAME_DEFAULT" "$DEVELOPER_GIT_REPO_NAME"
prompt_and_replace "DOCKER_NAME" "$DOCKER_NAME_DEFAULT" "$DOCKER_NAME"
prompt_and_replace "DOCKER_NAMESPACE" "$DOCKER_NAMESPACE_DEFAULT" "$DOCKER_NAMESPACE"
prompt_and_replace "VOL_DATA" "$VOL_DATA_DEFAULT" "$VOL_DATA"
prompt_and_replace "DEVELOPER_SSH_KEY_PRIVATE" "$DEVELOPER_SSH_KEY_PRIVATE_DEFAULT" "$DEVELOPER_SSH_KEY_PRIVATE"

DEVELOPER_SSH_KEY_PRIVATE=$(cat .env | grep DEVELOPER_SSH_KEY_PRIVATE | awk -F"=" '{print $2}')

DEVELOPER_SSH_KEY_PUBLIC=$DEVELOPER_SSH_KEY_PRIVATE.pub

if [[ -z "$DEVELOPER_SSH_KEY_PUBLIC" ]]; then
    replace_value_in_env "DEVELOPER_SSH_KEY_PUBLIC" ""
else
    replace_value_in_env "DEVELOPER_SSH_KEY_PUBLIC" "$DEVELOPER_SSH_KEY_PUBLIC"
    
    # dev-secrets folder : copy keys, change permissions
    test -e $DEVELOPER_SSH_KEY_PUBLIC && cp $DEVELOPER_SSH_KEY_PUBLIC ./dev-secrets/git-ssh-key.pub
    test -e $DEVELOPER_SSH_KEY_PRIVATE && cp $DEVELOPER_SSH_KEY_PRIVATE ./dev-secrets/git-ssh-key
    # change permissions
    chmod 700 ./dev-secrets
    chmod 600 ./dev-secrets/*

fi


#       remplacement auto FORCE_USER_ID/FURCE_GROUP_ID
MY_UID=$(id -u)
MY_GID=$(id -g)

sed_original="sed -i '' '/^FORCE_USER_ID=/s/=.*/='$MY_UID'/' .env"
test_zsh_4_sed "$sed_original"
sed_original="sed -i '' '/^FORCE_GROUP_ID=/s/=.*/='$MY_GID'/' .env"
test_zsh_4_sed "$sed_original"


# ███████╗███████╗████████╗    ███╗   ██╗███████╗██╗    ██╗     ██████╗ ██╗████████╗
# ██╔════╝██╔════╝╚══██╔══╝    ████╗  ██║██╔════╝██║    ██║    ██╔════╝ ██║╚══██╔══╝
# ███████╗█████╗     ██║       ██╔██╗ ██║█████╗  ██║ █╗ ██║    ██║  ███╗██║   ██║   
# ╚════██║██╔══╝     ██║       ██║╚██╗██║██╔══╝  ██║███╗██║    ██║   ██║██║   ██║   
# ███████║███████╗   ██║       ██║ ╚████║███████╗╚███╔███╔╝    ╚██████╔╝██║   ██║   
# ╚══════╝╚══════╝   ╚═╝       ╚═╝  ╚═══╝╚══════╝ ╚══╝╚══╝      ╚═════╝ ╚═╝   ╚═╝   

#       get the DEVELOPER_GIT_REPO_NAME variable 
DEVELOPER_GIT_REPO_NAME=$(cat .env | grep DEVELOPER_GIT_REPO_NAME | awk -F"=" '{print $2}')

echo    # (optional) move to a new line
read -p "Are you sure for the new git? " -n 1 -r

if [[ $REPLY =~ ^[Yy]$ ]]
then
    git init --initial-branch=main
    echo ------> new git repo ready.
    echo $DEVELOPER_GIT_REPO_NAME

fi

echo

#       Replace "DEVELOPER_GIT_REPO_NAME" by "GIT_RN_MODIFIED" in .env
config_name_std="pok.gitlab.com"
GIT_RN_MODIFIED="${DEVELOPER_GIT_REPO_NAME/git@gitlab.com/$config_name_std}"
replace_value_in_env "DEVELOPER_GIT_REPO_NAME" "$GIT_RN_MODIFIED"

#       Link local git and gitlab repo      
git remote add origin $GIT_RN_MODIFIED
echo new git remote add
git remote -v




# ██████╗ ██╗   ██╗██╗██╗     ██████╗ 
# ██╔══██╗██║   ██║██║██║     ██╔══██╗
# ██████╔╝██║   ██║██║██║     ██║  ██║
# ██╔══██╗██║   ██║██║██║     ██║  ██║
# ██████╔╝╚██████╔╝██║███████╗██████╔╝
# ╚═════╝  ╚═════╝ ╚═╝╚══════╝╚═════╝ 

#       Get the INSTALL_OPT variable
INSTALL_OPT=$(cat .env | grep INSTALL_OPT | awk -F"=" '{print $2}')
echo $INSTALL_OPT

#       source again in case of
source .env

#       If --use-cache option was not provided, use the cache as default
if [[ -z "$no_cache" ]]; then
    docker-compose build --force-rm --progress plain --no-cache
    echo "USE NO CACHE"
else
    docker-compose build --force-rm --progress plain
    echo "USE CACHE"
fi


#  ██████╗██╗  ██╗ █████╗ ███╗   ██╗ ██████╗ ███████╗    ████████╗ ██████╗                ██╗  ██╗      ██████╗  ██████╗██╗  ██╗
# ██╔════╝██║  ██║██╔══██╗████╗  ██║██╔════╝ ██╔════╝    ╚══██╔══╝██╔═══██╗               ╚██╗ ██║     ██╔═══██╗██╔════╝██║ ██╔╝
# ██║     ███████║███████║██╔██╗ ██║██║  ███╗█████╗         ██║   ██║   ██║    █████╗█████╗╚██╗██║     ██║   ██║██║     █████╔╝ 
# ██║     ██╔══██║██╔══██║██║╚██╗██║██║   ██║██╔══╝         ██║   ██║   ██║    ╚════╝╚════╝██╔╝██║     ██║   ██║██║     ██╔═██╗ 
# ╚██████╗██║  ██║██║  ██║██║ ╚████║╚██████╔╝███████╗       ██║   ╚██████╔╝               ██╔╝ ███████╗╚██████╔╝╚██████╗██║  ██╗
#  ╚═════╝╚═╝  ╚═╝╚═╝  ╚═╝╚═╝  ╚═══╝ ╚═════╝ ╚══════╝       ╚═╝    ╚═════╝                ╚═╝  ╚══════╝ ╚═════╝  ╚═════╝╚═╝  ╚═╝
replace_value_in_env "INSTALL_OPT" "LOCK"
