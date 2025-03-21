
#!/bin/bash

# pok.gitlab.com n'est pas dynamique, il faudrait une copie https et mettre la bp en publique

############################################################
# Help                                                     #
############################################################
Help()
{
    # Display Help
    echo "This program to assist the new project creation"
    echo
    echo "options:"
    echo "--help            Display this help"
    echo "--install-opt     Install option parameter : FRESH / PYSTARTER / LOCK "
    echo "--nameproject     Name the new project "
    echo "--location        Where the new project is created (absolute path)"
    echo "--pyversion       Choose the python version you want for this project"
    echo
}
############################################################
# END                                                    #
############################################################

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

INSTALL_OPT_DEFAULT="POETRY"
NEW_PROJECT_NAME_DEFAULT="new_project"
NEW_PROJECT_LOCATION_DEFAULT="."
PYTHON_VERSION_DEFAULT=3.10
DEVELOPER_GIT_REPO_NAME_DEFAULT=

#       source .env 
source .env


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
    sed -i '/^'$variable_name'=/s|^'$variable_name'=.*|'$variable_name'='$new_value'|' .env
}

#       Function to prompt and replace in .env
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

#       Function to prompt the install option
function prompt_install_opt {
    local variable_name="$1"
    local default_value="$2"
    local current_value="$3"

    # TODO: si on demande la DEVELOPER_SSH_KEY il faut absolument dire qu'il faut le path absolu
    read -p "Enter new value for $variable_name : (default=$default_value)" new_value

    case "$new_value" in
        "FRESH")
            out_value=FRESH
            ;;
        "FASTAPI")
            out_value=FASTAPI
            ;;
        "POETRY")
            out_value=POETRY
            ;;
        *)
            echo "wrong typing, please try again --> exiting project creation" 
            exit
            ;;
    esac
    # OUTPUT OF THE FUNCTION
    echo "$out_value"
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
        
#  ██████╗ ██████╗ ████████╗███████╗       ██╗        █████╗ ██████╗  ██████╗ ███████╗
# ██╔═══██╗██╔══██╗╚══██╔══╝██╔════╝       ██║       ██╔══██╗██╔══██╗██╔════╝ ██╔════╝
# ██║   ██║██████╔╝   ██║   ███████╗    ████████╗    ███████║██████╔╝██║  ███╗███████╗
# ██║   ██║██╔═══╝    ██║   ╚════██║    ██╔═██╔═╝    ██╔══██║██╔══██╗██║   ██║╚════██║
# ╚██████╔╝██║        ██║   ███████║    ██████║      ██║  ██║██║  ██║╚██████╔╝███████║
#  ╚═════╝ ╚═╝        ╚═╝   ╚══════╝    ╚═════╝      ╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝ ╚══════╝

#       bloc : passing options and arguments 
while [[ $# -gt 0 ]]; do
    case "$1" in
        --help)
            shift
            Help
            exit
            ;;
        --install-opt)
            shift
            install_opt="$1"
            ;;
        --nameproject)
            shift
            name_proj="$1"
            ;;
        --location)
            shift
            location="$1"
            ;;
        --pyversion)
            shift
            py_version="$1"
            ;;
        *)
            # Handle other arguments or options as needed
            ;;
    esac
    shift
done

#       Bloc to ask or get options from command line
if [[ -z "$install_opt" ]]; then
    INSTALL_OPT=$(prompt_install_opt "INSTALL_OPT" "$INSTALL_OPT_DEFAULT" "$INSTALL_OPT")
else
    INSTALL_OPT=$install_opt
fi
# ----------------
if [[ -z "$name_proj" ]]; then
    NEW_PROJECT_NAME=$(prompt_user "NEW_PROJECT_NAME" "$NEW_PROJECT_NAME_DEFAULT" "$NEW_PROJECT_NAME")
else
    NEW_PROJECT_NAME=$name_proj
fi
# ----------------
if [[ -z "$location" ]]; then
    NEW_PROJECT_LOCATION=$(prompt_user "NEW_PROJECT_LOCATION" "$NEW_PROJECT_LOCATION_DEFAULT" "$NEW_PROJECT_LOCATION")
else
    NEW_PROJECT_LOCATION=$location
fi
# ----------------
if [[ -z "$py_version" ]]; then
    PYTHON_VERSION=$(prompt_user "PYTHON_VERSION" "$PYTHON_VERSION_DEFAULT" "$PYTHON_VERSION")
else
    PYTHON_VERSION=$py_version
fi



#  ██████╗██╗      ██████╗ ███╗   ██╗███████╗     ██████╗ ██╗████████╗
# ██╔════╝██║     ██╔═══██╗████╗  ██║██╔════╝    ██╔════╝ ██║╚══██╔══╝
# ██║     ██║     ██║   ██║██╔██╗ ██║█████╗      ██║  ███╗██║   ██║   
# ██║     ██║     ██║   ██║██║╚██╗██║██╔══╝      ██║   ██║██║   ██║   
# ╚██████╗███████╗╚██████╔╝██║ ╚████║███████╗    ╚██████╔╝██║   ██║   
#  ╚═════╝╚══════╝ ╚═════╝ ╚═╝  ╚═══╝╚══════╝     ╚═════╝ ╚═╝   ╚═╝   

echo "- - ^"

#       Clone, checkout the mamba branch and remove the .git from boilerplate
git clone pok.gitlab.com:pokapok-projects/pkp0-reusable/python-boilerplates/py-boilerplate-base.git $NEW_PROJECT_LOCATION/$NEW_PROJECT_NAME --branch option-fastapi
git checkout option-fastapi
rm -rf $NEW_PROJECT_LOCATION/$NEW_PROJECT_NAME/.git

echo "- - ^"

# ███████╗███████╗████████╗    ███████╗███╗   ██╗██╗   ██╗██╗   ██╗███╗   ███╗██╗     
# ██╔════╝██╔════╝╚══██╔══╝    ██╔════╝████╗  ██║██║   ██║╚██╗ ██╔╝████╗ ████║██║     
# ███████╗█████╗     ██║       █████╗  ██╔██╗ ██║██║   ██║ ╚████╔╝ ██╔████╔██║██║     
# ╚════██║██╔══╝     ██║       ██╔══╝  ██║╚██╗██║╚██╗ ██╔╝  ╚██╔╝  ██║╚██╔╝██║██║     
# ███████║███████╗   ██║       ███████╗██║ ╚████║ ╚████╔╝██╗ ██║   ██║ ╚═╝ ██║███████╗
# ╚══════╝╚══════╝   ╚═╝       ╚══════╝╚═╝  ╚═══╝  ╚═══╝ ╚═╝ ╚═╝   ╚═╝     ╚═╝╚══════╝

#       copy the right environment.yml dependig of the INSTALL_OPT parameter
if [[ $INSTALL_OPT != "LOCK" ]]; then
    cp -f $script_dir/resources/$INSTALL_OPT/environment.yml $NEW_PROJECT_LOCATION/$NEW_PROJECT_NAME/app/.
    cp -f $script_dir/resources/$INSTALL_OPT/pyproject.toml $NEW_PROJECT_LOCATION/$NEW_PROJECT_NAME/app/.
else
    echo $INSTALL_OPT
fi

echo "- - ^"

#       replace python version in the environment.yml file
sed_original="sed -i '' "s/PYVERSION/${PYTHON_VERSION}/g" $NEW_PROJECT_LOCATION/$NEW_PROJECT_NAME/app/environment.yml"
test_zsh_4_sed "$sed_original"

#       replace poetry version in the environment.yml file
sed_original="sed -i '' "s/POETRYVERSION/${POETRY_VERSION}/g" $NEW_PROJECT_LOCATION/$NEW_PROJECT_NAME/app/environment.yml"
test_zsh_4_sed "$sed_original"

#       replace python version in the pyproject.toml file
sed_original="sed -i '' "s/PYVERSION/${PYTHON_VERSION}/g" $NEW_PROJECT_LOCATION/$NEW_PROJECT_NAME/app/pyproject.toml"
test_zsh_4_sed "$sed_original"


#  ██████╗  ██████╗ ██╗
# ██╔════╝ ██╔═══██╗██║
# ██║  ███╗██║   ██║██║
# ██║   ██║██║   ██║╚═╝
# ╚██████╔╝╚██████╔╝██╗
#  ╚═════╝  ╚═════╝ ╚═╝

#       go to the new project location
cd $NEW_PROJECT_LOCATION/$NEW_PROJECT_NAME
pwd
echo $INSTALL_OPT

./so-service-init.sh --change-opt $INSTALL_OPT --change-pyversion $PYTHON_VERSION
