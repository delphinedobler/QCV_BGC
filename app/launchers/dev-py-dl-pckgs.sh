#!/bin/bash

# Download python packages from the environment-dev.txt file in development mode


# be sure to have bashrc context
# source /home/serviceuser/.bashrc

# Disable Python bytecode compilation
export PYTHONDONTWRITEBYTECODE=1

APP_EXEC_MODE=${APP_EXEC_MODE:-DEV}

echo "must sudo.. enter your passwd"

echo "Application execution mode is $APP_EXEC_MODE"

if [ "$APP_EXEC_MODE" = "DEV" ]; then
    echo "Mode de développement détecté : installation des dépendances de développement"
    ENV_DEV_PATH=/app/environment-dev.txt
    echo "environment-dev.txt file is located at $ENV_DEV_PATH. Make sure that all the packages can be install (volumes mounted if necessary for example)"
    if [ -f "$ENV_DEV_PATH" ]; then

        # Array to store the paths
        declare -a extra_paths

        while IFS= read -r line; do
            # Check if the line starts with "-e" (editable mode)
            if [[ "$line" =~ ^-e ]]; then
                # Remove leading "-e" and extract the path
                path=$(echo "$line" | sed 's/^-e //')
                extra_paths+=("$path")
            fi
        done < "$ENV_DEV_PATH"

        # Display the paths collected
        echo "Collected paths:"
        for path in "${extra_paths[@]}"; do
            echo "$path"
        done
        # Change permissions to 777 for each path
        echo "Changing permissions to 777 for the collected paths..."
        for path in "${extra_paths[@]}"; do
            if [ -d "$path" ]; then
                sudo chmod 777 "$path"
                echo "Permissions changed for $path"
            else
                echo "Path $path does not exist or is not a directory"
            fi
        done

        sudo /opt/miniforge3/envs/env-develop/bin/python -m pip install --root-user-action=ignore --force-reinstall --no-cache-dir -r $ENV_DEV_PATH
    
        # Check if pip install was successful
        if [ $? -eq 0 ]; then
            echo "Installation des dépendances réussie !"

            # Disable Python bytecode compilation
            export PYTHONDONTWRITEBYTECODE=1

            # Ensure the file ends with a newline
            if [ "$(tail -c 1 "$ENV_DEV_PATH" | wc -l)" -eq 0 ]; then
                echo "" >> "$ENV_DEV_PATH"
            fi

            # Read the environment-dev.txt file and get the paths of the editable packages
            echo "---------------------------------------------------------------"
            echo "Fetching the paths from environment-dev.txt file..."
            echo "---------------------------------------------------------------"

            
            # Output the manual configuration instructions in case something goes wrong
            echo "---------------------------------------------------------------"
            echo "Manual Configuration Instructions for VSCode (if needed):"
            echo "---------------------------------------------------------------"
            echo "1. Open the Command Palette in VSCode (Ctrl+Shift+P)"
            echo "2. Search for 'Preferences: Open Workspace Settings (JSON)'"
            echo "3. Add the following configuration under 'settings':"
            echo "   \"settings\": {"
            echo "       \"python.analysis.extraPaths\": ["
            for path in "${extra_paths[@]}"; do
                echo "           \"$path\","
            done
            echo "       ]"
            echo "   }"
            echo "---------------------------------------------------------------"
            echo "Make sure to save the file. This will allow VSCode to recognize the local editable packages."
        else
            echo "L'installation des dépendances a échoué. Aucun changement n'a été effectué sur les paramètres de VSCode."
        fi
    else
        echo "No dev environment txt file found for $ENV_DEV_PATH"
    fi
else
    echo "Mode production : pas besoin d'installer des dépendances standards"
fi