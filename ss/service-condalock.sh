#!/bin/bash

# Script to export conda and pip dependencies in the explicit way with the mamba export 
# and a resume of really direct deps with environment_used.yml file. 
source /home/serviceuser/.bashrc
cd /home/serviceuser/service
mamba_exec="/home/serviceuser/miniforge3/condabin/mamba"
echo $mamba_exec
$mamba_exec --version
pwd

# Function to extract packages used by python scriupts in app
extract_python_imports() {
    # Use find to locate all Python files in the app/ folder
    find app/ -type f -name "*.py" -print0 | while IFS= read -r -d $'\0' file; do
        # Use grep and sed to extract the first import or from statement and write to used.tmp
        grep -E '(import|from) [a-zA-Z_][a-zA-Z_0-9]*' "$file" | sed -nE 's/(import|from) ([a-zA-Z_][a-zA-Z_0-9]*).*/\2/p' >> used.tmp
    done
    sort used.tmp | uniq -u > used.yml
    rm used.tmp

}

# search packages actually used > used.yml
extract_python_imports

#       Path to the export file
export_yml="app/environment.yml"
# export env again if needed
$mamba_exec env export -n env-develop > $export_yml

#       Path to the used.yml file
used_yml="used.yml"

#       Name of the new file to create
output_file="app/environment_used.yml"
# delete previous one 
rm $output_file

# Read the list of words from used.yml into a variable
words=$(cat "$used_yml")

# Use awk to filter the dependencies in isit.yml
awk -v words="$words" -v output_file="$output_file" -F= '

    # Store the list of allowed words in an array
    BEGIN {
        split(words, allowed, "\n")
        output_temp = output_file ".temp"
        output_final = output_file
            }

    # Function to check if a line should be kept
    function should_keep(line) {
        if (line ~ /[- ]pip:/) {
            return 1
        }
        if (line ~ /^[- ]pip:/) {
            return 1
        }
        if (line ~ / - python=/ || line ~ / - poetry=/ || line ~ / - conda-lock=/) {
            return 1
        }
        for (i in allowed) {
            if (line ~ "- " allowed[i] "=") {
                return 1
            }
        }
            return 0
    }

    # Process the lines
    {
        if (in_channels) {
            if (!added_nodefaults) {
                print "  - nodefaults" >> output_final
                added_nodefaults = 1
            }
        }
        if (in_dependencies && should_keep($0)) {
            print >> output_final
        } 
        else if (in_dependencies) {
            print >> output_temp
        } 
        else {
            if (!in_dependencies) {
                print >> output_final
            }
        }
        if ($0 ~ /channels:/) {
            in_channels = 1
        }
        if ($0 ~ /dependencies:/) {
            in_dependencies = 1
        }

    }
' "$export_yml"


#       delete unnecessarily files
rm $used_yml
rm "$output_file.temp"