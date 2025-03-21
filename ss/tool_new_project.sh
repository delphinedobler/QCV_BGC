#!/bin/bash


############################################################
# Help                                                     #
############################################################
Help()
{
    # Display Help
    echo "This program to assist the new project docker creation"
    echo
    echo "options:"
    echo "--help            Display this help"
    echo "--nameproject     Name the new project, no default, no name --> no new project"
    echo "--pathproject     path of the new project, default = ."
    echo "--namebranch      Name the boilerplate branch you wanna use . avail : main, simple_install"
    echo
}
############################################################
# END                                                    #
############################################################

# This part open .env in nanon neditor waiting for modifications
# Check if .env file exists
if [ -f .env ]; then
    nano .env
else
    echo ".env file not found."
fi


