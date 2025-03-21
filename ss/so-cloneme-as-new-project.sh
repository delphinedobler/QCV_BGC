#!/bin/bash

# TODO : script like loufiat !

start_project.sh
#######     CLONAGE BOILER PLATE     #######
# demande le nom du nouveau projet nouveauprojet
# demande le sous-repertoire cible : ou on veut copier la boilerplate 
# demande si il veut du $INSTALL_OPT = FRESH / PYSTARTER ; valeur par defaut = PYSTARTER
# demande la verison de python $PYTHON_VERSION; valeur par default = 3.10

# clone le git boiler plate dans le repertoire cible/nouveauprojet
# supprime le repertoire .git
# copie le environment.yml de resources/INSTALL_OPT dans app/
# sed le app/environment.yml pour remplacer la version de python 
# cd --> met d'utilisateur 
# lance le service-init.sh (init local) avec l'option --change-opt=$INSTALL_OPT --change-pyversion=$PYTHON_VERSION

service-init.sh
#######     INIT LOCAL     #######
# lire les options parameters (--change-opt ; --change-pyversion)
    
# si .env existe, source .env si non crée le .env à partir du .env.example 
# replace les default values ou current values

# demande les informations suivantes à l'utilisateur, affiche la current value, et demande si on la conserve
    # - NETWORK_FRONT_INTERNAL_PORT
    # - EXTERNAL_PORT 
    # - DEVELOPER_GIT_REPO_NAME ; default = ''
    # - DOCKER_NAME
    # - DOCKER_NAMESPACE
    # - nom de la clef ssh (sans le .pub); valeur par defaut = id_rsa
# fait les remplacements prompt depuis les demandes

# Si $INSTALL_OPT renseigné dans --change-opt ==> sed dans le .env la version de INSTALL_OPT
# Si $PYTHON_VERSION renseigné dans --change-pyversion ==> sed dans le .env la version de PYTHON

# fait les remplacements automatiques l'UID/GID

# copie les clefs ssh et securise le dossier secret

# build le docker
# change la variable d'environnement INSTALL_OPT=LOCK

