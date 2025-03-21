#!make
include .env
# export $(shell sed 's/=.*//' .env)

default: help


.PHONY: help
help: # Show help for each of the Makefile recipes.
	@grep -E '^[a-zA-Z0-9 -_]+:.*#'  Makefile | sort | while read -r l; do printf "\033[1;32m$$(echo $$l | cut -f 1 -d':')\033[00m:$$(echo $$l | cut -f 2- -d'#')\n"; done


# -------------------------------------------------------------------------------------------------
# - SERVICE TIDY UP
# -------------------------------------------------------------------------------------------------

# -- CLEANING
.PHONY: clean
clean: # reset git repo to a blank new one
	@./ss/service-clean-project.sh


# -------------------------------------------------------------------------------------------------
# - SERVICE GIT COMMANDS
# -------------------------------------------------------------------------------------------------


# -- GIT INIT
.PHONY: config_git_init
config_git_init: # reset git repo to a blank new one
	@./ss/service-git-init.sh


# -------------------------------------------------------------------------------------------------
# - SERVICE DOCKER COMMANDS
# -------------------------------------------------------------------------------------------------

# -- SERVICE INIT

.PHONY: config_init
config_init: # initialize the configuration with current user UID GID
	@./ss/service-init.sh

.PHONY: config_init_force
config_init_force: # initialize the configuration with current user UID GID
	@./ss/service-init.sh -f


# -- SERVICE RAZ

.PHONY: docker_raz
docker_raz:# kill all containers and remove image
	@docker kill ${DOCKER_NAME}
	@docker container prune
	@docker rmi $(DOCKER_NAME) || true


# -- SERVICE BUILD

.PHONY: docker_rebuild
docker_rebuild: draz dcb # clean and rebuild

.PHONY: dc_build
dc_build:# build docker images
	@./ss/service-build.sh

.PHONY: dcb
dcb:# build docker images
	@./ss/service-build.sh

.PHONY: dcb-no-cache
dcb-no-cache:# build docker images
	@./ss/service-build.sh --no-cache
# .PHONY: dc_down_and_save
# dc_down_and_save:# build docker images
# 	@docker-compose exec -u serviceuser pok_app /bin/bash -c /home/serviceuser/service/service-condalock.sh
# 	# @docker-compose down


# -- SERVICE UP / UPD

.PHONY: dc_up
dc_up:# run containers in background
	@./ss/service-up.sh	
.PHONY: dcu
dcu:# run containers in background
	@./ss/service-up.sh

.PHONY: dc_upd
dc_upd:# run containers in background
	@./ss/service-upd.sh
.PHONY: dcud
dcud:# run containers in background
	@./ss/service-upd.sh

.PHONY: services-up-fullstack
services-up-fullstack:# run containers in background
	@./ss/start-fullstack-services.sh


# -- SERVICE DOWN

.PHONY: dc_down
dc_down:# stop containers
	@./ss/service-down.sh
.PHONY: dcd
dcd:# stop containers
	@./ss/service-down.sh

.PHONY: services-down-fullstack
services-down-fullstack:# run containers in background
	@./ss/stop-fullstack-services.sh

# -- SERVICE PS

.PHONY: dc_ps
dc_ps:# show running containers
	@./ss/service-ps.sh
.PHONY: dcp
dcp:# show running containers
	@./ss/service-ps.sh

# -- SERVICE LOGS

.PHONY: dc_logs
dc_logs:# show containers logs
	@./ss/service-logs.sh
.PHONY: dcl
dcl:# show containers logs
	@./ss/service-logs.sh


# -- SERVICE EXEC

.PHONY: dc_exec_bash
dc_exec_bash:# enter bash into pok_app container
	@./ss/service-bash.sh
.PHONY: dceb
dceb:# enter bash into pok_app container
	@./ss/service-bash.sh






# -------------------------------------------------------------------------------------------------
# - DOCKER HUB
# -------------------------------------------------------------------------------------------------

# -- DOCKER PUSH

.PHONY: docker_push
docker_push: # docker push to repository
	@echo image name : ${DOCKER_REPONAME_V} 
	@docker push ${DOCKER_REPONAME_V}


.PHONY: docker_push_as_latest
docker_push_as_latest: # docker push to repository
	@echo image name : ${DOCKER_REPONAME_L}
	@docker tag ${DOCKER_REPONAME_V} ${DOCKER_REPONAME_L}
	@docker push ${DOCKER_REPONAME_L}


# -- DOCKER TAG

.PHONY: docker_tag
docker_tag: # docker tag
	@echo image name : ${DOCKER_REPONAME_V} 
	@docker tag ${DOCKER_REPONAME_V} ${DOCKER_REPONAME_V}


.PHONY: docker_tag_as_latest
docker_tag_as_latest: # docker tag
	@echo image name : ${DOCKER_REPONAME_L}
	@docker tag ${DOCKER_REPONAME_V} ${DOCKER_REPONAME_L}


# DOCKER CREDS / LOGIN

.PHONY: docker_rmcreds
docker_rmcreds: # docker registry, erase credentials
	@docker_logout
	@rm ~/.docker/config.json


.PHONY: docker_login
docker_login: # docker registry login
	@docker login


.PHONY: docker_logout
docker_logout:#  docker registry logout
	@docker logout


# -------------------------------------------------------------------------------------------------
# - DOCKER SYSTEM COMMANDS
# -------------------------------------------------------------------------------------------------

# -- RESTART DOCKER
.PHONY: docker_engine_restart
docker_engine_restart:#  docker engine restart
	@sudo service docker restart