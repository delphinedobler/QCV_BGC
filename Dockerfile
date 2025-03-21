#                        ██████╗  █████╗ ███████╗███████╗                    
#                        ██╔══██╗██╔══██╗██╔════╝██╔════╝                    
#    █████╗    █████╗    ██████╔╝███████║███████╗█████╗      █████╗    █████╗
#    ╚════╝    ╚════╝    ██╔══██╗██╔══██║╚════██║██╔══╝      ╚════╝    ╚════╝
#                        ██████╔╝██║  ██║███████║███████╗                    
#                        ╚═════╝ ╚═╝  ╚═╝╚══════╝╚══════╝                    


FROM ubuntu:24.04 as python-base


# ---- ARGS & ENV
ARG FORCE_USER_ID
ARG FORCE_GROUP_ID
ARG FORCE_USER_PWD

# ENV FORCE_USER_ID=${FORCE_USER_ID}
# ENV FORCE_GROUP_ID=${FORCE_GROUP_ID}

# ----  Static ENV for build

ENV DEBIAN_FRONTEND=noninteractive

ENV PYTHONNOUSERSITE=True
ENV PYTHONUNBUFFERED=1
ENV PYTHONDONTWRITEBYTECODE=1

ENV MAMBA_ROOT_PREFIX=/opt/miniforge3
ENV MAMBA_EXE=/opt/miniforge3/bin/mamba



# ------------------ ROOT USER ACTIONS ------------------
USER root

# ---- INSTALL UTILS & PYTHON
RUN apt-get update \
    && apt-get --yes install apt-utils \
    && apt-get --yes install curl wget nano vim \
    && apt-get --yes install openssh-client \
    && apt-get --yes install sudo \
    && apt-get --yes install gosu ncdu\
    && apt-get --yes install net-tools \
    && apt-get --yes install dnsutils \
    && apt-get --yes install iproute2 \
    && apt-get --yes install nmap \
    && apt-get --yes install iputils-ping \
    && apt-get --yes install git \
    && apt-get --yes install ftp lftp \
    && apt-get --yes install psmisc \
    && apt-get --yes install unixodbc unixodbc-dev libaio-dev odbcinst

RUN apt-get update
RUN apt-get --yes upgrade

# Configure locales --> needed for R
# Set the locale
RUN apt-get clean && apt-get update && apt-get install -y locales locales-all
RUN locale-gen en_US.UTF-8
RUN dpkg-reconfigure locales
RUN echo LC_ALL=en_US.UTF-8 >> /etc/environment
RUN echo LANG=en_US.UTF-8 >> /etc/environment
RUN update-locale LC_ALL=en_US.UTF-8
RUN update-locale LANG=en_US.UTF-8

# ---- SET UID/GID
# add a user with same GID and UID as the host user that owns the workspace files on the host (bind mount)
RUN groupadd -f servicegroup -g ${FORCE_GROUP_ID}
RUN useradd -s $(which bash) --uid ${FORCE_USER_ID} --gid ${FORCE_GROUP_ID} -m serviceuser
RUN echo serviceuser:${FORCE_USER_PWD}
RUN echo "serviceuser:${FORCE_USER_PWD}" | chpasswd
RUN usermod -aG sudo serviceuser
RUN echo "serviceuser ALL=(ALL) NOPASSWD:SETENV: /usr/local/bin/entrypoint.sh" >> /etc/sudoers


# ---- embedded resources
RUN mkdir -p /embedded-resources

# ---- RUNTIME folders
RUN mkdir -p /runtime/config
RUN mkdir -p /runtime/log
RUN mkdir -p /runtime/run
RUN mkdir -p /runtime/data-in
RUN mkdir -p /runtime/data-out


# ---- service user specific folders
RUN mkdir -p /home/serviceuser/.ssh
RUN mkdir -p /home/serviceuser/service
RUN mkdir -p /home/serviceuser/informations
RUN mkdir -p /home/serviceuser/service/dev-external-libs


# ---- GIVE PERMISSIONS
RUN chown -R serviceuser:servicegroup /home/serviceuser
RUN chown -R serviceuser:servicegroup /runtime
RUN chmod -R 777 /home/serviceuser/service/dev-external-libs





# ---- COPY & RIGHTS ENTRYPOINT
COPY ./docker-entrypoint/entrypoint.sh /usr/local/bin/entrypoint.sh
RUN chmod +x /usr/local/bin/entrypoint.sh


# ---- SET DEFAULT IMAGE ENTRYPOINT ( better to comment that for Bluecloud compatibility )
#          Related to - ENTRYPOINT_BYPASS - env var
# ENTRYPOINT ["sudo", "-E", "/usr/local/bin/entrypoint.sh"]



#                 ███╗   ███╗ █████╗ ███╗   ███╗██████╗  █████╗     ███████╗████████╗ █████╗  ██████╗ ███████╗                
#                 ████╗ ████║██╔══██╗████╗ ████║██╔══██╗██╔══██╗    ██╔════╝╚══██╔══╝██╔══██╗██╔════╝ ██╔════╝                
# █████╗█████╗    ██╔████╔██║███████║██╔████╔██║██████╔╝███████║    ███████╗   ██║   ███████║██║  ███╗█████╗      █████╗█████╗
# ╚════╝╚════╝    ██║╚██╔╝██║██╔══██║██║╚██╔╝██║██╔══██╗██╔══██║    ╚════██║   ██║   ██╔══██║██║   ██║██╔══╝      ╚════╝╚════╝
#                 ██║ ╚═╝ ██║██║  ██║██║ ╚═╝ ██║██████╔╝██║  ██║    ███████║   ██║   ██║  ██║╚██████╔╝███████╗                
#                 ╚═╝     ╚═╝╚═╝  ╚═╝╚═╝     ╚═╝╚═════╝ ╚═╝  ╚═╝    ╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚══════╝                
                                                                                                                            
FROM python-base as mamba_stage


# ---- ARGS & ENV ----
ARG DEVELOPER_GIT_USER_NAME
ARG DEVELOPER_GIT_USER_EMAIL
# ARG GIT_IFR_GITREPO
# ARG GIT_IFR_LOGIN
# ARG GIT_IFR_PASS

ENV DEVELOPER_GIT_USER_NAME=${DEVELOPER_GIT_USER_NAME}
ENV DEVELOPER_GIT_USER_EMAIL=${DEVELOPER_GIT_USER_EMAIL}


# ------------------ ROOT USER ------------------


# ---- COPY APP source code
COPY app /home/serviceuser/service/app

# ---- APP folder symlink
RUN ln -s /home/serviceuser/service/app /app

# ---- COPY embeded resources
COPY embedded-resources /embedded-resources


# ---- SSH EVAL configuration for serviceuser
RUN mkdir -p /home/serviceuser/.ssh
RUN chmod 700 /home/serviceuser/.ssh
COPY system-configs/SSH/config /home/serviceuser/.ssh/config
COPY dev-secrets-ssh/git-* /home/serviceuser/.ssh
RUN ls -la  /home/serviceuser/.ssh/
RUN chmod -R 600 /home/serviceuser/.ssh/*
RUN ssh-keyscan gitlab.com >> /home/serviceuser/.ssh/known_hosts
RUN eval "$(ssh-agent -s)" && \
    ssh-add /home/serviceuser/.ssh/git-ssh-key

RUN mkdir -p /home/serviceuser/bin
COPY system-configs/utilities/* /home/serviceuser/bin

# ---- SSH EVAL configuration for root
RUN mkdir -p /root/.ssh
RUN chmod 700 /root/.ssh
COPY system-configs/SSH/config /root/.ssh/config
COPY dev-secrets-ssh/git-* /root/.ssh
RUN ls -la  /root/.ssh/
RUN chmod -R 600 /root/.ssh/*
RUN ssh-keyscan gitlab.com >> /root/.ssh/known_hosts
RUN eval "$(ssh-agent -s)" && \
    ssh-add /root/.ssh/git-ssh-key

RUN mkdir -p /root/bin
COPY system-configs/utilities/* /root/bin

# ---- ensure service user ownership
RUN chown -R serviceuser:servicegroup /home/serviceuser
RUN chown -R serviceuser:servicegroup /runtime
RUN chown -R serviceuser:servicegroup /embedded-resources


# - - - - INSTALL MINIFORGE3 ( Globally : better UID:GID on the fly change )
WORKDIR /root
RUN curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh" 
RUN chmod +x Miniforge3-$(uname)-$(uname -m).sh
RUN bash Miniforge3-$(uname)-$(uname -m).sh -b -p $MAMBA_ROOT_PREFIX
RUN chmod 755 $MAMBA_ROOT_PREFIX
RUN rm -f Miniforge3-$(uname)-$(uname -m).sh

# - - - - ROOT Context : GIT additional credentials if any (using plaintext store...)
#   uncomment this line and add the following 3 ARGs to your env template & docker compose
# RUN /home/serviceuser/bin/git-add-creds.sh ${GIT_IFR_GITREPO} ${GIT_IFR_LOGIN} ${GIT_IFR_PASS}

# - - - - INSTALL MAMBA VIRTUAL ENV
WORKDIR /home/serviceuser/service/app
RUN $MAMBA_EXE update mamba
RUN $MAMBA_EXE env create -n env-develop -f environment.yml

# - - - -  Insert BASH content
RUN $MAMBA_EXE shell init --shell bash
RUN echo "PATH=${PATH}:/home/serviceuser/bin:/app/launchers"  >> ~/.bashrc
RUN echo "cd /app"  >> ~/.bashrc
RUN echo "mamba activate env-develop" >> ~/.bashrc


# ------------------ SERVICEUSER USER ------------------
USER serviceuser

WORKDIR /home/serviceuser


# - - - -  Git settings
RUN git config --global user.name "${DEVELOPER_GIT_USER_NAME}"
RUN git config --global user.email "${DEVELOPER_GIT_USER_EMAIL}"
# RUN git config --global pull.rebase true
RUN git config --global pull.ff only
RUN git config --global merge.ff false

# - - - - SERVICEUSER Context : GIT additional credentials if any (using plaintext store...)
#   uncomment this line and add the following 3 ARGs to your env template & docker compose
# RUN /home/serviceuser/bin/git-add-creds.sh ${GIT_IFR_GITREPO} ${GIT_IFR_LOGIN} ${GIT_IFR_PASS}

#       Default shell for next mamba envs installs
# SHELL ["$MAMBA_EXE", "run", "-n", "env-develop", "/bin/bash", "-c"]

# - - - -  Insert BASH content
RUN $MAMBA_EXE shell init --shell bash
RUN echo "PATH=${PATH}:/home/serviceuser/bin:/app/launchers"  >> ~/.bashrc
RUN echo "cd /app"  >> ~/.bashrc
RUN echo "mamba activate env-develop" >> ~/.bashrc



# ------------------ ROOT USER ------------------

USER root

# lazy rule for volume mounts
RUN chmod -R 777 /runtime

# - - - - clean user dev-secrets to ensure it is not pushed with the container image
RUN rm -f /home/serviceuser/.ssh/*
RUN rm -f /home/serviceuser/.git-credentials
RUN rm -f /root/.ssh/*
RUN rm -f /root/.git-credentials


# - - - - default startup command that allows to change user at container start
CMD ["/usr/local/bin/entrypoint.sh"]
