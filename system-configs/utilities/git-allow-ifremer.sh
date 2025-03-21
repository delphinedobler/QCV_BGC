git-add-creds.sh ${GIT_IFR_GITREPO} ${GIT_IFR_LOGIN} ${GIT_IFR_PASS}
if [ "$(id -u)" == 0  ]; then
    echo "copying from serviceuser"
    cp /home/serviceuser/.git-credentials /root
fi