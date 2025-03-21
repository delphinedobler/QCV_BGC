#!/bin/bash


# execution log
APP_LOGDIR=/home/serviceuser/log
APP_PGCLIENT_LOG=$APP_LOGDIR/pok_front_app.log

touch $APP_PGCLIENT_LOG


python /home/serviceuser/service/app/src_python/tests/test_pg_client.py
