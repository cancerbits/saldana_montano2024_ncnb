#!/bin/bash

source parse_yaml.sh
eval $(parse_yaml ../config.yaml CONF_)

docker build -t $CONF_project_docker -f $CONF_project_root_host/$CONF_project_name.Dockerfile .
