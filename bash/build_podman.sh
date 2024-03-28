#!/bin/bash

source parse_yaml.sh
eval $(parse_yaml ../config.yaml CONF_)

podman build --shm-size=4096m -t $CONF_project_docker - < $CONF_project_root_host/$CONF_project_name.Dockerfile
