#!/bin/bash

if [ $# -eq 0 ]
  then
    echo "Argument 1 should be a port number, argument 2 a string to use as password."
    exit 1
fi

# parse config parameters:
source bash/parse_yaml.sh
eval $(parse_yaml config.yaml CONF_)

docker run -d --rm \
  --name ${CONF_project_name}_${USER} \
  -p ${1}:8787 \
  -e PASSWORD=${2} \
  -e USERID=$(id -u) \
  -e GROUPID=$(id -g) \
  -e UMASK=002 \
  -e TZ=Europe/Vienna \
  --volume=$(realpath ${CONF_project_root_host%%*##*( )}):${CONF_project_root} \
  --volume=$(realpath ${CONF_data_root_host%%*##*( )}):${CONF_data_root}:ro \
  --volume=$(realpath ${CONF_out_root_host%%*##*( )}):${CONF_out_root} \
  --volume=$(realpath ${CONF_tmp_root_host%%*##*( )}):${CONF_tmp_root} \
  --volume=$(realpath ${CONF_oldsinglecell_root_host%%*##*( )}):${CONF_oldsinglecell_root} \
  --volume=$(realpath ${CONF_resource_root_host%%*##*( )}):${CONF_resource_root}:ro \
  ${CONF_project_docker}
