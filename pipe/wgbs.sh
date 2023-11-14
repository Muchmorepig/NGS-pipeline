#!/usr/bin/bash

script_path=$(dirname $(readlink -f $0))

source ${script_path}/../src/cmd_check.sh

command_exists trim_galore


source ${script_path}/../src/trim_galore_cmd.sh

galore ./rawdata clean