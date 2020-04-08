#! /bin/bash

blastnoutput=$1
outputfoler=$2

pickotus="/opt/miniconda3/envs/qiime1/bin/pick_otus.py"

source /opt/miniconda3/envs/qiime1/bin/activate qiime1

export LANG=en_US.UTF-8;

"$pickotus" -i "$blastnoutput" -o "$outputfoler" -s 0.97;

source /opt/miniconda3/envs/qiime1/bin/deactivate
