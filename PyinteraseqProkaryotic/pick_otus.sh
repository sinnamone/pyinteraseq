#! /bin/bash
source /opt/miniconda3/envs/qiime1/bin/activate /opt/miniconda3/envs/qiime1

blastnoutput=$1
outputfoler=$2

pickotus="/opt/miniconda3/envs/qiime1/bin/pick_otus.py"
"$pickotus" -i "$blastnoutput" -o "$outputfoler" -s 0.97