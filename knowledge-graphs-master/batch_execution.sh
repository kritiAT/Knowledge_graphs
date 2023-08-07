#!/bin/bash
# Execute avg. graph script for all features in one run
# 
# Usage:
# bash batch_execution.sh 220329_feature_list.txt 
#

# Maybe uncomment
# conda activate <env-name>

while read p; do
	python generate_graphs.py akgs 5000 ${p}_graph $p -t 25
done < $1
