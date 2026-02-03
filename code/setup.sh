#!/bin/bash
# This script is the basic setup script for the project
# run from the code directory with the following cmd:
# $ . ./setup.sh
# it will create the environment variables PROJECT_LOC & DATA_OUT in the current shell
cd ..
PROJECT_LOC=$(pwd)
export PROJECT_LOC
cd code
mkdir -p /mnt/speliotes-lab/Projects/UK_ATLAS/IndivProj/${USER}
DATA_OUT=/mnt/speliotes-lab/Projects/UK_ATLAS/IndivProj/${USER}/CRISPRai_NGS/
export DATA_OUT
mkdir -p /mnt/speliotes-lab/Projects/UK_ATLAS/IndivProj/${USER}/CRISPRai_NGS/
