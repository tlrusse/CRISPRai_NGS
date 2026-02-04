#!/bin/bash
# This script is the basic setup script for the project
# run from the code directory with the following cmd:
# $ . ./setup.sh
# it will create the environment variables PROJECT_LOC & DATA_OUT in the current shell

PROJ_ORIGINATOR="tlrusse"
export PROJ_ORIGINATOR

cd ..
PROJECT_LOC=$(pwd)
export PROJECT_LOC

DATA_IN="/nfs/corenfs/INTMED-Speliotes-data/Projects/15107-PP/"
export DATA_IN
mkdir -p /mnt/speliotes-lab/Projects/UK_ATLAS/IndivProj/${PROJ_ORIGINATOR}
DATA_OUT=/mnt/speliotes-lab/Projects/UK_ATLAS/IndivProj/${PROJ_ORIGINATOR}/CRISPRai_NGS/
export DATA_OUT
mkdir -p /mnt/speliotes-lab/Projects/UK_ATLAS/IndivProj/${PROJ_ORIGINATOR}/CRISPRai_NGS/


# if the following symlinks do not exist, create them
DATA_IN_LINK=./data_in
if [ ! -L $DATA_IN_LINK ]; then
    ln -s $DATA_IN $DATA_IN_LINK
fi
DATA_OUT_LINK=./data_out
if [ ! -L $DATA_OUT_LINK ]; then
    ln -s $DATA_OUT $DATA_OUT_LINK
fi

cd code