#!/bin/bash

# Source all the GENIE stuff
source ~/ssd/lamp/environment_setup.sh
export LIBXML2_LIB=/usr/lib64
export LIBXML2_INC=/usr/include


# Source the NEUT stuff
export NEUT_ROOT=/vols/build/t2k/cvw09/neut_5.3.3_Autumn2016
export CERN=/vols/build/t2k/cvw09/CERNLIB
export CERN_LEVEL=2005
export LD_LIBRARY_PATH=${NEUT_ROOT}/src/reweight:${LD_LIBRARY_PATH}



export NUISANCE=$(pwd -P)
