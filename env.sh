#!/bin/bash

if [ $# -ne 1 ] ; then
    echo "Usage: ./env.sh env_name"
    echo "       env_name = { bbtt | hh | lcg }"
    exit 1
fi

if ! [ -z ${ENV_NAME+x} ] ; then
    echo "Another environment is already loaded. Please, exit from the current environment first."
    exit 2
fi

function run_cmd {
    "$@"
    RESULT=$?
    if [ $RESULT -ne 0 ] ; then
        echo "Error while rinning '$@'"
        exit 1
    fi
}

function run_source {
    source setup.sh
    RESULT=$?
    if [ $RESULT -ne 0 ] ; then
        echo "Error while rinning 'source setup.sh'"
        exit 1
    fi
}

export ENV_NAME=$1

if [ $ENV_NAME = "bbtt" ] ; then
    cmssw_ver=CMSSW_10_2_13
    if ! [ -f $cmssw_ver/src/.installed ] ; then
        echo "Installing an environment for the bbtt shape createion..."
        export SCRAM_ARCH=slc7_amd64_gcc700
        scram project CMSSW $cmssw_ver
        run_cmd cd $cmssw_ver/src
        run_cmd eval `scramv1 runtime -sh`
        run_cmd git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
        run_cmd cd HiggsAnalysis/CombinedLimit
        run_cmd git checkout v8.1.0
        run_cmd cd ../..
        run_cmd git clone https://github.com/cms-analysis/CombineHarvester.git CombineHarvester
        run_cmd scram b -j4
        run_cmd touch .installed
    else
        echo "Loading an environment for the bbtt shape createion..."
        run_cmd cd $cmssw_ver/src
        run_cmd eval `scramv1 runtime -sh`
    fi
    run_cmd cd ../..
    export PYTHONPATH=$PWD/$PYTHONPATH
elif [ $ENV_NAME = "hh" ] ; then
    if ! [ -f inference/.installed ] ; then
        echo "Installing an environment for the HH statistical inference..."
        run_cmd git clone --recursive ssh://git@gitlab.cern.ch:7999/hh/tools/inference.git
        run_cmd cd inference
        run_source
        run_cmd touch .installed
    else
        echo "Loading an environment for the HH statistical inference..."
        run_cmd cd inference
        run_source
    fi
    run_cmd cd ..
elif [ $ENV_NAME = "lcg" ] ; then
    source /cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_99 x86_64-centos7-gcc10-opt
else
    echo "Unknown environment '$ENV_NAME'. Supported environments: bbtt hh"
    exit 3
fi

echo "Starting new $SHELL ..."
$SHELL

echo "Back to the original environment."
