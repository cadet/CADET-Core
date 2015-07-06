#!/bin/bash

MATLAB_ROOT=/usr/local/MATLAB/R2014a

if [ $# -ne 1 ]
then
    echo "Runs the already built MEX file in Matlab environment using a wrapper."
    echo "Usage: buildRunMEXdebug.sh inputFile.mat"
    echo "   Where inputFile.mat is a MAT file that contains a struct 'data' which holds"
    echo "   the file format structs (i.e., it contains a struct 'input' which corresponds"
    echo "   to the input group of the file format specification)."
    echo ""
    echo "Note: Matlab path is set to "
    echo "      $MATLAB_ROOT" 
    echo "      is this correct? You can change it by opening this file and editing the"
    echo "      MATLAB_ROOT variable accordingly."
    exit 1
fi


# Build wrapper
$MATLAB_ROOT/bin/mex -g -f $MATLAB_ROOT/bin/engopts.sh -ldl -lmat CadetMexDebugLinux.cpp

# Run wrapper with already built MEX file and given input
OMP_NUM_THREADS=1 PATH=$MATLAB_ROOT/bin:$PATH LD_LIBRARY_PATH=$MATLAB_ROOT/bin/glnxa64/:$MATLAB_ROOT/sys/os/glnxa64/:$LD_LIBRARY_PATH ./CadetMexDebugLinux ../src/cadet-mi/CadetMex.mexa64 $1

