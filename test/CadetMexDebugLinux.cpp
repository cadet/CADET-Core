// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2008-2015: Eric von Lieres¹, Joel Andersson¹,
//                         Andreas Puettmann¹, Sebastian Schnittert¹,
//                         Samuel Leweke¹
//                                      
//    ¹ Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

// This code is partly taken from angainor's answer to the question "How do I profile a MEX-function in Matlab" 
// (http://stackoverflow.com/questions/11220250/how-do-i-profile-a-mex-function-in-matlab) on Stack Overflow.
// See http://stackoverflow.com/a/12405131 for the answer containing the code.
// 
// In order to avoid the Matlab engine, the input data is read using Matlab's MAT interface.
// This very interface is also used to save the computed results.
//
// Command for compiling (assuming Matlab R2014a is installed in the default path):
//
// /usr/local/MATLAB/R2014a/bin/mex -g -f /usr/local/MATLAB/R2014a/bin/engopts.sh -ldl -lmat CadetMexDebugLinux.cpp
//
// Commands for running:
// 
// export OMP_NUM_THREADS=1       optional, overrides thread settings given in data
// export PATH=/usr/local/MATLAB/R2014a/bin:$PATH 
// export LD_LIBRARY_PATH=/usr/local/MATLAB/R2014a/bin/glnxa64/:/usr/local/MATLAB/R2014a/sys/os/glnxa64/:$LD_LIBRARY_PATH 
// ./CadetMexDebugLinux /path/to/CadetMex.mexa64 /path/to/Input.mat
// 
// The input MAT file contains the "data" struct as created by run(obj, task) function
// in the Simulator.m class, i.e., it contains the input group specified in the file format.
// The result is written to the file "result-XYZ.mat", where XYZ denotes the PID of the
// CadetMexDebugLinux process.

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <dlfcn.h>
#include <unistd.h>
#include "mat.h"

typedef void (*mexFunction_t)(int nargout, mxArray *pargout [ ], int nargin, const mxArray *pargin[]);

int main(int argc, const char *argv[])
{
    printf("Start\n");

    // Check input args
    if (argc < 3)
    {
        fprintf(stderr, "Error. Give full path to the MEX and the input file as input parameters.\n");
        return -1;
    }

    // Open MEX file
    void *handle = dlopen(argv[1], RTLD_NOW);
    if (!handle)
    {
        fprintf(stderr, "A dynamic linking error occurred: (%s)\n", dlerror());
        fprintf(stderr, "Error loading MEX file: %s\n", strerror(errno));
        return -1;
    }

    // Load MEX function handle from shared object
    mexFunction_t mexfunction = (mexFunction_t)dlsym(handle, "mexFunction");
    if (!mexfunction)
    {
        fprintf(stderr, "MEX file does not contain mexFunction\n");
        return -1;
    }

    printf("Load input\n");

    // Load input from MAT-file
    MATFile* pmat = matOpen(argv[2], "r");
    mxArray* arg1 = matGetVariable(pmat, "data");
    mxArray *pargout[1] = {0};
    const mxArray *pargin[1] = {arg1};

    // Call MEX function
    printf("Execute MEX\n");
    mexfunction(1, pargout, 1, pargin);

    char saveString[50];
    snprintf(saveString, 50, "result-%d.mat", getpid());

    printf("Save results\n");
    MATFile* resmat = matOpen(saveString, "w");
    matPutVariable(resmat, "data", pargout[0]);
    matClose(resmat);
    matClose(pmat);

    printf("Finished\n");

    // Clean up
    mxDestroyArray(pargout[0]);
    dlclose(handle);

    return 0;
}
