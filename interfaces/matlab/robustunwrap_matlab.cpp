/*
 This is the seperated Matlab MEX interface for

 3D Unwrapping algorithm
 Rhodri Cusack 2000-2006
 Algorithm described in 
 Cusack, R. & Papadakis, N. (2002) 
    "New robust 3-D phase unwrapping algorithms: application to magnetic field mapping and undistorting echoplanar images."
    Neuroimage. 2002 Jul;16(3 Pt 1):754-64.
 Distributed under MIT License
 Comments to cusackrh@tcd.ie
 Version 3.00 Oct 2006 adapted for matlab

 Refactoring done by Podranski, K. (Sep 2023)

 Compile within Matlab with:
 mex robustunwrap_matlab.cpp robustunwrap.cpp raiseerror_matlab.cpp -output robustunwrap
 */



#include "mex.h"
#include "matrix.h"
#include "../../include/robustunwrap.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    const mwSize ndims = 3;

    /* Check for proper number of arguments. */
    if (nrhs < 3) {
        mexErrMsgTxt("Robust unwrapping algorithm as in Cusack & Papadakis (2002).\nExpect at least 3 input arguments, unwrappeddata=unwrap([seed coords, 3 vector],[phase data, 3D double array],[magnitude data, 3D double array],{[numunwrapbins]}) ");
    }
    if (nlhs > 1) {
        mexErrMsgTxt("Too many output arguments");
    }

    if (mxGetNumberOfDimensions(prhs[0]) != 2 || mxGetDimensions(prhs[0])[0] != 1 || mxGetDimensions(prhs[0])[1] != 3) {
        mexErrMsgTxt("Seed coords must be line vector with 3 elements.");
    }
    if (mxIsComplex(prhs[0])) {
        mexErrMsgTxt("Voxel dimensions of seed coords should not be complex.");
    }

    /* The input must be a noncomplex double 3D matrix. */
    {
        mwSize ndims_phs = mxGetNumberOfDimensions(prhs[1]);
        if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || ndims_phs != ndims) {
            mexErrMsgTxt("Data to be unwrapped must be a double 3D matrix.");
        }
    }
    {
        mwSize ndims_mag = mxGetNumberOfDimensions(prhs[2]);
        if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) || ndims_mag != ndims) {
            mexErrMsgTxt("Data to be unwrapped must be a double 3D matrix.");
        }
    }

    mxDouble *voxdims = mxGetPr(prhs[0]);
    ptrdiff_t seedx = ptrdiff_t(voxdims[0]) - 1;
    ptrdiff_t seedy = ptrdiff_t(voxdims[1]) - 1;
    ptrdiff_t seedz = ptrdiff_t(voxdims[2]) - 1;

    ptrdiff_t dim[ndims];
    const mwSize *dims = mxGetDimensions(prhs[1]);
    const mwSize *dims_mag = mxGetDimensions(prhs[2]);
    size_t sze = 1;
    for (size_t i = 0; i < ndims; i++) {
        if (dims[i] != dims_mag[i]) {
            mexErrMsgTxt("Phase data and magnitude data should have the same size.");
        }
        sze *= dims[i];
        dim[i] = dims[i];
    };

    /* define strides */
    ptrdiff_t m_bsx = 1;
    ptrdiff_t m_bsy = dims[0];
    ptrdiff_t m_bsz = dims[0] * dims[1];

    if (seedx < 0 || seedx >= dims[0] || seedy < 0 || seedy >= dims[1] || seedz < 0 || seedz >= dims[2]) {
        mexErrMsgTxt("The seed specified was outside the matrix bounds.");
    }

    mxDouble *phase = mxGetPr(prhs[1]); // TODO: eventually replace with mxGetDoubles

    /* Negate input as low polefield values unwrapped first */
    mxDouble *maginput = mxGetPr(prhs[2]); // TODO: eventually replace with mxGetDoubles
    mxDouble *mag = new mxDouble[sze];
    for(size_t i = 0; i < sze; i++)
        mag[i] = -maginput[i];

    int numunwrapbins = 10000;
    if (nrhs == 4) {
        if (!mxIsDouble(prhs[3]) ||  mxIsComplex(prhs[3]))
            mexErrMsgTxt("Number of unwrapping bins should be a non-complex double.");
        numunwrapbins = int(*mxGetPr(prhs[3]));  // TODO: eventually replace with mxGetDoubles
    }

    /* Create matrix for the return argument. */
    plhs[0] = mxCreateNumericArray(ndims, dims, mxDOUBLE_CLASS, mxREAL);
    mxDouble *unwrapped = mxGetPr(plhs[0]);

    /* Assign pointers to each input and output. */
    unwrap_helper(seedx, seedy, seedz, numunwrapbins, dim, sze, m_bsx, m_bsy, m_bsz, phase, mag, unwrapped);

    /* cleanup */
    delete mag;
}
