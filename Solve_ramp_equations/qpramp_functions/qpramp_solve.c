#include "mex.h"
#include "qpramp.h"

static int initialized = 0;

static void cleanup(void) {
    qpramp_cleanup();
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if(nrhs!=5) {
        mexErrMsgTxt("Five inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgTxt("One output required.");
    }

    for (int i = 0; i < nrhs; ++i) {
        if( !mxIsDouble(prhs[i]) || 
             mxIsComplex(prhs[i])) {
            mexErrMsgTxt("Input matrix must be type double.");
        }
    }

    // Note that transposed matrices are needed because of MATLABs column-major memory layout.
    size_t c = mxGetM(prhs[0]);
    size_t n = mxGetM(prhs[1]);
    size_t p = mxGetM(prhs[3]);

    if(mxGetM(prhs[0])!=c || mxGetN(prhs[0])!=c) {
        mexErrMsgTxt("Input neg_g_invh_gt must be c by c.");
    }
    if(mxGetM(prhs[1])!=n || mxGetN(prhs[1])!=c) {
        mexErrMsgTxt("Input neg_s must be n by c.");
    }
    if(mxGetM(prhs[2])!=1 || mxGetN(prhs[2])!=c) {
        mexErrMsgTxt("Input neg_w must be 1 by c.");
    }
    if(mxGetM(prhs[3])!=p || mxGetN(prhs[3])!=c) {
        mexErrMsgTxt("Input neg_g_invh must be p by c.");
    }
    if(mxGetM(prhs[4])!=1 || mxGetN(prhs[4])!=n) {
        mexErrMsgTxt("Input x must be 1 by n.");
    }

    double *neg_g_invh_gt = mxGetPr(prhs[0]);
    double *neg_s = mxGetPr(prhs[1]);
    double *neg_w = mxGetPr(prhs[2]);
    double *neg_g_invh = mxGetPr(prhs[3]);
    double *x = mxGetPr(prhs[4]);

    plhs[0] = mxCreateDoubleMatrix(1,(mwSize)p,mxREAL);
    double *z = mxGetPr(plhs[0]);

    if (!initialized) {
        qpramp_init(c, p);
        qpramp_enable_infeasibility_error(1e-12, 1e12); // Remove comment to enable infeasibility errors
        initialized = 1;
    }
    int err = qpramp_solve(c, n, p, (double(*)[])neg_g_invh_gt, (double(*)[])neg_s, neg_w, (double(*)[])neg_g_invh, x, z);
    switch (err) {
        case (QPRAMP_ERROR_INFEASIBLE):
            mexErrMsgTxt("Problem is infeasible.");
            break;
        case (QPRAMP_ERROR_RANK_2_UPDATE):
            mexErrMsgTxt("Unable to perform rank 2 update.");
            break;
        default:
            break;
    }
    mexAtExit(cleanup);
}
 
