#include "mex.h"
#include "raiseerror.h"

void raiseerror_api(const char *msg)
{
    mexErrMsgTxt(msg);
};