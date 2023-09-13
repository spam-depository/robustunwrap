#include "mex.h"
#include "../../include/raiseerror.h"

void raiseerror_api(const char *msg)
{
    mexErrMsgTxt(msg);
};