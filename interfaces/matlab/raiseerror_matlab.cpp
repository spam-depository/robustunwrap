#include "mex.h"
#include "raiseerror.h"

void raiseerror(const char *msg)
{
    mexErrMsgTxt(msg);
};