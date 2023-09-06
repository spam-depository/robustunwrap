#include <stdio.h>
#include <string.h>
#include "math.h"
#include "getopt.h"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

/*
 3D Unwrapping algorithm
 Rhodri Cusack 2000-2006
 Algorithm described in 
 Cusack, R. & Papadakis, N. (2002) 
    "New robust 3-D phase unwrapping algorithms: application to magnetic field mapping and undistorting echoplanar images."
    Neuroimage. 2002 Jul;16(3 Pt 1):754-64.
 Distributed under MIT License
 Comments to cusackrh@tcd.ie
 Version 3.00 Oct 2006 adapted for matlab
 This version adapted for Python/Numpy by K. Podranski Sep 2023
 
 manual compilation with
 c++ -O3 -Wall -shared -std=c++11 -fPIC $(python3 -m pybind11 --includes) robustunwrap_python.cpp -o robustunwrap$(python3-config --extension-suffix)
 */

#define PI 3.14159265358979
#define FILEFORMAT_CAMRES 1
#define FILEFORMAT_ANALYZE 2
#define MAXFILENAMESIZE 1024

int errcode;
long dim[3];
long sze;
long polesze;
const double *phase;
double *mag;
double *unwrapped;
bool *flag;
long m_bsx,m_bsy,m_bsz;

struct FIELD
{
    double d;	/* value of field */
    long p[3];  /* at this offset */
    int x,y,z;
};

struct QUEUEENTRY
{
    int x,y,z;
    long p;
    double v;
};

struct FIELD_2
{
    double d;	/*value of field*/
    long p;		/*at this offset*/
    int x,y,z;
};

void raiseerror(std::string msg)
{
    //mexErrMsgTxt(msg);
    //errcode=1;
    throw std::runtime_error(msg);
};

const long NUMQUEUES=10000;
const long DEFAULTBLOCK=100;
const long BLOCKINCREMENT=500;

long m_bot[NUMQUEUES], m_top[NUMQUEUES];
long m_size[NUMQUEUES], m_chunk[NUMQUEUES], m_sizeb[NUMQUEUES];
char *m_q[NUMQUEUES];

void InitQueue(int queuenum,int chunk)
{
    m_chunk[queuenum]=chunk;
    m_bot[queuenum]=0;
    m_top[queuenum]=0;
    m_size[queuenum]=DEFAULTBLOCK;
    m_sizeb[queuenum]=m_size[queuenum]*m_chunk[queuenum];
    
    m_q[queuenum]=0;
};

void TerminateQueue(int queuenum)
{
    if (m_q[queuenum]) delete m_q[queuenum];
    m_q[queuenum]=0;
};

int Push(int queuenum,void *x)
{
    if (!m_q[queuenum]) m_q[queuenum]=new char[m_sizeb[queuenum]];
    if (!m_q[queuenum])
    {
        printf("Out of memory - could not generate new point queue.");
        return(-1);
    };
    
    int stacksize;
    stacksize=m_top[queuenum]-m_bot[queuenum];
    if (stacksize<0) stacksize+=m_size[queuenum];
    memcpy(m_q[queuenum]+m_top[queuenum]*m_chunk[queuenum],x,m_chunk[queuenum]);
    m_top[queuenum]++;
    
    if (m_top[queuenum]==m_size[queuenum]) m_top[queuenum]=0;
    
    if (m_top[queuenum]==m_bot[queuenum])
    {
        char *newq;
        long newsize, newsizeb, abovebot_b;
        // previously out of memory - now auto-expand
        newsize=m_size[queuenum]+BLOCKINCREMENT;
        newsizeb=newsize*m_chunk[queuenum];
        newq=new char[newsizeb];
        if (!newq)
        {
            printf("Out of memory - point queue full.");
            return(-1);
        };
        // While we're shifting circular buffer, it is actually easier to re-origin it
        // to zero.
        // first, copy top bit from m_bot upwards
        abovebot_b=(m_size[queuenum]-m_bot[queuenum])*m_chunk[queuenum];
        if (abovebot_b>0) memcpy(newq,(char *)m_q[queuenum]+m_bot[queuenum]*m_chunk[queuenum],abovebot_b);
        // then, do m_top downwards
        if (m_top[queuenum]!=0) memcpy((char *) newq+abovebot_b,m_q[queuenum],m_top[queuenum]*m_chunk[queuenum]);
        m_bot[queuenum]=0;
        m_top[queuenum]=m_size[queuenum];
        m_size[queuenum]=newsize;
        m_sizeb[queuenum]=newsizeb;
        delete m_q[queuenum]; // recover old memory
        m_q[queuenum]=newq;
    };
    
    return(0);
    
};

int Pop(int queuenum,void *x)
{
    if (m_bot[queuenum]==m_top[queuenum]) return(1);
    
    memcpy(x,m_q[queuenum]+m_bot[queuenum]*m_chunk[queuenum],m_chunk[queuenum]);
    m_bot[queuenum]++;
    if (m_bot[queuenum]==m_size[queuenum]) m_bot[queuenum]=0;
    
    return(0);
};

void Check(int primaryqueuenum, QUEUEENTRY* qe,long offp, int offx, int offy, int offz)
{
    /* first check bounds */
    QUEUEENTRY nqe;
    nqe.x=qe->x+offx;
    if ((nqe.x<0)||(nqe.x>=dim[0])) return;
    nqe.y=qe->y+offy;
    if ((nqe.y<0)||(nqe.y>=dim[1])) return;
    nqe.z=qe->z+offz;
    if ((nqe.z<0)||(nqe.z>=dim[2])) return;
    
    
    nqe.p=qe->p+offp;
    
    if (flag[nqe.p]) return; /* Already been here */
    
    
    /* Actually do unwrap */
    int wholepis;
    wholepis=int((phase[nqe.p]-qe->v)/PI);
    
    if (wholepis>=1)
        nqe.v=(double) phase[nqe.p]-2*PI*int((wholepis+1)/2);
    else if (wholepis<=-1)
        nqe.v=(double) phase[nqe.p]+2*PI*int((1-wholepis)/2);
    else nqe.v=phase[nqe.p];
    
    unwrapped[nqe.p]=nqe.v;
    flag[nqe.p]=true;
    
    if (Push(primaryqueuenum,&nqe))
        raiseerror("Out of memory!");
    
    
};

void unwrap(int seedx, int seedy, int seedz,long UNWRAPBINS)
{
    long i;
    long seedp;
    
    errcode=0;
    
    /* Minimum number of unwrapping bins is 2 */
    if (UNWRAPBINS<2) UNWRAPBINS=2;
    /* Find min and max */
    double min,max;
    min=(double) 1e38;
    max=0;
    
    for (i=0; i<sze; i++)
    {
        min=mag[i]<min?mag[i]:min;
        max=mag[i]>max?mag[i]:max;
    };
    
    
    double diff;
    diff=(double) 1.00001*(max-min);
    
    seedp=seedx+dim[0]*(seedy+dim[1]*seedz);
    
    flag=new bool[sze];
    
    if (!flag) raiseerror("Out of memory. ");
    for (i=0; i<sze; i++)
        flag[i]=false;
    
    
    for (i=0; i<UNWRAPBINS; i++)
        InitQueue(i,sizeof(QUEUEENTRY));
    QUEUEENTRY qe;
    
    qe.p=seedp;
    qe.x=seedx;
    qe.y=seedy;
    qe.z=seedz;
    unwrapped[85]=phase[85];
    unwrapped[seedp]=qe.v=phase[seedp];
    flag[seedp]=true;
    
/* push seed */
    Push(0,&qe);
    
    
    //int opc,pc;
    //opc=-1;
    
/* First, work out pole field threshold that we're going to use */
    double *polefieldthresholds;
    int ind=0;
    int voxdone,voxdef;
    
    polefieldthresholds=new double [UNWRAPBINS];
    voxdef=voxdone=0;
    
    for (i=0; i<(UNWRAPBINS); i++)
        polefieldthresholds[ind++]=min+diff*i/(UNWRAPBINS-1); 
    //opc=-10;
    
    for (i=0; i<UNWRAPBINS ; i++)
    {
        voxdef=0;
        while (!errcode && !Pop(i,&qe))
        {
            if (mag[qe.p]>polefieldthresholds[i])  /* too close to a scary pole, so just defer by pushing to other stack */
            {
                ind=i;
                while (mag[qe.p]>polefieldthresholds[++ind]);
                Push(ind,&qe);										/* just defer by pushing to relevant stack */
                voxdef++;
            }
            else
            {
                Check(i,&qe,+m_bsz,0,0,1);
                Check(i,&qe,-m_bsz,0,0,-1);
                Check(i,&qe,+m_bsy,0,1,0);
                Check(i,&qe,-m_bsy,0,-1,0);
                Check(i,&qe,+m_bsx,1,0,0);
                Check(i,&qe,-m_bsx,-1,0,0);
                voxdone++;
            };
        };
        TerminateQueue(i);	/* done with this Queue so free up memory */
        
        if (errcode) break;
    };
    
    delete flag;
};


py::array_t<double, py::array::f_style | py::array::forcecast> robustunwrap(py::list seed, py::array_t<double, py::array::f_style | py::array::forcecast> phaseData, py::array_t<double, py::array::f_style | py::array::forcecast> magnitudeData)
{
    //long i,j,k;
    int seedx,seedy,seedz;
    //double *voxdims;
    //int ndims;
    int numunwrapbins;

    if (seed.size() != 3) {
        throw std::runtime_error("seed must be an index into the 3D data array consisting of 3 ints.");
    }
    
    const long int *dims;
    //double *maginput;
    dims = magnitudeData.shape();

    if (magnitudeData.ndim() != 3 || phaseData.ndim() != 3) {
        throw std::runtime_error("Number of dimensions must be three");
    //} else if (magnitudeData.shape() != phaseData.shape()) {
    //} else if (!std::equal(std::begin(magnitudeData.shape()), std::end(magnitudeData.shape()), std::begin(phaseData.shape()))) {
    } else if (!std::equal(dims, dims + 3, phaseData.shape())) {
        throw std::runtime_error("magnitudeData and phaseData shapes must match");
    }
    // else if (magnitudeData.dtype != ??? || phaseData.dtype != ???) { ...

    /* The input must be a noncomplex scalar double.*/
    //ndims = magnitudeData.ndim();

    //global
    sze=1;
    for (long i=0; i<3; i++)
    {
        sze *= dims[i];
        dim[i] = dims[i];
    };
    
  /* Create matrix for the return argument. */
    //plhs[0] = mxCreateNumericArray(ndims,dims, mxDOUBLE_CLASS,mxREAL);
    //auto unwrappedArray = py::array_t<double>(phaseData.size());
    auto unwrappedArray = py::array_t<double, py::array::f_style | py::array::forcecast>({dim[0], dim[1], dim[2]}); // not very elegant

    //voxdims=mxGetPr(prhs[0]);
    seedx = py::cast<int>(seed[0]); // maybe a bit naive. might not work that way.
    seedy = py::cast<int>(seed[1]);
    seedz = py::cast<int>(seed[2]);
    
    if (seedx<0 || seedx>=dims[0] || seedy<0 || seedy>=dims[1] || seedz<0 || seedz>=dims[2])
        throw std::runtime_error("The seed specified was outside the matrix bounds.");
    
    //double *phase = static_cast<double *>(phaseData.ptr);
    //const double *phase = phaseData.data();
    phase = phaseData.data();
    
    /* Negate input as low polefield values unwrapped first */
    //double *maginput = static_cast<double *>(magnitudeData.ptr);
    const double *maginput = magnitudeData.data();
    mag=new double[sze];
    for(long i=0; i<sze; i++) // along the flattened ndarray
        mag[i]=-maginput[i];
    
    //double *unwrapped = static_cast<double *>(unwrappedArray.request().ptr);
    unwrapped = static_cast<double *>(unwrappedArray.request().ptr);
    /* TODO
    if (nrhs==4)
    {
        if (!mxIsDouble(prhs[3]) ||  mxIsComplex(prhs[3]))
        {
            mexErrMsgTxt("Number of unwrapping bins should be a non-complex double.");
        }
        numunwrapbins=int(*mxGetPr(prhs[3]));
    }
    else
    */
        numunwrapbins=10000;
    m_bsx=1;
    m_bsy=dims[0];
    m_bsz=dims[0]*dims[1];
  /* Assign pointers to each input and output. */
    unwrap(seedx,seedy,seedz,numunwrapbins);
    
    return unwrappedArray;
}

py::array_t<double> add_arrays(py::array_t<double> input1, py::array_t<double> input2) {
    py::buffer_info buf1 = input1.request(), buf2 = input2.request();

    if (buf1.ndim != 1 || buf2.ndim != 1)
        throw std::runtime_error("Number of dimensions must be one");

    if (buf1.size != buf2.size)
        throw std::runtime_error("Input shapes must match");

    /* No pointer is passed, so NumPy will allocate the buffer */
    auto result = py::array_t<double>(buf1.size);

    py::buffer_info buf3 = result.request();

    double *ptr1 = static_cast<double *>(buf1.ptr);
    double *ptr2 = static_cast<double *>(buf2.ptr);
    double *ptr3 = static_cast<double *>(buf3.ptr);

    for (size_t idx = 0; idx < (unsigned long) buf1.shape[0]; idx++)
        ptr3[idx] = ptr1[idx] + ptr2[idx];

    return result;
}

void show_arr(py::array_t<double, py::array::f_style | py::array::forcecast> data) {

    int ndim = data.ndim();
    const long int *shp = data.shape();
    const double *dataptr = data.data();
    py::print(ndim);
    py::print(shp);
    for(long i=0; i<data.size(); i++) {
        py::print(dataptr + i);
    }
}

PYBIND11_MODULE(robustunwrap, m) {
    m.doc() = "attempt to interface robustunwrap (Cusack & Papadakis 2002) with pybind";
    m.def("add_arrays", &add_arrays, "Add two NumPy arrays");
    m.def("show_arr", &show_arr, "print array shape and content");
    m.def("robustunwrap", &robustunwrap, "robustunwrap (Cusack & Papadakis 2002)");
}
