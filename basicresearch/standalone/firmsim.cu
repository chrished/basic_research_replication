#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <cuda.h>
#include <cublas.h>
#include <mex.h>

#include "mersennetwister/MersenneTwister.h"
#include "mersennetwister/MersenneTwister_kernel.cu"

using namespace std;

///////////////////////////////////////////////////////////////////////////////
// Constants
///////////////////////////////////////////////////////////////////////////////
//#define nF 16384
#define nF 32768
//#define nF 65536

#define R_BURN_INIT 10000
#define R_BURN_MED 2000

#define T_PERIODS 1

#define RAND_CHUNKS 5
#define N_QDISTS 3

#define nN 64
#define nM 10

#define nS 8
#define nSpow 3
#define STATE_NOTHING      0
#define STATE_LOSE         1
#define STATE_GAIN_APP     2
#define STATE_GAIN_BAS_1   3
#define STATE_GAIN_BAS_2   4
#define STATE_GAIN_FREE    5
#define STATE_EXPAND       6
#define STATE_EXIT         7

#define BLOCK_SIZE_F 16

///////////////////////////////////////////////////////////////////////////////
// Common host and device function 
///////////////////////////////////////////////////////////////////////////////

//ceil(a / b)
extern "C" int iDivUp(int a, int b){
    return ((a % b) != 0) ? (a / b + 1) : (a / b);
}

//floor(a / b)
extern "C" int iDivDown(int a, int b){
    return a / b;
}

//Align a to nearest higher multiple of b
extern "C" int iAlignUp(int a, int b){
    return ((a % b) != 0) ?  (a - a % b + b) : a;
}

//Align a to nearest lower multiple of b
extern "C" int iAlignDown(int a, int b){
    return a - a % b;
}

void Check_CUDA_Error(const char *message)
{
   cudaError_t error = cudaGetLastError();
   if(error!=cudaSuccess) {
      fprintf(stderr,"ERROR: %s: %s\n", message, cudaGetErrorString(error) );
      //exit(-1);
   }
}

///////////////////////////////////////////////////////////////////////////////
// Device kernels
///////////////////////////////////////////////////////////////////////////////

// binary search (only powers of 2 sizes)
__device__ inline int dsamp2(double* pbeg, int pow2, double r)
{
  int len = 1<<(pow2-1);
  int pos = len-1;

  double val;
  for (int s = 0; s < (pow2-1); s++) {
    len /= 2;
    val = pbeg[pos];
    pos += (r < val) ? -len : len;
  }

  val = pbeg[pos];
  pos += (r > val) ? 1 : 0;

  return pos;
}

__global__ void fqsim(double* fquals, double qdec, double qmin)
{
  int q = blockDim.x*blockIdx.x + threadIdx.x;
  int f = blockDim.y*blockIdx.y + threadIdx.y;

  if ((f < nF) && (q < nN)) {
    int i = f*nN+q;
    double qval = fquals[i]*qdec;
    fquals[i] = (qval > qmin) ? qval : qmin;
  } 
}

__global__ void fzsim(int* fnprod, int* fmind, int* fage, double* fquals, int* fpinds, int* fexited,
                      double* qdists, double* qbins, double* tvecs, double* fnrand, int nBpow, int nB, int primetime)
{
  int f = blockDim.x*blockIdx.x + threadIdx.x;

  if (f < nF) {
    int n = fnprod[f];
    int m = fmind[f];
    int age = fage[f];

    int exited = 0;
    if (primetime == 1) {
      exited = fexited[f];
    }

    int* pinds = fpinds + nN*f;
    double* quals = fquals + nN*f;

    int nind;
    int qind;
    int mind;

    // requests
    double* qreq1 = NULL;
    double* qreq2 = NULL;

    // randoms
    double r0 = fnrand[f];
    double r1 = fnrand[nF+f];
    double r2 = fnrand[2*nF+f];
    double r3 = fnrand[3*nF+f];
    double r4 = fnrand[4*nF+f];

    // the sampler
    double* tvec = tvecs + (m-1)*(nN+1)*nS + n*nS;
    int s = dsamp2(tvec,nSpow,r0);
    switch (s) {
      case STATE_EXIT:
        n = 0;
        break;
      case STATE_GAIN_APP:
        qreq1 = qdists;
        break;
      case STATE_GAIN_FREE:
        qreq1 = qdists+2*nB;
        break;
      case STATE_GAIN_BAS_1:
        qreq1 = qdists+nB;
        break;
      case STATE_GAIN_BAS_2:
        qreq1 = qdists+nB;
        qreq2 = qdists+nB;
        break;
      case STATE_LOSE:
        if (n > 1) {
          nind = floor(r1*n);

          quals[nind] = quals[n-1];
          pinds[nind] = pinds[n-1];
        }
        n--;
        break;
      case STATE_EXPAND:
        if (m < nM) {
          qind = dsamp2(qdists,nBpow,r1);

          quals[n] = qbins[qind];
          pinds[n] = m;

          m++;
          n++;
        }
        break;
    }

    age++;

    if (n == 0) {
      exited = 1;
      n = 1;
      m = 1;
      age = 0;

      qind = dsamp2(qdists,nBpow,r3);

      quals[0] = qbins[qind];
      pinds[0] = 0;
    }

    if ((qreq1 != NULL) && (n < nN)) {
      qind = dsamp2(qreq1,nBpow,r1);
      mind = floor(r3*m);

      quals[n] = qbins[qind];
      pinds[n] = mind;

      n++;
    }

    if ((qreq2 != NULL) && (n < nN)) {
      qind = dsamp2(qreq2,nBpow,r2);
      mind = floor(r4*m);

      quals[n] = qbins[qind];
      pinds[n] = mind;

      n++;
    }

    fnprod[f] = n;
    fmind[f] = m;
    fage[f] = age;

    if (primetime == 1) {
      fexited[f] = exited;
    }
  }
}

__global__ void calc_qpow(int* fnprod, double* fquals, double* fqpow1, double pow1, double qmin)
{
  int f = blockDim.x*blockIdx.x + threadIdx.x;

  if (f < nF) {
    int n = fnprod[f];
    double* quals = fquals + nN*f;

    double qp1 = 0.0;
    for (int i = 0; i < n; i++) {
      qp1 += powf(quals[i],pow1);
    }

    fqpow1[f] = qp1;
  }
}

__global__ void calc_mpos(int* fnprod, int* fpinds, int* fmpos)
{
  int f = blockDim.x*blockIdx.x + threadIdx.x;

  if (f < nF) {
    int n = fnprod[f];
    int* pinds = fpinds + nN*f;

    int mind;
    int mpres[nM] = {0,0,0,0,0,0,0,0,0,0};
    for (int i = 0; i < n; i++) {
      mind = pinds[i];
      mpres[mind] = 1;
    }

    int npos = 0;
    for (mind = 0; mind < nM; mind++) {
      if (mpres[mind] == 1) npos++;
    }

    fmpos[f] = npos;
  }
}

__global__ void sampn(int* fnprod, int* fmind, int* fage, double* fquals, int* fpinds,
                      double* qdist, double* qbins, double* fnrand, int nBpow, int nB)
{
  int f = blockDim.x*blockIdx.x + threadIdx.x;

  if (f < nF) {
    int* pinds = fpinds + nM*f;
    double* quals = fquals + nN*f;

    int qind = dsamp2(qdist+2*nB,nBpow,fnrand[f]);
    double qval = qbins[qind];

    fnprod[f] = 1;
    fmind[f] = 1;
    fage[f] = 0;
    quals[0] = qval;
    pinds[0] = 0;
  }
}

///////////////////////////////////////////////////////////////////////////////
// persistent memory
///////////////////////////////////////////////////////////////////////////////

// Memory sizes
static int N_PER_RNG = 0;
static int RAND_N = 0;

static int nBpow = 0;
static int nB = 0;

static int R_BURN = 0;
static int R_SIM = 0;
static int R_TOT = 0;
static int R_PER_T = 0;

static int nel_ft = 0;
static int nel_fq = 0;
static int nel_fm = 0;
static int nel_em = 0;

static size_t fqsize = 0;
static size_t fmsize = 0;
static size_t ftsize = 0;
static size_t nrsize = 0;
static size_t qdsize = 0;
static size_t qbsize = 0;
static size_t ttsize = 0;
static size_t emsize = 0;

// Host memory
static double* h_ttable = NULL;

// Device memory
static int* d_fnprod = NULL;
static int* d_fmind = NULL;
static int* d_fpinds = NULL;
static int* d_fage = NULL;
static double* d_fqpow1 = NULL;
static double* d_fquals = NULL;
static int* d_fexited = NULL;

static int* d_fmpos = NULL;
static double* d_fnrand = NULL;
static double* d_qdists = NULL;
static double* d_qbins = NULL;
static double* d_ttable = NULL;

// Flags
static int initialized = 0;

// Mersenne info
const char *dat_path = "standalone/mersennetwister/data/MersenneTwister.dat";
static unsigned int SEED = 191871;

//MTRand* mt = NULL;

///////////////////////////////////////////////////////////////////////////////
// initialization/deinitialization code
///////////////////////////////////////////////////////////////////////////////

void cleanup() {
  if (initialized == 1) {
    printf("Deinitializing firmsim_cu.\n");

    // Free device memory
    cudaFree(d_fnprod);
    cudaFree(d_fmind);
    cudaFree(d_fpinds);
    cudaFree(d_fage);
    cudaFree(d_fqpow1);
    cudaFree(d_fquals);
    cudaFree(d_fexited);

    cudaFree(d_fmpos);
    cudaFree(d_fnrand);
    cudaFree(d_qdists);
    cudaFree(d_qbins);
    cudaFree(d_ttable);

    // Free host memory
    free(h_ttable);

    // Set initialized flag
    initialized = 0;
  }
}

void initialize(int nBpow_in, int R_PER_T_IN) {
  if (initialized == 0) {
    printf("Initializing firmsim_cu.\n");

    //mt = new MTRand(SEED);

    nBpow = nBpow_in;
    nB = 1<<nBpow;

    R_PER_T = R_PER_T_IN;
    R_SIM = T_PERIODS*R_PER_T;

    // Initialize MersenneTwister
    N_PER_RNG = iAlignUp(iDivUp(nF,MT_RNG_COUNT),2);
    RAND_N = MT_RNG_COUNT*N_PER_RNG;
    loadMTGPU(dat_path);
    seedMTGPU(SEED);

    // element counts
    nel_ft = nF;
    nel_fq = nN*nF;
    nel_fm = nN*nF;
    nel_em = nF;

    // Set up memory sizes
    fqsize = sizeof(double)*nel_fq;
    fmsize = sizeof(double)*nel_fm;
    ftsize = sizeof(int)*nel_ft;
    nrsize = sizeof(double)*RAND_CHUNKS*RAND_N;
    qbsize = sizeof(double)*nB;
    qdsize = sizeof(double)*N_QDISTS*nB;
    emsize = sizeof(double)*nel_em;
    ttsize = sizeof(double)*nM*(nN+1)*nS;

    // Allocate host memory
    h_ttable = (double*)malloc(ttsize);

    // Allocate device memory
    cudaMalloc((void**)&d_fnprod,ftsize);
    cudaMalloc((void**)&d_fmind,ftsize);
    cudaMalloc((void**)&d_fpinds,fmsize);
    cudaMalloc((void**)&d_fage,ftsize);
    cudaMalloc((void**)&d_fqpow1,emsize);
    cudaMalloc((void**)&d_fquals,fqsize);
    cudaMalloc((void**)&d_fexited,ftsize);

    cudaMalloc((void**)&d_fmpos,ftsize);
    cudaMalloc((void**)&d_fnrand,nrsize);
    cudaMalloc((void**)&d_qdists,qdsize);
    cudaMalloc((void**)&d_qbins,qbsize);
    cudaMalloc((void**)&d_ttable,ttsize);

    Check_CUDA_Error("Failed cudaMalloc.");

    // Set initialized flag
    mexAtExit(cleanup);
    initialized = 1;
  }
}

///////////////////////////////////////////////////////////////////////////////
// MEX code
///////////////////////////////////////////////////////////////////////////////
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  // read in arguments
  if (nrhs != 15) {
    printf("Too few input arguments: nrhs = %i.\n",nrhs);
    return;
  }

  if (nlhs != 11) {
    printf("Too few output arguments: nlhs = %i.\n",nlhs);
    return;
  }

  // Handle matlab data
  int nBpow_in = mxGetScalar(prhs[0]);
  double* m_qbins = (double*)mxGetData(prhs[1]);
  double* m_qdists = (double*)mxGetData(prhs[2]);
  double* m_xa = (double*)mxGetData(prhs[3]);
  double* m_xb = (double*)mxGetData(prhs[4]);
  double xe = mxGetScalar(prhs[5]);
  double* m_rho = (double*)mxGetData(prhs[6]);
  double* m_epr = (double*)mxGetData(prhs[7]);
  double tau = mxGetScalar(prhs[8]);
  double g = mxGetScalar(prhs[9]);
  double r = mxGetScalar(prhs[10]);
  double eps = mxGetScalar(prhs[11]);
  double kappa = mxGetScalar(prhs[12]);
  double qmin = mxGetScalar(prhs[13]);
  int* m_r_per_t = (int*)mxGetData(prhs[14]);
  int R_PER_T_IN = *m_r_per_t;

  int qbN = mxGetN(prhs[1]);
  int qbM = mxGetM(prhs[1]);
  int qbL = max(qbN,qbM);
  if (qbL != (1<<nBpow_in)) {
    printf("nBpow wrong.\n");
    return;
  }

  if ((initialized == 1) && ((nBpow_in != nBpow) || (R_PER_T_IN != R_PER_T))) cleanup();

  // Initialize
  if (initialized == 0) initialize(nBpow_in,R_PER_T_IN);

  R_BURN = R_BURN_INIT;
  R_TOT = R_BURN+R_SIM;

  int gen_init = 1;
  seedMTGPU(SEED);
  //seedMTGPU(mt->randInt());

  double delt = 1.0/R_PER_T;
  double qdec = 1.0/(1.0+delt*g);

  // Make transition table
  double* svec;
  double ssum;
  double xat;
  double xbt;
  double xet;
  double rhot;
  double eprt;
  for (int m = 0; m < nM; m++) {
    xat = m_xa[m];
    xbt = m_xb[m];
    xet = xe;
    rhot = m_rho[m];
    eprt = m_epr[m];

    for (int n = 0; n < nN+1; n++) {
      svec = h_ttable + m*(nN+1)*nS + n*nS;

      ssum = 1.0;
      svec[STATE_EXIT] = ssum;
      ssum -= delt*(kappa);
      svec[STATE_EXPAND] = ssum;
      ssum -= delt*(xet*eprt);
      svec[STATE_GAIN_FREE] = ssum;
      ssum -= delt*(n*kappa);
      svec[STATE_GAIN_BAS_2] = ssum;
      ssum -= delt*(n*xbt*rhot);
      svec[STATE_GAIN_BAS_1] = ssum;
      ssum -= delt*(n*xbt*(1.0-rhot));
      svec[STATE_GAIN_APP] = ssum;
      ssum -= delt*(n*xat+xet*(1.0-eprt));
      svec[STATE_LOSE] = ssum;
      ssum -= delt*(n*tau);
      svec[STATE_NOTHING] = ssum;
    }
  }

  if (ssum <= 0.0) printf("delt too large.");

  // powers
  double pow1 = eps-1.0;

  // output arrays
  mxArray* m_nprod = mxCreateNumericMatrix(nF,1,mxINT32_CLASS,mxREAL);
  mxArray* m_mind = mxCreateNumericMatrix(nF,1,mxINT32_CLASS,mxREAL);
  mxArray* m_age = mxCreateNumericMatrix(nF,1,mxINT32_CLASS,mxREAL);
  mxArray* m_mpos = mxCreateNumericMatrix(nF,1,mxINT32_CLASS,mxREAL);
  mxArray* m_qpow1 = mxCreateDoubleMatrix(nF,1,mxREAL);
  mxArray* m_nprod_zero = mxCreateNumericMatrix(nF,1,mxINT32_CLASS,mxREAL);
  mxArray* m_mind_zero = mxCreateNumericMatrix(nF,1,mxINT32_CLASS,mxREAL);
  mxArray* m_age_zero = mxCreateNumericMatrix(nF,1,mxINT32_CLASS,mxREAL);
  mxArray* m_mpos_zero = mxCreateNumericMatrix(nF,1,mxINT32_CLASS,mxREAL);
  mxArray* m_qpow1_zero = mxCreateDoubleMatrix(nF,1,mxREAL);
  mxArray* m_exited = mxCreateNumericMatrix(nF,1,mxINT32_CLASS,mxREAL);

  int* h_nprod = (int*)mxGetData(m_nprod);
  int* h_mind = (int*)mxGetData(m_mind);
  int* h_age = (int*)mxGetData(m_age);
  int* h_mpos = (int*)mxGetData(m_mpos);
  double* h_qpow1 = (double*)mxGetData(m_qpow1);
  int* h_nprod_zero = (int*)mxGetData(m_nprod_zero);
  int* h_mind_zero = (int*)mxGetData(m_mind_zero);
  int* h_age_zero = (int*)mxGetData(m_age_zero);
  int* h_mpos_zero = (int*)mxGetData(m_mpos_zero);
  double* h_qpow1_zero = (double*)mxGetData(m_qpow1_zero);
  int* h_exited = (int*)mxGetData(m_exited);

  // Block size info
  int block_size_q = 32;
  int block_size_f = 16;

  int n_blocks_q = iDivUp(nN,block_size_q);
  int n_blocks_f = iDivUp(nF,block_size_f);

  dim3 block_size(block_size_q,block_size_f);
  dim3 n_blocks(n_blocks_q,n_blocks_f);

  cudaMemcpy(d_qdists,m_qdists,qdsize,cudaMemcpyHostToDevice);
  cudaMemcpy(d_qbins,m_qbins,qbsize,cudaMemcpyHostToDevice);
  cudaMemcpy(d_ttable,h_ttable,ttsize,cudaMemcpyHostToDevice);

  // Firm inital sizes
  if (gen_init == 1) {
    RandomGPU<<<32,128>>>(d_fnrand,N_PER_RNG);
    Check_CUDA_Error("Failed RandomGPU.");

    sampn<<<n_blocks_f,block_size_f>>>(d_fnprod,d_fmind,d_fage,d_fquals,d_fpinds,d_qdists,d_qbins,d_fnrand,nBpow,nB);
    Check_CUDA_Error("Failed sampn.");
  }

  int primet = 0;
  for (int rep = 0; rep < R_TOT; rep++) {
    if (rep == R_BURN) {
      primet = 1;

      calc_qpow<<<n_blocks_f,block_size_f>>>(d_fnprod,d_fquals,d_fqpow1,pow1,qmin);
      Check_CUDA_Error("Failed calc_qpow.");

      calc_mpos<<<n_blocks_f,block_size_f>>>(d_fnprod,d_fpinds,d_fmpos);
      Check_CUDA_Error("Failed calc_mpos.");

      cudaMemcpy(h_nprod_zero,d_fnprod,ftsize,cudaMemcpyDeviceToHost);
      cudaMemcpy(h_mind_zero,d_fmind,ftsize,cudaMemcpyDeviceToHost);
      cudaMemcpy(h_age_zero,d_fage,ftsize,cudaMemcpyDeviceToHost);
      cudaMemcpy(h_mpos_zero,d_fmpos,ftsize,cudaMemcpyDeviceToHost);
      cudaMemcpy(h_qpow1_zero,d_fqpow1,emsize,cudaMemcpyDeviceToHost);
      cudaMemset((void*)d_fexited,0,ftsize);
    }

    RandomGPU<<<32,128>>>(d_fnrand,RAND_CHUNKS*N_PER_RNG);
    Check_CUDA_Error("Failed RandomGPU.");

    fzsim<<<n_blocks_f,block_size_f>>>(d_fnprod,d_fmind,d_fage,d_fquals,d_fpinds,d_fexited,d_qdists,d_qbins,d_ttable,d_fnrand,nBpow,nB,primet);
    Check_CUDA_Error("Failed fzsim.");

    fqsim<<<n_blocks,block_size>>>(d_fquals,qdec,qmin);
    Check_CUDA_Error("Failed fqsim.");

    if (rep == R_TOT-1) {
      calc_qpow<<<n_blocks_f,block_size_f>>>(d_fnprod,d_fquals,d_fqpow1,pow1,qmin);
      Check_CUDA_Error("Failed calc_qpow.");

      calc_mpos<<<n_blocks_f,block_size_f>>>(d_fnprod,d_fpinds,d_fmpos);
      Check_CUDA_Error("Failed calc_mpos.");

      cudaMemcpy(h_nprod,d_fnprod,ftsize,cudaMemcpyDeviceToHost);
      cudaMemcpy(h_mind,d_fmind,ftsize,cudaMemcpyDeviceToHost);
      cudaMemcpy(h_age,d_fage,ftsize,cudaMemcpyDeviceToHost);
      cudaMemcpy(h_mpos,d_fmpos,ftsize,cudaMemcpyDeviceToHost);
      cudaMemcpy(h_qpow1,d_fqpow1,emsize,cudaMemcpyDeviceToHost);
      cudaMemcpy(h_exited,d_fexited,ftsize,cudaMemcpyDeviceToHost);
    }
  }

  plhs[0] = m_nprod;
  plhs[1] = m_mind;
  plhs[2] = m_age;
  plhs[3] = m_mpos;
  plhs[4] = m_qpow1;
  plhs[5] = m_nprod_zero;
  plhs[6] = m_mind_zero;
  plhs[7] = m_age_zero;
  plhs[8] = m_mpos_zero;
  plhs[9] = m_qpow1_zero;
  plhs[10] = m_exited;
}

