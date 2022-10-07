#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <cuda.h>
#include <cublas.h>
#include <mex.h>

using namespace std;

const int nM = 10;

///////////////////////////////////////////////////////////////////////////////
// Utility functions
///////////////////////////////////////////////////////////////////////////////

/*
inline int topBlock(int nel, int bs)
{
  return (nel+(bs-nel%bs))/bs;
}
*/

int topBlock(int nel, int bs)
{
  return ((nel%bs) != 0) ? (nel/bs+1) : (nel/bs);
}

void checkCudaError(cudaError_t err)
{
  if(!err) return;

  printf("%s\n", cudaGetErrorString(err));

  exit(EXIT_FAILURE);
}

///////////////////////////////////////////////////////////////////////////////
// Device kernels
///////////////////////////////////////////////////////////////////////////////

__global__ void spmv(double* vin, double* vout, double pout, double pdown, 
                     double* pup1, double* pup2, double* mup0, double* mup1, int nN)
{
  volatile __shared__ double vsum[nM];
  volatile __shared__ double l_pup1[nM];
  volatile __shared__ double l_pup2[nM];
  volatile __shared__ double l_mup0[nM];
  volatile __shared__ double l_mup1[nM];

  int tm = threadIdx.x;
  int tn = threadIdx.y;

  int m = blockDim.x*blockIdx.x + tm;
  int n = blockDim.y*blockIdx.y + tn;

  if ((m < nM) && (n < nN)) {
    if (tn == 0) {
      l_pup1[m] = pup1[m];
      l_pup2[m] = pup2[m];
      l_mup0[m] = mup0[m];
      l_mup1[m] = mup1[m];
    }
    __syncthreads();

    double vtmp = 0.0;
    int i = m*nN+n;

    // stay
    double pstay = 1.0-pout-(n+1)*pdown;
    if (n<nN-2) pstay -= (n+1)*l_pup2[m];
    if (n<nN-1) pstay -= (n+1)*l_pup1[m]+l_mup0[m];
    if ((m<nM-1)&&(n<nN-1)) pstay -= l_mup1[m];
    vtmp += vin[i]*pstay;

    // transitions
    if (n<nN-1) vtmp += vin[i+1]*(n+2)*pdown;
    if (n>0) vtmp += vin[i-1]*(n*l_pup1[m]+l_mup0[m]);
    if (n>1) vtmp += vin[i-2]*(n-1)*l_pup2[m];
    if ((m>0)&&(n>0)) vtmp += vin[i-nN-1]*l_mup1[m-1];

    // pesky first row, just do it
    if (n == 0) {
      vsum[m] = vin[m*nN];
      if (m < 5)  { vsum[m] += vsum[m+5]; }
      if (m < 2)  { vsum[m] += vsum[m+2]; }
      if (m == 0) { vtmp += pout + (vsum[0] + vsum[1] + vsum[4])*pdown; }
    }

    vout[i] = vtmp;
  }
}

__global__ void dinit(double* ev, double norm, int nN)
{
  int m = blockDim.x*blockIdx.x + threadIdx.x;
  int n = blockDim.y*blockIdx.y + threadIdx.y;

  if ((m < nM) && (n < nN)) {
    int i = m*nN+n;
    ev[i] = norm/double((n+1)*(m+1));
  }
}

__global__ void vec_diff(double* vin, double* vout, int N)
{
  int i = blockDim.x*blockIdx.x + threadIdx.x;

  if (i < N) {
    vout[i] -= vin[i];
  }
}

__global__ void vec_max(double* vin, int N)
{
  __shared__ double sdata[512];

  unsigned int tid = threadIdx.x;
  unsigned int i = blockIdx.x*blockDim.x + threadIdx.x;

  sdata[tid] = (i < N) ? vin[i] : 0.0;
  __syncthreads();

  // do reduction in shared mem
  if (tid < 256) { sdata[tid] = fmaxf(sdata[tid], sdata[tid + 256]); } __syncthreads();
  if (tid < 128) { sdata[tid] = fmaxf(sdata[tid], sdata[tid + 128]); } __syncthreads();
  if (tid <  64) { sdata[tid] = fmaxf(sdata[tid], sdata[tid +  64]); } __syncthreads();
    
  if (tid < 32)
  {
      sdata[tid] = fmaxf(sdata[tid], sdata[tid + 32]);
      sdata[tid] = fmaxf(sdata[tid], sdata[tid + 16]);
      sdata[tid] = fmaxf(sdata[tid], sdata[tid +  8]);
      sdata[tid] = fmaxf(sdata[tid], sdata[tid +  4]);
      sdata[tid] = fmaxf(sdata[tid], sdata[tid +  2]);
      sdata[tid] = fmaxf(sdata[tid], sdata[tid +  1]);
  }

  // write result for this block to global mem 
  if (tid == 0) vin[blockIdx.x] = sdata[0];
}

///////////////////////////////////////////////////////////////////////////////
// persistent memory
///////////////////////////////////////////////////////////////////////////////

static double* d_vold = NULL;
static double* d_vnew = NULL;
static double* d_pup1 = NULL;
static double* d_pup2 = NULL;
static double* d_mup0 = NULL;
static double* d_mup1 = NULL;
static int initialized = 0;
static int nN = 0;
static double norm = 0.0;

void cleanup(void) {
  printf("Deinitializing eigsim.\n");
  if (initialized == 1) {
    cudaFree(d_vold);
    cudaFree(d_vnew);
    cudaFree(d_pup1);
    cudaFree(d_pup2);
    cudaFree(d_mup0);
    cudaFree(d_mup1);
    
    initialized = 0;
  }
}

///////////////////////////////////////////////////////////////////////////////
// main code
///////////////////////////////////////////////////////////////////////////////
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int n,m;

  // read in params
  if (nrhs != 7) {
    printf("Too few arguments.\n");
    return;
  }

  double* h_pup1 = (double*)mxGetData(prhs[0]);
  double* h_pup2 = (double*)mxGetData(prhs[1]);
  double* h_mup0 = (double*)mxGetData(prhs[2]);
  double* h_mup1 = (double*)mxGetData(prhs[3]);
  double h_pdown = mxGetScalar(prhs[4]);
  double h_pout = mxGetScalar(prhs[5]);
  int nNin = mxGetScalar(prhs[6]);

  if ((initialized == 1) && (nNin != nN)) cleanup();
  nN = nNin;

  int nTM = nM*nN;
  size_t msize = sizeof(double)*nM;
  size_t evsize = sizeof(double)*nTM;

  // Create matlab output array
  mxArray* m_ev = mxCreateDoubleMatrix(1,nTM,mxREAL);
  double* h_ev =  (double*)mxGetData(m_ev);
  plhs[0] = m_ev;

  cudaError_t err;
  double* d_vtmp;

  dim3 threadsPerBlock(16,16);
  dim3 numBlocks(topBlock(nM,threadsPerBlock.x),topBlock(nN,threadsPerBlock.y));

  dim3 tpb_diff(64,1);
  dim3 nbl_diff(topBlock(nTM,tpb_diff.x),1);

  dim3 tpb_max(512,1);
  dim3 nbl_max(topBlock(nTM,tpb_max.x),1);
  int nred = nbl_max.x;

  cublasStatus stat;

  if (initialized == 0) {
    printf("Initializing eigsim.\n");

    // Initialize CUBLAS
    cublasInit();
    stat = cublasGetError();
    if (stat != CUBLAS_STATUS_SUCCESS) printf("Error in cublasInit! code = %i\n",stat);

    double sumn = 0.0;
    for (n = 0; n < nN; n++) {
      sumn += 1.0/double(n+1);
    }
    double summ = 0.0;
    for (m = 0; m < nM; m++) {
      summ += 1.0/double(m+1);
    }
    sumn = 1.0/sumn;
    summ = 1.0/summ;
    norm = sumn*summ;

    err = cudaMalloc((void**)&d_vold,evsize);
    checkCudaError(err);
    err = cudaMalloc((void**)&d_vnew,evsize);
    checkCudaError(err);
    err = cudaMalloc((void**)&d_pup1,msize);
    checkCudaError(err);
    err = cudaMalloc((void**)&d_pup2,msize);
    checkCudaError(err);
    err = cudaMalloc((void**)&d_mup0,msize);
    checkCudaError(err);
    err = cudaMalloc((void**)&d_mup1,msize);
    checkCudaError(err);

    dinit<<<numBlocks,threadsPerBlock>>>(d_vold,norm,nN);

    mexAtExit(cleanup);
    initialized = 1;
  } else {
    d_vtmp = d_vold;
    d_vold = d_vnew;
    d_vnew = d_vtmp;
  }

  double* h_vold = (double*)malloc(evsize);
  double* h_vnew = (double*)malloc(evsize);
  double* sum_red = (double*)malloc(nred*sizeof(double));

  err = cudaMemcpy(d_pup1,h_pup1,msize,cudaMemcpyHostToDevice);
  checkCudaError(err);
  err = cudaMemcpy(d_pup2,h_pup2,msize,cudaMemcpyHostToDevice);
  checkCudaError(err);
  err = cudaMemcpy(d_mup0,h_mup0,msize,cudaMemcpyHostToDevice);
  checkCudaError(err);
  err = cudaMemcpy(d_mup1,h_mup1,msize,cudaMemcpyHostToDevice);
  checkCudaError(err);

  int nchecks = 0;
  int maxiter = nN*2000+1;
  double diff;
  int t;
  for (t = 1; t < maxiter; t++) {
    spmv<<<numBlocks,threadsPerBlock>>>(d_vold,d_vnew,h_pout,h_pdown,d_pup1,d_pup2,d_mup0,d_mup1,nN);

    if (t%5000 == 0) {
      vec_diff<<<nbl_diff,tpb_diff>>>(d_vnew,d_vold,nTM);
      vec_max<<<nbl_max,tpb_max>>>(d_vold,nTM);
      cudaMemcpy(sum_red,d_vold,nred*sizeof(double),cudaMemcpyDeviceToHost);

      diff = 0.0;
      for (int i = 0; i < nred; i++) {
        diff = max(diff,sum_red[i]);
      }

      nchecks++;
      if (abs(diff) < 1e-12) {
        break;
      }
    }

    d_vtmp = d_vold;
    d_vold = d_vnew;
    d_vnew = d_vtmp;
  }

  //printf("t = %i, nchecks = %i\n",t,nchecks);
  if (t == maxiter) {
    printf("Hit max iterations! diff = %15.12f\n",diff);
  }

  err = cudaMemcpy(h_ev,d_vnew,evsize,cudaMemcpyDeviceToHost);
  checkCudaError(err);

  /*
  for (m = 0; m < nM; m++) {
    for (n = 0; n < nN; n++) {
      i = m*nN+n;
      sum1 += (n+1)*h_ev[i];
    }
  }
  */

  /*
  double sum1 = 0.0;
  for (int i = 0; i < nTM; i++) {
    sum1 += h_ev[i];
  }
  printf("sum1 = %12.10f\n\n",sum1);
  */

  /*
  for (int tmi = 0; tmi < 10; tmi++) {
    printf("%12.10f\n",h_ev[tmi]);
  }
  printf("\n");
  */
}

