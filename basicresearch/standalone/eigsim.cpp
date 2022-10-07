#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <mex.h>
#include <algorithm>

using namespace std;

const int nM = 10;

///////////////////////////////////////////////////////////////////////////////
// Device kernels
///////////////////////////////////////////////////////////////////////////////

void spmv(double* vin, double* vout, double pout, double pdown,
          double* pup0, double* mup0, double* mup1,
          int nN, int nK)
{
  int k;

  double vfirst = 0.0;
  for (int im = 0; im < nM; im++) {
    vfirst += vin[im*nN];
  }

  for (int im = 0; im < nM; im++) {
    double* pup = pup0 + im*nK;
    for (int in = 0; in < nN; in++) {
      int i = im*nN + in;

      int m = im + 1;
      int n = in + 1;

      double vtmp = 0.0;

      // stay
      double pstay = 1.0-pout-n*pdown;
      for (int ik = 0; ik+1 <= min(nK, nN-n); ik++) {
        pstay -= n*pup[ik];
      }
      if (n<nN) pstay -= mup0[im];
      if ((m<nM)&&(n<nN)) pstay -= mup1[im];
      vtmp += vin[i]*pstay;

      // products
      if (n<nN) vtmp += vin[i+1]*(n+1)*pdown;
      if (n>1) vtmp += vin[i-1]*mup0[im];
      for (int ik = 0; ik+1 <= min(nK, n-1); ik++) {
        k = ik + 1;
        vtmp += vin[i-k]*(n-k)*pup[ik];
      }

      // spillover
      if ((m>1)&&(n>1)) vtmp += vin[i-nN-1]*mup1[im-1];
      // pesky first row
      if ((n==1)&&(m==1)) vtmp += pout + vfirst*pdown;

      vout[i] = vtmp;
    }
  }
}

void dinit(double* ev, double norm, int nN)
{
  for (int m = 0; m < nM; m++) {
    for (int n = 0; n < nN; n++) {
      int i = m*nN+n;
      ev[i] = norm/double((n+1)*(m+1));
    }
  }
}

void vec_diff(double* vin, double* vout, int N)
{
  for (int i = 0; i < N; i++) {
    vout[i] -= vin[i];
  }
}

double vec_max(double* vin, int N)
{
  double ret = 0.0;
  for (int i = 0; i < N; i++) {
    ret = fmax(ret,vin[i]);
  }
  return ret;
}

double vec_sum(double* vin, int N)
{
  double ret = 0.0;
  for (int i = 0; i < N; i++) {
    ret += vin[i];
  }
  return ret;
}

///////////////////////////////////////////////////////////////////////////////
// persistent memory
///////////////////////////////////////////////////////////////////////////////

static double* h_vold = NULL;
static double* h_vnew = NULL;
static int initialized = 0;
static int nN = 0;
static int nK = 0;
static double norm = 0.0;

void cleanup(void) {
  printf("Deinitializing eigsim.\n");
  free(h_vold);
  free(h_vnew);
  initialized = 0;
}

///////////////////////////////////////////////////////////////////////////////
// main code
///////////////////////////////////////////////////////////////////////////////
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int n,m;

  // read in params
  if (nrhs != 8) {
    printf("MEX: Incorrect number of arguments.\n");
    return;
  }

  double* h_pup = (double*)mxGetData(prhs[0]);
  double* h_mup0 = (double*)mxGetData(prhs[1]);
  double* h_mup1 = (double*)mxGetData(prhs[2]);
  double h_pdown = mxGetScalar(prhs[3]);
  double h_pout = mxGetScalar(prhs[4]);
  int nNin = mxGetScalar(prhs[5]);
  int nKin = mxGetScalar(prhs[6]);
  int reset = mxGetScalar(prhs[7]);

  if ((initialized == 1) && (nNin != nN)) cleanup();
  nN = nNin;
  nK = nKin;

  int nTM = nM*nN;
  size_t msize = sizeof(double)*nM;
  size_t evsize = sizeof(double)*nTM;

  // Create matlab output array
  mxArray* m_ev = mxCreateDoubleMatrix(1,nTM,mxREAL);
  double* h_ev =  (double*)mxGetData(m_ev);
  plhs[0] = m_ev;

  double* h_vtmp;

  if (reset == 1) {
      cleanup();
  }

  if (initialized == 0) {
    printf("Initializing eigsim.\n");

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

    h_vold = (double*)malloc(evsize);
    h_vnew = (double*)malloc(evsize);

    dinit(h_vold,norm,nN);

    initialized = 1;
  } else {
    h_vtmp = h_vold;
    h_vold = h_vnew;
    h_vnew = h_vtmp;
  }

  int nchecks = 0;
  int maxiter = nN*2000;
  double diff;
  double vsum;

  // vsum = vec_sum(h_vold,nTM);
  // printf("%d: sum = %15.12f\n",0,vsum);

  int t;
  for (t = 1; t-1 < maxiter; t++) {
    spmv(h_vold,h_vnew,h_pout,h_pdown,h_pup,h_mup0,h_mup1,nN,nK);

    if (t%50 == 0) {
      vec_diff(h_vnew,h_vold,nTM);
      diff = vec_max(h_vold,nTM);
      vsum = vec_sum(h_vnew,nTM);

      //printf("t = %d, diff = %12.10f, sum = %12.10f\n", t, diff, sum1);

      nchecks++;
      if (abs(diff) < 1e-12) {
        break;
      }
    }

    h_vtmp = h_vold;
    h_vold = h_vnew;
    h_vnew = h_vtmp;
  }

  // printf("t = %i, nchecks = %i\n",t,nchecks);
  if (t == maxiter) {
    printf("Hit max iterations! diff = %15.12f\n",diff);
  }

  memcpy(h_ev,h_vnew,evsize);

  /*
  double sum1 = vec_sum(h_vnew,nTM);
  printf("sum1 = %12.10f\n\n",sum1);
  */

  /*
  for (int tmi = 0; tmi < 10; tmi++) {
    printf("%12.10f\n",h_ev[tmi]);
  }
  printf("\n");
  */
}
