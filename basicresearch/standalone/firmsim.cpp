#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <mex.h>

#include "MersenneTwister.h"

using namespace std;

///////////////////////////////////////////////////////////////////////////////
// Constants
///////////////////////////////////////////////////////////////////////////////
//#define nF 1024
//#define nF 4096
//#define nF 8192
//#define nF 16384
#define nF 32768

#define R_BURN 10000

#define T_PERIODS 1

#define nN 128
#define nNpow 6
#define nM 10
#define nK 11

#define nS 8
#define nSpow 3
#define STATE_NOTHING1     0
#define STATE_NOTHING2     1
#define STATE_LOSE         2
#define STATE_GAIN_APP     3
#define STATE_GAIN_BAS     4
#define STATE_GAIN_FREE    5
#define STATE_EXPAND       6
#define STATE_EXIT         7

///////////////////////////////////////////////////////////////////////////////
// globals
///////////////////////////////////////////////////////////////////////////////

static int R_SIM = 0;
int R_PER_T = 0;

///////////////////////////////////////////////////////////////////////////////
// random chunker
///////////////////////////////////////////////////////////////////////////////

#define CHUNK_SIZE 1024

MTRand* mt = NULL;
MTRand* mt_s = NULL;
MTRand* mt_f = NULL;
MTRand* mt_nf = NULL;
int num_left = 0;
static double* rand_vec = NULL;

void init_chunker()
{
  rand_vec = (double*)malloc(sizeof(double)*CHUNK_SIZE);
}

void deinit_chunker()
{
  free(rand_vec);
}

void seed_chunker(int seed_in)
{
  mt->seed(seed_in);
  num_left = 0;
}

void fill_random(double* rout, int nR, MTRand* mt_in)
{
  for (int i = 0; i < nR; i++) {
    rout[i] = mt_in->rand();
  }
}

inline double get_rand()
{
  //return mt->rand();

  if (num_left == 0) {
    fill_random(rand_vec,CHUNK_SIZE,mt);
    num_left = CHUNK_SIZE;
  }

  num_left--;

  return rand_vec[num_left];
}

///////////////////////////////////////////////////////////////////////////////
// Device kernels
///////////////////////////////////////////////////////////////////////////////

// sequential search for small len
inline int dsamp1(double* pbeg, int len, double r)
{
  double* ppos = pbeg;
  double* pend = pbeg + len;
  int s = 0;
  while (ppos < pend) {
    if (r < (*ppos)) {
      break;
    }
    s++;
    ppos++;
  }
  return s;
}

// binary search (only powers of 2 sizes)
inline int dsamp2(double* pbeg, int pow2, double r)
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

void fsim_single(double* tvecs, double* qdists, double* qbins, double* pupa, double* pupb, double* srand, double qdec, double qmin, double eps, int* nprod_out, int* mind_out, int* fage_out, int* mpos_out, double* empl_out, int* nprod_zero_out, int* mind_zero_out, int* fage_zero_out, int* mpos_zero_out, double* empl_zero_out, int* exited_out, int nBpow, int nB, int disp, int f)
{
  int s;
  int qind;
  int mind;
  int k;
  int i;
  double qval;

  int r_sub = 0;
  double emp;
  double prod;

  double empl_zero;
  int mind_zero;

  double* tvec;

  int age = 0;
  int exited = -1;
  int m = 0;
  int n = 0;
  double quals[nN];
  int pinds[nN];

  for (int rep = 0; rep < R_BURN+R_SIM+1; rep++) {
    // handle exits and initialization
    if (n == 0) {
      seed_chunker(mt_nf->randInt());

      age = 0;
      exited = (exited == -1) ? 0 : 1;
      m = 1;
      n = 1;

      qind = dsamp2(qdists,nBpow,get_rand());
      mind = floor(get_rand()*m);
      qval = qbins[qind];
      quals[0] = qval;
      pinds[0] = mind;
    }

    // aggregate statistics
    if (rep >= R_BURN) {
      if (rep == R_BURN) {
        r_sub = 0;

        emp = 0.0;
        for (qind = 0; qind < n; qind++) {
          qval = quals[qind];
          emp += pow(qval,eps-1.0);
        }

        int mpos = 0;
        int mpres[nM] = {0,0,0,0,0,0,0,0,0,0};
        for (qind = 0; qind < n; qind++) {
          mind = pinds[qind];
          mpres[mind] = 1;
        }
        for (mind = 0; mind < nM; mind++) {
          if (mpres[mind] == 1) mpos++;
        }

        empl_zero_out[f] = emp;
        mind_zero_out[f] = m;
        nprod_zero_out[f] = n;
        fage_zero_out[f] = age;
        mpos_zero_out[f] = mpos;

        exited = 0;
      }

      if (r_sub == R_PER_T) {
        emp = 0.0;
        for (qind = 0; qind < n; qind++) {
          qval = quals[qind];
          emp += pow(qval,eps-1.0);
        }

        int mpos = 0;
        int mpres[nM] = {0,0,0,0,0,0,0,0,0,0};
        for (qind = 0; qind < n; qind++) {
          mind = pinds[qind];
          mpres[mind] = 1;
        }
        for (mind = 0; mind < nM; mind++) {
          if (mpres[mind] == 1) mpos++;
        }

        empl_out[f] = emp;
        mind_out[f] = m;
        nprod_out[f] = n;
        exited_out[f] = exited;
        fage_out[f] = age;
        mpos_out[f] = mpos;

        //printf("f = %i: n = %i, emp = %f, prod = %f\n",f,n,emp,prod);
        //if (f == 281) printf("r_sub = %i, rep = %i\n",r_sub,rep);

        return;
      }

      r_sub++;
    }

    // the sampler
    tvec = tvecs + (m-1)*(nN+1)*nS + n*nS;
    s = dsamp1(tvec,nS,srand[rep]);

    //if ((f == 281)&&(s>0)) {
    //  printf("rep = %i, s = %i: m = %i, n = %i, r = %15.12f\n",rep,s,m,n,srand[rep]);
    //}

    switch (s) {
      case STATE_GAIN_APP:
        k = dsamp1(pupa+(m-1)*nK,nK,get_rand()) + 1;
        for (i = 0; i < k; i++) {
          if (n < nN) {
            qind = dsamp2(qdists,nBpow,get_rand());
            mind = floor(get_rand()*m);
            quals[n] = qbins[qind];
            pinds[n] = mind;
            n++;
          }
        }
        break;
      case STATE_GAIN_BAS:
        k = dsamp1(pupb+(m-1)*nK,nK,get_rand()) + 1;
        for (i = 0; i < k; i++) {
          if (n < nN) {
            qind = dsamp2(qdists+nB,nBpow,get_rand());
            mind = floor(get_rand()*m);
            quals[n] = qbins[qind];
            pinds[n] = mind;
            n++;
          }
        }
        break;
        case STATE_GAIN_FREE:
          if (n < nN) {
            qind = dsamp2(qdists+2*nB,nBpow,get_rand());
            mind = floor(get_rand()*m);
            quals[n] = qbins[qind];
            pinds[n] = mind;
            n++;
          }
          break;
      case STATE_LOSE:
        if (n > 1) {
          qind = floor(get_rand()*n);
          quals[qind] = quals[n-1];
          pinds[qind] = pinds[n-1];
        }
        n--;
        break;
      case STATE_EXPAND:
        if (m < nM) {
          qind = dsamp2(qdists,nBpow,get_rand());

          quals[n] = qbins[qind];
          pinds[n] = m;

          m++;
          n++;
        }
        break;
      case STATE_EXIT:
        n = 0;
        break;
    }

    // decrement by growth
    for (qind = 0; qind < n; qind++) {
      qval = qdec*quals[qind];
      quals[qind] = (qval > qmin) ? qval : qmin;
    }

    // increment age
    age++;
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

static int R_TOT = 0;

static int nel_tt = 0;
static int nel_qd = 0;
static int nel_qb = 0;

static size_t fqsize = 0;
static size_t ftsize = 0;
static size_t emsize = 0;
static size_t ttsize = 0;
static size_t qdsize = 0;
static size_t qbsize = 0;
static size_t srsize = 0;

static int EMPL_N_BLOCKS = 0;
static size_t empl_red_size = 0;

// Device memory
static double* h_ttable = NULL;
static double* h_srand = NULL;

// Flags
static int initialized = 0;

// Mersenne info
static unsigned int SEED = 191871;
//static unsigned int SEED = 384092348; // orig
//static unsigned int SEED_F = 91283982;
//static unsigned int SEED_NF = 1429043024;

///////////////////////////////////////////////////////////////////////////////
// initialization/deinitialization code
///////////////////////////////////////////////////////////////////////////////

void cleanup() {
  if (initialized == 1) {
    printf("Deinitializing firmsim_cpp.\n");

    deinit_chunker();

    // Free device memory
    free(h_ttable);
    free(h_srand);

    // Set initialized flag
    initialized = 0;
  }
}

void initialize(int nBpow_in, int R_PER_T_IN) {
  if (initialized == 0) {
    printf("Initializing firmsim_cpp.\n");

    init_chunker();

    nBpow = nBpow_in;
    nB = 1<<nBpow;

    R_PER_T = R_PER_T_IN;
    R_SIM = R_PER_T;
    R_TOT = R_BURN+R_SIM;

    // element counts
    nel_tt = nM*(nN+1)*nS;
    nel_qd = 3*nB;
    nel_qb = nB;

    // Set up memory sizes
    fqsize = sizeof(double)*nN*nF;
    ftsize = sizeof(int)*nF;
    emsize = sizeof(double)*nF;
    ttsize = sizeof(double)*nel_tt;
    qdsize = sizeof(double)*nel_qd;
    qbsize = sizeof(double)*nel_qb;
    srsize = sizeof(double)*R_TOT;

    // Allocate memory
    h_ttable = (double*)malloc(ttsize);
    h_srand = (double*)malloc(srsize);

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
  if (nrhs != 16) {
    printf("Too few input arguments.\n");
    return;
  }

  if (nlhs != 11) {
    printf("Too few output arguments.\n");
    return;
  }

  // Handle matlab data
  int nBpow_in = mxGetScalar(prhs[0]);
  double* m_qbins = (double*)mxGetData(prhs[1]);
  double* m_qdists = (double*)mxGetData(prhs[2]);
  double* m_xa = (double*)mxGetData(prhs[3]);
  double* m_xb = (double*)mxGetData(prhs[4]);
  double xe = mxGetScalar(prhs[5]);
  double* m_epr = (double*)mxGetData(prhs[6]);
  double* m_pupa = (double*)mxGetData(prhs[7]);
  double* m_pupb = (double*)mxGetData(prhs[8]);
  double tau = mxGetScalar(prhs[9]);
  double g = mxGetScalar(prhs[10]);
  double r = mxGetScalar(prhs[11]);
  double eps = mxGetScalar(prhs[12]);
  double kappa = mxGetScalar(prhs[13]);
  double qmin = mxGetScalar(prhs[14]);
  int* m_r_per_t = (int*)mxGetData(prhs[15]);
  int R_PER_T_IN = *m_r_per_t;

  int qbN = mxGetN(prhs[1]);
  int qbM = mxGetM(prhs[1]);
  int qbL = max(qbN,qbM);
  if (qbL != (1<<nBpow_in)) {
    printf("nBpow wrong.\n");
    return;
  }

  if ((initialized == 1) && ((nBpow_in != nBpow) || (R_PER_T_IN != R_PER_T))) cleanup();
  if (initialized == 0) initialize(nBpow_in,R_PER_T_IN);

  // reinit random every time
  if (mt) delete mt;
  if (mt_s) delete mt_s;
  if (mt_f) delete mt_f;
  if (mt_nf) delete mt_nf;

  mt = new MTRand(SEED);
  mt_s = new MTRand(mt->randInt());
  mt_f = new MTRand(mt->randInt());
  mt_nf = new MTRand(mt->randInt());
  num_left = 0;

  // input values
  double delt = 1.0/R_PER_T;
  double qdec = 1.0/(1.0+delt*g);

  // Make transition table
  double* svec;
  double ssum;
  double xat;
  double xbt;
  double xet;
  double eprt;
  for (int m = 0; m < nM; m++) {
    xat = m_xa[m];
    xbt = m_xb[m];
    xet = xe;
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
      svec[STATE_GAIN_BAS] = ssum;
      ssum -= delt*(n*xbt);
      svec[STATE_GAIN_APP] = ssum;
      ssum -= delt*(n*xat+xet*(1.0-eprt));
      svec[STATE_LOSE] = ssum;
      ssum -= delt*(n*tau);
      svec[STATE_NOTHING2] = ssum;
      svec[STATE_NOTHING1] = ssum;
    }

    if (ssum <= 0.0) printf("delt too large.");
  }

  /*
  double qmean = 0.0;
  int qrep = 10000;
  int qind;
  for (int i = 0; i < qrep; i++) {
    qind = dsamp2(m_qdists+2*nB,nBpow,mt->rand());
    qmean += m_qbins[qind];
  }
  qmean /= (double)qrep;
  printf("qmean = %f\n",qmean);
  */

  // output arrays
  mxArray* m_nprod = mxCreateNumericMatrix(nF,1,mxINT32_CLASS,mxREAL);
  mxArray* m_mind = mxCreateNumericMatrix(nF,1,mxINT32_CLASS,mxREAL);
  mxArray* m_age = mxCreateNumericMatrix(nF,1,mxINT32_CLASS,mxREAL);
  mxArray* m_mpos = mxCreateNumericMatrix(nF,1,mxINT32_CLASS,mxREAL);
  mxArray* m_empl = mxCreateDoubleMatrix(nF,1,mxREAL);
  mxArray* m_nprod_zero = mxCreateNumericMatrix(nF,1,mxINT32_CLASS,mxREAL);
  mxArray* m_mind_zero = mxCreateNumericMatrix(nF,1,mxINT32_CLASS,mxREAL);
  mxArray* m_age_zero = mxCreateNumericMatrix(nF,1,mxINT32_CLASS,mxREAL);
  mxArray* m_mpos_zero = mxCreateNumericMatrix(nF,1,mxINT32_CLASS,mxREAL);
  mxArray* m_empl_zero = mxCreateDoubleMatrix(nF,1,mxREAL);
  mxArray* m_exited = mxCreateNumericMatrix(nF,1,mxINT32_CLASS,mxREAL);

  int* h_nprod = (int*)mxGetData(m_nprod);
  int* h_mind = (int*)mxGetData(m_mind);
  int* h_age = (int*)mxGetData(m_age);
  int* h_mpos = (int*)mxGetData(m_mpos);
  double* h_empl = (double*)mxGetData(m_empl);
  int* h_nprod_zero = (int*)mxGetData(m_nprod_zero);
  int* h_mind_zero = (int*)mxGetData(m_mind_zero);
  int* h_age_zero = (int*)mxGetData(m_age_zero);
  int* h_mpos_zero = (int*)mxGetData(m_mpos_zero);
  double* h_empl_zero = (double*)mxGetData(m_empl_zero);
  int* h_exited = (int*)mxGetData(m_exited);

  for (int f = 0; f < nF; f++) {
    fill_random(h_srand,R_TOT+1,mt_s);

    seed_chunker(mt_f->randInt());

    fsim_single(h_ttable,m_qdists,m_qbins,m_pupa,m_pupb,h_srand,qdec,qmin,eps,h_nprod,h_mind,h_age,h_mpos,h_empl,h_nprod_zero,h_mind_zero,h_age_zero,h_mpos_zero,h_empl_zero,h_exited,nBpow,nB,0,f);
  }

  plhs[0] = m_nprod;
  plhs[1] = m_mind;
  plhs[2] = m_age;
  plhs[3] = m_mpos;
  plhs[4] = m_empl;
  plhs[5] = m_nprod_zero;
  plhs[6] = m_mind_zero;
  plhs[7] = m_age_zero;
  plhs[8] = m_mpos_zero;
  plhs[9] = m_empl_zero;
  plhs[10] = m_exited;
}
