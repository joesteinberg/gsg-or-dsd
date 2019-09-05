#ifndef __CALIBRATE_H__
#define __CALIBRATE_H__

#include "globals_1sector.h"

typedef struct
{
  // ......................................................................
  // final demand
  // ......................................................................

  // combine final demand from different countries into sector-specific bundles
  // f = H * [sum_{j=1}^{NC} (theta_j * f_j)^sig]^(1/sig)
  double sig[NC];
  double theta[NC][NC];
  double H[NC];

  // ......................................................................
  // gross output parameters
  // ......................................................................

  // combine intermediates from different countries into sector-specific bundles
  // M = C * [sum_{j=1}^{NC} (mu_j * m_j)^zeta]^(1/zeta)
  double zeta[NC];
  double mu[NC][NC];
  double M[NC];

  // combine value added and intermediate bundles from different sectors... Leontief
  // Gross output = min[VA/lam_va, M/lam]
  double lam_va[NC];
  double lam[NC];

  // value added... Cobb-Douglas
  // VA = A * [k^alpha * (gam * ell)^(1-alpha)]
  double alpha[NC];
  double A[NC];
  double ga[NC];
  double ga_bgp;
  double conv_speed;

  // ......................................................................
  // households
  // ......................................................................
  
  // capital formation
  double delta; // depreciation rate
  double tauk[NC];  // capital tax rate
  double rss; // steady-state real interest rate
  
  // household preferences
  double beta[NC]; // discount factors
  double psi; // intertemporal elasticity
  double phi[NC]; // consumption share

  // endowments
  double kk0[NC];
  double b0[NC];
  double gpopt[NC];
  double gpopw[NC];

  // ......................................................................
  // time series parameters
  // ......................................................................

  // exogenous params
  double popt_ts[NT+1][NC];
  double popw_ts[NT+1][NC];
  double pope_ts[NT+1][NC];
  double lbar_ts[NT+1][NC];
  double a_ts[NT+1][NC];
  double ga_ts[NT+1][NC];

  // wedges
  double wedge_speed;

  // data to match
  double RGDP_ts[NMATCH][NC];
  double C_frac_ts[NMATCH][NC];
  double TB_frac_ts[NMATCH][NC];
  double I_frac_ts[NMATCH][NC];
  double RER_ts[NMATCH+2];
  double RIR_ts[NMATCH];

  // ......................................................................
  // base-period equilibrium values
  // ......................................................................
  double iomat[NS*NC+2][NS*NC + NF*NC + 1];
  double r0[NC];
  double ii0[NC];
  double ll0[NC];
  double y0[NC];
  double va0[NC];
  double k0[NC];
  double l0[NC];
  double md0[NC];
  double ex0[NC][NC];
  double im0[NC][NC];
  double nx0[NC][NC];
  double c0[NC];
  double i0[NC];
  double m0[NC];
  double m02[NC][NC];
  double q0[NC];
  double q02[NC][NC];
  double im02[NC][NC];
  double ex02[NC][NC];
  double nx02[NC][NC];
  double lbar0[NC];

  // ......................................................................
  // adjustment cost parameters
  // ......................................................................
  double etaM;
  double etaF;
  double etaK;
 
}params;

params ppp0[NTH];
params ppp1[NTH];
params ppp2[NTH];

void set_nontargeted_params(params * p);
uint load_iomat(params * p);
uint load_ts_params(params * p);
uint store_base_period_values(params * p);
uint calibrate_prod_params(params * p);
uint calibrate_fin_params(params * p);
uint calibrate_hh_params(params * p);
uint calibrate(params * p);
uint write_params(const params * p);

static inline double prod_go(double va, const double md, double lam_va, const double lam)
{
  return fmin( va/lam_va, md/lam );
}

static inline double prod_va(double k, double l, double A, double alpha, double a)
{
  return A * pow(k,alpha) * pow(a*l,(1.0 - alpha));
}

static inline double prod_m(const double m2[], double M, const double mu[], double zeta)
{
  double tmp = 0.0;
  uint i;
  for(i=0; i<NC; i++)
    {
      tmp = tmp + mu[i] * pow(m2[i],zeta);
    }
  return M * pow(tmp, 1.0/zeta );
}

static inline double prod_q(const double q2[], double H, const double theta[], double sig)
{
  double tmp = 0.0;
  uint i;
  for(i=0; i<NC; i++)
    {
      tmp = tmp + theta[i] * pow(q2[i],sig);
    }
  return H * pow(tmp, 1.0/sig );
}

static inline double muc(const double c, double l, double lbar, double ne, double nw, double phi, double psi)
{
  double leisure;
  if(lbar-l > 0.0001)
    {
      leisure = lbar-l;
    }
  else
    {
      leisure = 0.0001 / log(0.0001-(lbar-l));
    }

  return phi * pow(c/ne, psi*phi-1.0 ) * pow(leisure/nw,(1.0-phi)*psi) / ne;
}

static inline double mul(const double c, double l, double lbar, double ne, double nw, double phi, double psi)
{
  double leisure;
  if(lbar-l > 0.0001)
    {
      leisure = lbar-l;
    }
  else
    {
      leisure = 0.0001 / log(0.0001-(lbar-l));
    }

  return (1.0-phi) * pow(c/ne, psi*phi ) * pow(leisure/nw,(1.0-phi)*psi - 1.0) / nw;
}

static inline double phiK(double x, double delta, double ga, double etaK)
{
  return (pow(delta+ga-1.0,1.0-etaK) * pow(x,etaK) - (1.0-etaK)*(delta+ga-1.0))/etaK;
}

static inline double dphiK(double x, double delta, double ga, double etaK)
{
  return pow(delta+ga-1.0,1.0-etaK) * pow(x,etaK-1.0);
}

#endif
