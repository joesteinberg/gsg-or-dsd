#ifndef __CALIBRATE_H__
#define __CALIBRATE_H__

#include "globals.h"

typedef struct
{
  // ......................................................................
  // final demand
  // ......................................................................

  //combine final demand from different sectors... cons is CES, inv is Cobb-Douglas
  // Cons = [eps * goods^rho + (1-eps) * services^rho]^(1/rho)
  // Inv = G * [goods^eps * services^(1-eps)]
  double rho;
  double eps[NC][NF][NS];
  double G[NC];

  // combine final demand from different countries into sector-specific bundles
  // f_s = H * [sum_{j=1}^{NC} (theta_j * f_j)^sig]^(1/sig)
  double sig[NC][NS-1];
  double theta[NC][NS-1][NC];
  double H[NC][NS-1];

  // ......................................................................
  // gross output parameters
  // ......................................................................

  // combine intermediates from different countries into sector-specific bundles
  // M_s = C * [sum_{j=1}^{NC} (mu_j * m_j)^zeta]^(1/zeta)
  double zeta[NC][NS-1];
  double mu[NC][NS-1][NC];
  double M[NC][NS-1];

  // combine value added and intermediate bundles from different sectors... Leontief
  // Gross output = min[VA/lam_va, M_goods/lam_goods, M_svcs/lam_svcs]
  double lam_va[NC][NS];
  double lam[NC][NS][NS-1];

  // value added... Cobb-Douglas
  // VA = B * [k^alpha * (gam * ell)^(1-alpha)]
  double alpha[NC][NS];
  double A[NC][NS];
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
  double y0[NC][NS];
  double va0[NC][NS];
  double k0[NC][NS];
  double l0[NC][NS];
  double md0[NC][NS][NS-1];
  double ex0[NC][NC];
  double im0[NC][NC];
  double nx0[NC][NC];
  double c0[NC][NS-1];
  double i0[NC][NS];
  double m0[NC][NS-1];
  double m02[NC][NS-1][NC];
  double q0[NC][NS-1];
  double q02[NC][NS-1][NC];
  double im02[NC][NS-1][NC];
  double ex02[NC][NS-1][NC];
  double nx02[NC][NS-1][NC];
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

static inline double prod_go(double va, const double md[], double lam_va, const double lam[])
{
  double mm = HUGE_VAL;
  uint i;
  for(i=0; i<NS-1; i++)
    {
      mm = fmin(mm,md[i]/lam[i]);
    }
  return fmin( va/lam_va, mm );
}

static inline double prod_va(double k, double l, double A, double alpha, double a)
{
  return A * pow(k,alpha) * pow(a*l,(1.0 - alpha));
}

static inline double prod_inv(const double x[], const double eps[], double G)
{
  return G * pow(x[0],eps[0]) * pow(x[1],eps[1]) * pow(x[2],eps[2]);
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

static inline double muc(const double c[], double l, double lbar, double ne, double nw, const double eps[], double rho, double phi, double psi, uint s)
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

  return phi * eps[s] * pow(c[s]/ne,rho-1.0) * 
    pow( eps[0]*pow(c[0]/ne,rho) + eps[1]*pow(c[1]/ne,rho), psi*phi/rho-1.0 ) * 
    pow(leisure/nw,(1.0-phi)*psi) / ne;
}

static inline double mul(const double c[], double l, double lbar, double ne, double nw, const double eps[], double rho, double phi, double psi)
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

  return (1.0-phi) * 
    pow( eps[0]*pow(c[0]/ne,rho) + eps[1]*pow(c[1]/ne,rho), psi*phi/rho )  * 
    pow(leisure/nw,(1.0-phi)*psi - 1.0) / nw;
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
