#ifndef __CALIBRATE_C__
#define __CALIBRATE_C__

#include "calibrate_1sector.h"

void set_nontargeted_params(params * p)
{
  uint i;

  // parameters common across countries
  p->rss = 0.025;
  p->ga_bgp = 1.02;
  p->conv_speed = 0.75;
  p->delta = 0.06;

  p->etaK = 0.8;

  p->etaM = 0.00;
  p->etaF = 0.00;
  p->psi = -1.0;

  p->wedge_speed = 0.65;

  SET_ALL_V(p->alpha,NC,0.34);

  //SET_ALL_V(p->r0,NC,0.025);
  //SET_ALL_V(p->tauk,NC,0.25);
  p->tauk[0] = 0.39; // statutory rate for USA
  p->tauk[1] = 0.34; // average statutory rate for non-US countries

  if(sym_g_flag==1)
    {
      for(i=0; i<NC; i++)
	{
	  p->ga[i] = p->ga_bgp;
	}
    }
  else
    {
      p->ga[0] = 1.0201075;
      p->ga[1] = 1.0309967;
    }

  for(i=0; i<NC; i++)
    {
      if(noio_flag==1)
	{
	  p->sig[i] = 1.0-1.0/2.0;
	}
      else
	{
	  p->sig[i] = 1.0-1.0/1.5;
	}

      p->zeta[i] = 1.0-1.0/2.0;
    }

  return;
}

uint load_iomat(params * p)
{
  uint i, j, got;
  double tmp;
  FILE * file;
  if(noio_flag==1)
    {
      file = fopen("../python/output/iomat_noio_2c.txt","rb");
    }
  else
    {
      file = fopen("../python/output/iomat_2c.txt","rb");
    }
  if(file)
    {
      for(i=0; i<(NS*NC+2); i++)
	{
	  for(j=0; j<(NS*NC+NF*NC+1); j++)
	    {
	      got = fscanf(file,"%lf",&tmp);
	      if(got != 1)
		{
		  fprintf(logfile,KRED "Error reading IO matrix!\n" RESET);
		  fclose(file);
		  return 1;
		}
	      p->iomat[i][j] = tmp;
	    }
	}
      fclose(file);
      return 0;
    }
  else
    {
      fprintf(logfile,KRED "Error loading IO matrix!\n" RESET);
      return 1;
    }
}

uint load_ts_params(params * p)
{
  double tmp, tmp1, tmp2;
  FILE * file;
  uint t, got, cnt;

  //  USA data
  file = fopen("../python/output/rgdppc_USA.txt","rb");
  for(t=0; t<NMATCH; t++)
    {
      got = fscanf(file,"%lf",&tmp);
      if(got != 1)
	{
	  fprintf(logfile,KRED "Error reading text file!\n" RESET);
	  fclose(file);
	  return 1;
	}
      p->RGDP_ts[t][0] = tmp;
    }
  fclose(file);

  file = fopen("../python/output/C_frac.txt","rb");
  for(t=0; t<NMATCH; t++)
    {
      got = fscanf(file,"%lf",&tmp);
      if(got != 1)
	{
	  fprintf(logfile,KRED "Error reading text file!\n" RESET);
	  fclose(file);
	  return 1;
	}
      p->C_frac_ts[t][0] = tmp;
    }
  fclose(file);

  file = fopen("../python/output/TB_frac.txt","rb");
  for(t=0; t<NMATCH; t++)
    {
      got = fscanf(file,"%lf",&tmp);
      if(got != 1)
	{
	  fprintf(logfile,KRED "Error reading text file!\n" RESET);
	  fclose(file);
	  return 1;
	}
      p->TB_frac_ts[t][0] = tmp;
    }
  fclose(file);

  file = fopen("../python/output/I_frac.txt","rb");
  for(t=0; t<NMATCH; t++)
    {
      got = fscanf(file,"%lf",&tmp);
      if(got != 1)
	{
	  fprintf(logfile,KRED "Error reading text file!\n" RESET);
	  fclose(file);
	  return 1;
	}
      p->I_frac_ts[t][0] = tmp;
    }
  fclose(file);

  file = fopen("../python/output/rer.txt","rb");
  for(t=0; t<NMATCH+2; t++)
    {
      got = fscanf(file,"%lf",&tmp);
      if(got != 1)
	{
	  fprintf(logfile,KRED "Error reading text file!\n" RESET);
	  fclose(file);
	  return 1;
	}
      p->RER_ts[t] = tmp;
    }
  fclose(file);

  file = fopen("../python/output/rdata.csv","rb");
  for(t=0; t<NMATCH; t++)
    {
      got = fscanf(file,"%lf",&tmp);
      if(got != 1)
	{
	  fprintf(logfile,KRED "Error reading text file!\n" RESET);
	  fclose(file);
	  return 1;
	}
      p->RIR_ts[t] = tmp;
    }
  fclose(file);

  file = fopen("../python/output/demo_USA.txt","rb");
  cnt=0;
  for(t=0; t<NT+1; t++)
    {
      got = fscanf(file,"%lf %lf",&tmp1,&tmp2);
      if(got != 2)
	{
	  fprintf(logfile,KRED "Error reading text file!\n" RESET);
	  fclose(file);
	  return 1;
	}
      if(t<TCONV0)
	{
	  p->popt_ts[t][0] = tmp1;
	  p->popw_ts[t][0] = tmp2;
	}
      else if(t<TCONV1)
	{
	  double frac = ((double)cnt)/(TCONV1-TCONV0);
	  p->popt_ts[t][0] = p->popt_ts[t-1][0]*(1.0+frac*(tmp1/p->popt_ts[t-1][0]-1.0));
	  p->popw_ts[t][0] = p->popw_ts[t-1][0]*(1.0+frac*(tmp2/p->popw_ts[t-1][0]-1.0));
	  cnt=cnt+1;
	}
      else
	{
	  p->popt_ts[t][0] = p->popt_ts[t-1][0];
	  p->popw_ts[t][0] = p->popw_ts[t-1][0];
	}
    }
  fclose(file);

  cnt=0;
  file = fopen("../python/output/lp_USA.txt","rb");
  for(t=0; t<NT+1; t++)
    {
      if(t<NMATCH)
	{
	  got = fscanf(file,"%lf",&tmp);
	  if(got != 1)
	    {
	      fprintf(logfile,KRED "Error reading text file!\n" RESET);
	      fclose(file);
	      return 1;
	    }
	  p->a_ts[t][0] = tmp;
	}
      else if(t<TCONV0)
	{
	  p->a_ts[t][0] = p->a_ts[t-1][0]*p->ga[0];
	}
      else if(t<TCONV1)
	{
	  double frac = ((double)cnt)/(TCONV1-TCONV0);
	  double g = frac*p->ga_bgp + (1.0-frac)*p->ga[0];
	  p->a_ts[t][0] = p->a_ts[t-1][0]*g;
	  cnt=cnt+1;
	}
      else
	{
	  p->a_ts[t][0] = p->a_ts[t-1][0]*p->ga_bgp;
	}
    }
  fclose(file);

  //  other country data
  file = fopen("../python/output/rgdppc_ROW_gdp.txt","rb");
  for(t=0; t<NMATCH; t++)
    {
      got = fscanf(file,"%lf",&tmp);
      if(got != 1)
	{
	  fprintf(logfile,KRED "Error reading text file!\n" RESET);
	  fclose(file);
	  return 1;
	}
      p->RGDP_ts[t][1] = tmp;
    }
  fclose(file);

  file = fopen("../python/output/C_frac_ROW.txt","rb");
  for(t=0; t<NMATCH; t++)
    {
      got = fscanf(file,"%lf",&tmp);
      if(got != 1)
	{
	  fprintf(logfile,KRED "Error reading text file!\n" RESET);
	  fclose(file);
	  return 1;
	}
      p->C_frac_ts[t][1] = tmp;
    }
  fclose(file);

  file = fopen("../python/output/I_frac_ROW.txt","rb");
  for(t=0; t<NMATCH; t++)
    {
      got = fscanf(file,"%lf",&tmp);
      if(got != 1)
	{
	  fprintf(logfile,KRED "Error reading text file!\n" RESET);
	  fclose(file);
	  return 1;
	}
      p->I_frac_ts[t][1] = tmp;
    }
  fclose(file);

  file = fopen("../python/output/demo_ROW_gdp.txt","rb");
  cnt=0;
  for(t=0; t<NT+1; t++)
    {
      got = fscanf(file,"%lf %lf",&tmp1,&tmp2);
      if(got != 2)
	{
	  fprintf(logfile,KRED "Error reading text file!\n" RESET);
	  fclose(file);
	  return 1;
	}
      if(t<TCONV0)
	{
	  p->popt_ts[t][1] = tmp1;
	  p->popw_ts[t][1] = tmp2;
	}
      else if(t<TCONV1)
	{
	  double frac = ((double)cnt)/(TCONV1-TCONV0);
	  p->popt_ts[t][1] = p->popt_ts[t-1][1]*(1.0+frac*(tmp1/p->popt_ts[t-1][1]-1.0));
	  p->popw_ts[t][1] = p->popw_ts[t-1][1]*(1.0+frac*(tmp2/p->popw_ts[t-1][1]-1.0));
	  cnt=cnt+1;
	}
      else
	{
	  p->popt_ts[t][1] = p->popt_ts[t-1][1];
	  p->popw_ts[t][1] = p->popw_ts[t-1][1];
	}
    }
  fclose(file);

  cnt=0;
  file = fopen("../python/output/lp_ROW_gdp.txt","rb");
  for(t=0; t<NT+1; t++)
    {
      if(t<NMATCH)
	{
	  got = fscanf(file,"%lf",&tmp);
	  if(got != 1)
	    {
	      fprintf(logfile,KRED "Error reading text file!\n" RESET);
	      fclose(file);
	      return 1;
	    }
	  p->a_ts[t][1] = tmp;
	}
      else if(t<TCONV0)
	{
	  p->a_ts[t][1] = p->a_ts[t-1][1]*p->ga[1];
	}
      else if(t<TCONV1)
	{
	  double frac = ((double)cnt)/(TCONV1-TCONV0);
	  double g = frac*p->ga_bgp + (1.0-frac)*p->ga[1];
	  p->a_ts[t][1] = p->a_ts[t-1][1]*g;
	  cnt=cnt+1;
	}
      else
	{
	  p->a_ts[t][1] = p->a_ts[t-1][1]*p->ga_bgp;
	}
    }
  fclose(file);

  uint i;
  for(i=0; i<NC; i++)
    {
      for(t=0; t<NT+1; t++)
	{
	  p->pope_ts[t][i] = (p->popw_ts[t][i]*0.5608+0.5*(p->popt_ts[t][i]-p->popw_ts[t][i]*0.5608))/1.5608;
	}
      double tmp = p->pope_ts[0][i];
      for(t=0; t<NT+1; t++)
	{
	  p->pope_ts[t][i] = p->pope_ts[t][i]/tmp;
	}
    }

  if(no_demo_flag)
    {
      for(i=0; i<NC; i++)
	{
	  for(t=0; t<NT+1; t++)
	    {
	      p->popt_ts[t][i] = 1.0;
	      p->popw_ts[t][i] = 1.0;
	      p->pope_ts[t][i] = 1.0;
	    }
	}
    }
  
  if(sym_g_flag)
    {
      for(i=0; i<NC; i++)
	{
	  double ah=1.0;
	  for(t=0; t<NT+1; t++)
	    {
	      if(t>0)
		{
		  ah = ah*p->ga_bgp;
		}
	      p->a_ts[t][i] = ah;
	    }
	}
    }

  if(homotopy)
    {
      for(i=0; i<NC; i++)
	{
	  double ah=1.0;
	  for(t=0; t<NT+1; t++)
	    {
	      if(t>0)
		{
		  ah = ah*p->ga_bgp;
		}
	      p->a_ts[t][i] = ah*hfrac + p->a_ts[t][i]*(1.0-hfrac);
	    }
	}      
    }

  for(i=0; i<NC; i++)
      {
	p->ga_ts[0][i] = p->ga_bgp;
	for(t=1; t<NT+1; t++)
	  {
	    p->ga_ts[t][i] = p->a_ts[t][i]/p->a_ts[t-1][i];
	  }
      }

  return 0;
}

uint store_base_period_values(params * p)
{
  double mkt_clear_tol = 1.0e-7;
  uint varow = NC*NS;
  uint gorow = NC*NS+1;

  SET_ALL_V(p->y0,NC,0.0);
  SET_ALL_V(p->va0,NC,0.0);
  SET_ALL_V(p->k0,NC,0.0);
  SET_ALL_V(p->l0,NC,0.0);
  SET_ALL_V(p->md0,NC,0.0);
  SET_ALL_V(p->m0,NC,0.0);
  SET_ALL_V(p->m02,NC*NC,0.0);
  SET_ALL_V(p->q0,NC,0.0);
  SET_ALL_V(p->q02,NC*NC,0.0);
  SET_ALL_V(p->ex0,NC*NC,0.0);
  SET_ALL_V(p->im0,NC*NC,0.0);
  SET_ALL_V(p->nx0,NC*NC,0.0);
  SET_ALL_V(p->c0,NC,0.0);
  SET_ALL_V(p->i0,NC,0.0);
  SET_ALL_V(p->ii0,NC,0.0);  

  uint i, s, j, r;

  for(i=0; i<NC; i++)
    {
      uint ccol = NC*NS+i;
      uint icol = NC*NS+NC+i;

      for(s=0; s<NS; s++)
	{  
	  // first get value added and factors
	  uint scol = i*NS + s;
	  p->y0[i] = p->y0[i] + p->iomat[gorow][scol];
	  p->va0[i] = p->va0[i] + p->iomat[varow][scol];

	  // now get demand for products from different source countries and sectors 
	  for(j=0; j<NC; j++)
	    {
	      p->i0[i] = p->i0[i] + p->iomat[j*NS+s][icol];
	      p->c0[i] = p->c0[i] + p->iomat[j*NS+s][ccol];
	      p->q02[i][j] = p->q02[i][j] + p->iomat[j*NS+s][ccol] + p->iomat[j*NS+s][icol];

	      for(r=0; r<NS; r++)
		{
		  uint rcol = i*NS + r;
		  if(noio_flag)
		    {
		      p->m02[i][j] = 0.0;
		      p->md0[i] = 0.0;
		    }
		  else
		    {
		      p->m02[i][j] = p->m02[i][j] + p->iomat[j*NS+s][rcol];
		      p->md0[i] = p->md0[i] + p->iomat[j*NS+s][rcol];
		    }
		}
	    }
	}

      p->q0[i] = sum(p->q02[i],NC);
      p->m0[i] = sum(p->m02[i],NC);
      p->l0[i] = p->l0[i] + (1.0 - p->alpha[i]) * p->va0[i];
      p->ll0[i] = p->l0[i];
      p->ii0[i] = p->i0[i];

    }

  // initial capital stocks and tax rates
  p->kk0[0] = 281.288;
  p->kk0[1] = 306.257 * p->va0[1] / 100.0; // trade-weighted: 299.7659
  //p->kk0[1] = 250. * sum(p->va0[1],NS)/100.0;

  for(i=0; i<NC; i++)
    {
      double rdky = p->alpha[i] * p->va0[i];
      double dky = p->delta * p->kk0[i];
      p->r0[i] = ((1.0-p->tauk[i])*rdky-dky)/p->kk0[i];
      p->k0[i] = p->kk0[i];
    }

  for(i=0; i<NC; i++)
    {
      for(j=0; j<NC; j++)
	{
	  if(j != i)
	    {
	      p->im02[i][j] = p->q02[i][j] + p->m02[i][j];
	      p->ex02[i][j] = p->q02[j][i] + p->m02[j][i];
	      p->nx02[i][j] = p->ex02[i][j] - p->im02[i][j];
	    }
	  else
	    {
	      p->im02[i][j] = 0.0;
	      p->ex02[i][j] = 0.0;
	      p->nx02[i][j] = 0.0;
	    }
	  p->im0[i][j] = p->im02[i][j];
	  p->ex0[i][j] = p->ex02[i][j];
	  p->nx0[i][j] = p->nx02[i][j];
	}
    }

  double tmp=0.0;
  for(i=0; i<NC; i++)
    {
      tmp = (p->va0[i] - (p->q0[i] + sum(p->ex0[i],NC) - sum(p->im0[i],NC)))/p->va0[i];
      if(fabs(tmp)>mkt_clear_tol)
	{
	  fprintf(logfile,KRED "GDP != C+I+NX for country %d, error = %f\n" RESET,i,tmp);
	  return 1;
	}

      tmp = p->y0[i];
      for(j=0; j<NC; j++)
	{
	  tmp = tmp - p->q02[j][i] - p->m02[j][i];
	}
      tmp = tmp/p->y0[i];
      if(fabs(tmp)>mkt_clear_tol)
	{
	  fprintf(logfile,KRED "supply != demand for country %d, error = %f\n" RESET,i,tmp);
	  return 1;
	}

      tmp = (p->y0[i] - (p->va0[i] + p->md0[i]))/p->y0[i];
      if(fabs(tmp)>mkt_clear_tol)
	{
	  fprintf(logfile,KRED "go != va + m for country %d, error = %f\n" RESET,i,tmp);
	  return 1;
	}

      tmp = p->m0[i] - SUM(p->m02[i],NC);
      if(fabs(tmp)>mkt_clear_tol)
	{
	  fprintf(logfile,KRED "m != sum(m2) for country %d, error = %f\n" RESET,i,tmp);
	  return 1;
	}
    }

  // from Lane and Milesi-Ferretti (2007)
  double usa_nfa_frac = -0.072;
  //double chn_nfa_frac = -0.086;
  
  p->b0[0] = 100*usa_nfa_frac;
  p->b0[1] = -p->b0[0];

  return 0;

}

uint calibrate_prod_params(params * p)
{
  uint i;
  double tmp;
  
  SET_ALL_V(p->lam,NC,0.0);
  SET_ALL_V(p->A,NC,0.0);

  for(i=0; i<NC; i++)
    {
      p->A[i] = p->va0[i] / ( pow(p->k0[i],p->alpha[i])*pow(p->l0[i],1.0-p->alpha[i]) );	 
      p->lam_va[i] = p->va0[i] / p->y0[i];
      p->lam[i] = p->md0[i] / p->y0[i];
    }

  for(i=0; i<NC; i++)
    {
      if(noio_flag==0)
	{
	  tmp = prod_go(p->va0[i],p->md0[i],p->lam_va[i],p->lam[i]) - p->y0[i];
	  if(fabs(tmp)>TINY)
	    {
	      fprintf(logfile,KRED "prod_go != y0 for country %d, error = %f" RESET,i,tmp);
	      return 1;
	    }
	}

      tmp = prod_va(p->k0[i],p->l0[i],p->A[i],p->alpha[i],1.0) - p->va0[i];
      if(fabs(tmp)>TINY)
	{
	  fprintf(logfile,KRED "prod_va != va0 for country/sector %d, error = %f" RESET,i,tmp);
	  return 1;
	}

      tmp = p->y0[i] - p->va0[i] - p->md0[i];
      if(fabs(tmp)>TINY)
	{
	  fprintf(logfile,KRED "nonzero profits for country %d, error = %f" RESET,i,tmp);
	  return 1;
	}

      tmp = (1.0-p->lam[i]) * (1.0-p->alpha[i]) * p->A[i]/p->lam_va[i] *
	pow(p->k0[i],p->alpha[i]) * pow(p->l0[i],-(p->alpha[i])) - 1.0;
      if(fabs(tmp)>TINY)
	{
	  fprintf(logfile,KRED "labor FOC for country %d, error = %f" RESET,i,tmp);
	  return 1;
	}

      tmp = (1.0-p->lam[i]) * (p->alpha[i]) * p->A[i]/p->lam_va[i] *
	pow(p->k0[i],p->alpha[i]-1.0) * pow(p->l0[i],1.0-p->alpha[i]) - 
	(p->r0[i] + p->delta)/(1.0-p->tauk[i]);
      if(fabs(tmp)>TINY)
	{
	  fprintf(logfile,KRED "capital FOC for country %d, error = %f" RESET,i,tmp);
	  return 1;
	}
    }

  return 0;
}

uint calibrate_fin_params(params * p)
{
  uint i,j,jj,cnt,idx;
  double tmp;
  double tmp1[NC];

  SET_ALL_V(p->mu,NC*NC,0.0);
  SET_ALL_V(p->M,NC,0.0);
  SET_ALL_V(p->H,NC,0.0);
  SET_ALL_V(p->theta,NC*NC,0.0);

  if(noio_flag==0)
    {
      for(i=0; i<NC; i++)
	{
	  idx=i;
	  for(j=0; j<NC; j++)
	    {
	      tmp1[j] = pow(p->m02[i][j]/p->m02[i][idx],1.0 - p->zeta[i]);
	    }
	  p->mu[i][idx] = 1.0/sum(tmp1,NC);
	  cnt=0;
	  for(j=0; j<NC; j++)
	    {
	      cnt=cnt+1;
	      if(j != idx)
		{
		  if(cnt<NC)
		    {
		      p->mu[i][j] = p->mu[i][idx]*tmp1[j];
		    }
		  else
		    {
		      p->mu[i][j] = 1.0 - sum(p->mu[i],NC-1);
		    }
		}
	    
	    }
	  tmp = pow( DOT_PROD_EX(p->m02[i],p->mu[i],NC,p->zeta[i]), 1.0/p->zeta[i] );
	  p->M[i] = p->m0[i]/tmp;
	}
    }

  for(i=0; i<NC; i++)
    {
      idx=i;
      for(j=0; j<NC; j++)
	{
	  tmp1[j] = pow(p->q02[i][j]/p->q02[i][idx],1.0 - p->sig[i]);
	}
      p->theta[i][idx] = 1.0/sum(tmp1,NC);
      cnt=0;
      for(j=0; j<NC; j++)
	{
	  cnt=cnt+1;
	  if(j != idx)
	    {
	      if(cnt<NC)
		{
		  p->theta[i][j] = p->theta[i][idx]*tmp1[j];
		}
	      else
		{
		  p->theta[i][j] = 1.0 - sum(p->theta[i],NC-1);
		}
	    }
	    
	}
      tmp = pow( DOT_PROD_EX(p->q02[i],p->theta[i],NC,p->sig[i]), 1.0/p->sig[i] );
      p->H[i] = p->q0[i]/tmp;
    }

  for(i=0; i<NC; i++)
    {
      if(noio_flag==0)
	{
	  tmp = p->m0[i] - prod_m(p->m02[i],p->M[i],p->mu[i],p->zeta[i]);
	  if(fabs(tmp)>TINY)
	    {
	      fprintf(logfile,KRED "Intermediate Armington production function for country %d, error = %f" RESET,i,tmp);
	      return 1;
	    }
	}

      tmp = p->q0[i] - prod_q(p->q02[i],p->H[i],p->theta[i],p->sig[i]);
      if(fabs(tmp)>TINY)
	{
	  fprintf(logfile,KRED "Final Armington production function for country/sector %d, error = %f" RESET,i,tmp);
	  return 1;
	}

      for(j=0; j<NC; j++)
	{
	  if(noio_flag==0)
	    {
	      tmp = 1.0 - p->mu[i][j] * pow(p->M[i],p->zeta[i]) * 
		pow(p->m0[i]/p->m02[i][j],1.0-p->zeta[i]);

	      if(fabs(tmp)>TINY)
		{
		  fprintf(logfile,KRED "Intermediate Armington FOC for country/country %d/%d, error = %f" RESET, i,j,tmp);
		  return 1;
		}
	    }

	  tmp = 1.0 - p->theta[i][j] * pow(p->H[i],p->sig[i]) * 
	    pow(p->q0[i]/p->q02[i][j],1.0-p->sig[i]);

	  if(fabs(tmp)>TINY)
	    {
	      fprintf(logfile,KRED "Final Armington FOC for country/country %d/%d, error = %f" RESET, i,j,tmp);
	      return 1;
	    }

	}

      for(j=0; j<NC; j++)
	{
	  for(jj=0; jj<NC; jj++)
	    {
	      if(noio_flag==0)
		{
		  tmp = 1.0 - (p->mu[i][j] / p->mu[i][jj]) * 
		    pow(p->m02[i][jj]/p->m02[i][j],1.0-p->zeta[i]);
		  if(fabs(tmp)>TINY)
		    {
		      fprintf(logfile,KRED "Intermediate Armington FOC v2 for country/country/country %d/%d/%d, error = %f" RESET,
			      i,j,jj,tmp);
		      return 1;
		    }
		}

	      tmp = 1.0 - (p->theta[i][j] / p->theta[i][jj]) * 
		pow(p->q02[i][jj]/p->q02[i][j],1.0-p->sig[i]);
	      if(fabs(tmp)>TINY)
		{
		  fprintf(logfile,KRED "Final Armington FOC v2 for country/country/country %d/%d/%d, error = %f" RESET,
			  i,j,jj,tmp);
		  return 1;
		}
	    }
	}
    }

  return 0;
  
}

uint calibrate_hh_params(params * p)
{
  uint i, t;
  double tmp;
  double tmp1;
  
  for(i=0; i<NC; i++)
    {
      p->lbar0[i] = 3.0 * p->ll0[i];

      if(fixl==1)
	{
	  p->phi[i]=1.0;
	}
      else
	{
	  tmp = p->c0[i]/(p->lbar0[i] - p->ll0[i]);
	  p->phi[i] = tmp/(1.0+tmp);
	}

      if(fixl==0)
	{
	  tmp = muc(p->c0[i],p->ll0[i],p->lbar0[i],1.0,1.0,p->phi[i], p->psi) /
	    mul(p->c0[i],p->ll0[i],p->lbar0[i],1.0,1.0,p->phi[i],p->psi) - 1.0;
	  if(fabs(tmp)>TINY)
	    {
	      fprintf(logfile,KRED "HH intratemp FOC 2 for country %d, error = %f" RESET,i,tmp);
	      return 1;
	    }
	}

      tmp1 = p->c0[i]*p->ga_bgp;
      // note: beta = (1.0+rss)/gbgp^(phi*psi-1) > 1
      // but, as long as beta*gbgp^(phi*psi) < 1, we can calculate welfare just fine
      p->beta[i] = muc(p->c0[i],p->ll0[i],p->lbar0[i],1.0,1.0,p->phi[i], p->psi) /
	muc(tmp1,p->ll0[i],p->lbar0[i],1.0,1.0,p->phi[i], p->psi) / (1.0 + p->rss);
    }

  // population growth growth
  for(i=0; i<NC; i++)
    {
      for(t=0; t<NT+1; t++)
	{
	  p->lbar_ts[t][i] = p->lbar0[i] * p->popw_ts[t][i];
	}
    }

  return 0;
}

uint calibrate(params * p)
{
  set_nontargeted_params(p);

  if(load_iomat(p))
    {
      return 1;
    }

  load_ts_params(p);

  if(store_base_period_values(p))
    {
      return 1;
    }

  if(calibrate_prod_params(p))
    {
      return 1;
    }

  if(calibrate_fin_params(p))
    {
      return 1;
    }

  if(calibrate_hh_params(p))
    {
      return 1;
    }

  return 0;
}

uint write_params(const params * p)
{
  
  FILE * file;
  file = fopen("output/params_2c_1sector.txt","wb");

  if(file)
    {
      fprintf(file,"Scalar parameters:\n");
      fprintf(file,"delta: %0.4f\n",p->delta);
      fprintf(file,"rss: %0.4f\n",p->rss);
      fprintf(file,"psi: %0.4f\n",p->psi);
      
      fprintf(file,"\nVECTOR PARAMETERS (1 x NC):\n");

      fprintf(file,"beta:");
      fprintf_vec(file,p->beta,NC);

      fprintf(file,"tauk:");
      fprintf_vec(file,p->tauk,NC);

      fprintf(file,"r0:");
      fprintf_vec(file,p->r0,NC);

      fprintf(file,"phi:");
      fprintf_vec(file,p->phi,NC);

      fprintf(file,"lbar:");
      fprintf_vec(file,p->lbar0,NC);

      fprintf(file,"kk0:");
      fprintf_vec(file,p->kk0,NC);

      fprintf(file,"b0:");
      fprintf_vec(file,p->b0,NC);

      fprintf(file,"\nVECTOR PARAMETERS (1x  NC):\n\n");
      
      fprintf(file,"sig:\n");
      fprintf_vec(file,p->sig,NC);

      fprintf(file,"H:\n");
      fprintf_vec(file,p->H,NC);

      fprintf(file,"zeta:\n");
      fprintf_vec(file,p->zeta,NC);

      fprintf(file,"M:\n");
      fprintf_vec(file,p->M,NC);

      fprintf(file,"lam_va:\n");
      fprintf_vec(file,p->lam_va,NC);

      fprintf(file,"lam:\n");
      fprintf_vec(file,p->lam,NC);

      fprintf(file,"alpha:\n");
      fprintf_vec(file,p->alpha,NC);

      fprintf(file,"A:\n");
      fprintf_vec(file,p->A,NC);

      fprintf(file,"\nMATRIX PARAMETERS:\n\n");

      fprintf(file,"theta (NC x NC):\n");
      fprintf_mat_3(file,p->theta,NC);

      fprintf(file,"mu (NC x NC):\n");
      fprintf_mat_3(file,p->mu,NC);

      fclose(file);
      return 0;
    }
  else
    {
      fprintf(logfile,KRED "Error opening file to write parameters!\n" RESET);
      return 1;
    }
}

#endif
