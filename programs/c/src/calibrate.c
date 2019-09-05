#ifndef __CALIBRATE_C__
#define __CALIBRATE_C__

#include "calibrate.h"

void set_nontargeted_params(params * p)
{
  uint i;

  // parameters common across countries

  if(low_r_flag)
    {
      p->rss=0.005;
    }
  else
    {
      p->rss = 0.025;
    }
  p->ga_bgp = 1.02;
  p->conv_speed = 0.75;
  p->delta = 0.06;

  if(high_etaK_flag)
    {
      p->etaK = 0.9;
    }
  else if(low_etaK_flag)
    {
      p->etaK = 0.7;
    }
  else
    {
      p->etaK = 0.8;
    }

  p->etaM = 0.00;
  p->etaF = 0.00;
  p->rho = 1.0-1.0/0.65;
  p->psi = -1.0;

  if(high_rhow_flag)
    {
      p->wedge_speed = 0.8;
    }
  else if(low_rhow_flag)
    {
      p->wedge_speed = 0.5;
    }
  else
    {
      p->wedge_speed = 0.65;
    }


  SET_ALL_V(p->alpha,NC*NS,0.34);

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
	  p->sig[i][0] = 1.0-1.0/6.0;
	  p->sig[i][1] = 1.0-1.0/1.01;
	}
      else if(high_arm_flag)
	{
	  p->sig[i][0] = 1.0-1.0/4.0;
	  p->sig[i][1] = 1.0-1.0/1.01; 
	}
      else if(low_arm_flag)
	{
	  p->sig[i][0] = 1.0-1.0/1.5;
	  p->sig[i][1] = 1.0-1.0/1.01; 
	}
      else
	{
	  p->sig[i][0] = 1.0-1.0/2.0;
	  p->sig[i][1] = 1.0-1.0/1.01;
	}

      if(high_arm_flag)
	{
	  p->zeta[i][0] = 1.0-1.0/6.0;
	  p->zeta[i][1] = 1.0-1.0/1.01;
	}
      else if(low_arm_flag)
	{
	  p->zeta[i][0] = 1.0-1.0/2.01;
	  p->zeta[i][1] = 1.0-1.0/1.01;      	  
	}
      else
	{
	  p->zeta[i][0] = 1.0-1.0/3.0;
	  p->zeta[i][1] = 1.0-1.0/1.01;      
	}
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

  SET_ALL_V(p->y0,NC*NS,0.0);
  SET_ALL_V(p->va0,NC*NS,0.0);
  SET_ALL_V(p->k0,NC*NS,0.0);
  SET_ALL_V(p->l0,NC*NS,0.0);
  SET_ALL_V(p->md0,NC*NS*(NS-1),0.0);
  SET_ALL_V(p->m0,NC*(NS-1),0.0);
  SET_ALL_V(p->m02,NC*(NS-1)*NC,0.0);
  SET_ALL_V(p->q0,NC*(NS-1),0.0);
  SET_ALL_V(p->q02,NC*(NS-1)*NC,0.0);
  SET_ALL_V(p->ex0,NC*NC,0.0);
  SET_ALL_V(p->im0,NC*NC,0.0);
  SET_ALL_V(p->nx0,NC*NC,0.0);
  SET_ALL_V(p->c0,NC*(NS-1),0.0);
  SET_ALL_V(p->i0,NC*NS,0.0);
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
	  p->y0[i][s] = p->iomat[gorow][scol];
	  p->va0[i][s] = p->iomat[varow][scol];

	  p->l0[i][s] = (1.0 - p->alpha[i][s]) * p->va0[i][s];

	  // now get demand for products from different source countries and sectors 
	  for(j=0; j<NC; j++)
	    {
	      p->i0[i][s] = p->i0[i][s] + p->iomat[j*NS+s][icol];
	      if(s<NS-1)
		{
		  p->c0[i][s] = p->c0[i][s] + p->iomat[j*NS+s][ccol];
		  p->q02[i][s][j] = p->iomat[j*NS+s][ccol] + p->iomat[j*NS+s][icol];

		  for(r=0; r<NS; r++)
		    {
		      uint rcol = i*NS + r;
		      if(noio_flag)
			{
			  p->m02[i][s][j] = 0.0;
			  p->md0[i][r][s] = 0.0;
			}
		      else
			{
			  p->m02[i][s][j] = p->m02[i][s][j] + p->iomat[j*NS+s][rcol];
			  p->md0[i][r][s] = p->md0[i][r][s] + p->iomat[j*NS+s][rcol];
			}
		    }
		  p->q0[i][s] = sum(p->q02[i][s],NC);
		  p->m0[i][s] = sum(p->m02[i][s],NC);
		}
	    }
	}

      p->ll0[i] = sum(p->l0[i],NS);
      p->ii0[i] = sum(p->i0[i],NS);

    }

  // initial capital stocks and tax rates
  p->kk0[0] = 281.288;
  p->kk0[1] = 306.257 * sum(p->va0[1],NS)/100.0; // trade-weighted: 299.7659
  //p->kk0[1] = 250. * sum(p->va0[1],NS)/100.0;

  for(i=0; i<NC; i++)
    {
      double rdky = p->alpha[i][0] * sum(p->va0[i],NS);
      double dky = p->delta * p->kk0[i];
      //double rky = rdky - dky;
      p->r0[i] = ((1.0-p->tauk[i])*rdky-dky)/p->kk0[i];
      //p->tauk[i] = 1.0 - (p->r0[i] + p->delta)*p->kk0[i]/rdky;

      for(s=0; s<NS; s++)
	{
	  p->k0[i][s] = p->alpha[i][s] * p->va0[i][s] / ((p->r0[i] + p->delta) / (1.0 - p->tauk[i]));
	}

      double tmp = fabs(p->kk0[i] - sum(p->k0[i],NS));
      if(fabs(tmp)>mkt_clear_tol)
	{
	  fprintf(logfile,KRED "kk0 != sum(k0), for i = %d, error = %f\n" RESET,i,tmp);
	  return 1;
	}
    }

  for(i=0; i<NC; i++)
    {
      for(j=0; j<NC; j++)
	{
	  for(s=0; s<NS-1; s++)
	    {
	      if(j != i)
		{
		  p->im02[i][s][j] = p->q02[i][s][j] + p->m02[i][s][j];
		  p->ex02[i][s][j] = p->q02[j][s][i] + p->m02[j][s][i];
		  p->nx02[i][s][j] = p->ex02[i][s][j] - p->im02[i][s][j];
		}
	      else
		{
		  p->im02[i][s][j] = 0.0;
		  p->ex02[i][s][j] = 0.0;
		  p->nx02[i][s][j] = 0.0;
		}
	    }
	  p->im0[i][j] = p->im02[i][0][j] + p->im02[i][1][j];
	  p->ex0[i][j] = p->ex02[i][0][j] + p->ex02[i][1][j];
	  p->nx0[i][j] = p->nx02[i][0][j] + p->nx02[i][1][j];
	}
    }

  double tmp=0.0;
  for(i=0; i<NC; i++)
    {
      tmp = (sum(p->va0[i],NS) - (sum(p->q0[i],NS-1) + 
				  p->i0[i][NS-1] + 
				  sum(p->ex0[i],NC) - 
				  sum(p->im0[i],NC)))/sum(p->va0[i],NS);
      if(fabs(tmp)>mkt_clear_tol)
	{
	  fprintf(logfile,KRED "GDP != C+I+NX for country %d, error = %f\n" RESET,i,tmp);
	  return 1;
	}

      for(s=0; s<NS-1; s++)
	{
	  tmp = p->y0[i][s];
	  for(j=0; j<NC; j++)
	    {
	      tmp = tmp - p->q02[j][s][i] - p->m02[j][s][i];
	    }
	  tmp = tmp/p->y0[i][s];
	  if(fabs(tmp)>mkt_clear_tol)
	    {
	      fprintf(logfile,KRED "supply != demand for country/sector %d/%d, error = %f\n" RESET,i,s,tmp);
	      return 1;
	    }
	}

      s=NS-1;
      tmp = (p->y0[i][NS-1] - p->i0[i][NS-1])/p->y0[i][NS-1];
      if(fabs(tmp)>mkt_clear_tol)
	{
	  fprintf(logfile,KRED "supply != demand for country/sector %d/%d, error = %f\n" RESET,i,s,tmp);
	  return 1;
	}

      for(s=0; s<NS; s++)
	{
	  tmp = p->y0[i][s] - (p->va0[i][s] + SUM(p->md0[i][s],NS-1));
	  if(fabs(tmp)>mkt_clear_tol)
	    {
	      fprintf(logfile,KRED "go != va + m for country/sector %d/%d, error = %f\n" RESET,i,s,tmp);
	      return 1;
	    }

	  if(s<NS-1)
	    {
	      tmp = p->m0[i][s] - SUM(p->m02[i][s],NC);
	      if(fabs(tmp)>mkt_clear_tol)
		{
		  fprintf(logfile,KRED "m != sum(m2) for country/sector %d/%d, error = %f\n" RESET,i,s,tmp);
		  return 1;
		}
	    }
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
  uint i,s,r;
  double tmp;
  
  SET_ALL_V(p->lam,NC*NS*(NS-1),0.0);
  SET_ALL_V(p->A,NC*NS,0.0);

  for(i=0; i<NC; i++)
    {
      for(s=0; s<NS; s++)
	{
	  p->A[i][s] = p->va0[i][s] / ( pow(p->k0[i][s],p->alpha[i][s])*pow(p->l0[i][s],1.0-p->alpha[i][s]) );	 
	  p->lam_va[i][s] = p->va0[i][s] / p->y0[i][s];
	  for(r=0; r<NS-1; r++)
	    {
	      p->lam[i][s][r] = p->md0[i][s][r] / p->y0[i][s];
	    }
	}
    }

  for(i=0; i<NC; i++)
    {
      for(s=0; s<NS; s++)
	{
	  if(noio_flag==0)
	    {
	      tmp = prod_go(p->va0[i][s],p->md0[i][s],p->lam_va[i][s],p->lam[i][s]) - p->y0[i][s];
	      if(fabs(tmp)>TINY)
		{
		  fprintf(logfile,KRED "prod_go != y0 for country/sector %d/%d, error = %f" RESET,i,s,tmp);
		  return 1;
		}
	    }

	  tmp = prod_va(p->k0[i][s],p->l0[i][s],p->A[i][s],p->alpha[i][s],1.0) - p->va0[i][s];
	  if(fabs(tmp)>TINY)
	    {
	      fprintf(logfile,KRED "prod_va != va0 for country/sector %d/%d, error = %f" RESET,i,s,tmp);
	      return 1;
	    }

	  tmp = p->y0[i][s] - p->va0[i][s] - sum(p->md0[i][s],NS-1);
	  if(fabs(tmp)>TINY)
	    {
	      fprintf(logfile,KRED "nonzero profits for country/sector %d/%d, error = %f" RESET,i,s,tmp);
	      return 1;
	    }

	  tmp = (1.0-sum(p->lam[i][s],NS-1)) * (1.0-p->alpha[i][s]) * p->A[i][s]/p->lam_va[i][s] *
	    pow(p->k0[i][s],p->alpha[i][s]) * pow(p->l0[i][s],-(p->alpha[i][s])) - 1.0;
	  if(fabs(tmp)>TINY)
	    {
	      fprintf(logfile,KRED "labor FOC for country/sector %d/%d, error = %f" RESET,i,s,tmp);
	      return 1;
	    }

	  tmp = (1.0-sum(p->lam[i][s],NS-1)) * (p->alpha[i][s]) * p->A[i][s]/p->lam_va[i][s] *
	    pow(p->k0[i][s],p->alpha[i][s]-1.0) * pow(p->l0[i][s],1.0-p->alpha[i][s]) - 
	    (p->r0[i] + p->delta)/(1.0-p->tauk[i]);
	  if(fabs(tmp)>TINY)
	    {
	      fprintf(logfile,KRED "capital FOC for country/sector %d/%d, error = %f" RESET,i,s,tmp);
	      return 1;
	    }
	}
    }

  return 0;
}

uint calibrate_fin_params(params * p)
{
  uint i,s,j,jj,cnt,idx;
  double tmp;
  double tmp1[NC];

  SET_ALL_V(p->mu,NC*(NS-1)*NC,0.0);
  SET_ALL_V(p->M,NC*(NS-1),0.0);
  SET_ALL_V(p->G,NC,0.0);
  SET_ALL_V(p->H,NC*(NS-1),0.0);
  SET_ALL_V(p->eps,NC*NF*NS,0.0);
  SET_ALL_V(p->theta,NC*(NS-1)*NC,0.0);

  if(noio_flag==0)
    {
      for(i=0; i<NC; i++)
	{
	  for(s=0; s<NS-1; s++) // no construction here to worry about
	    {
	      idx=i;
	      for(j=0; j<NC; j++)
		{
		  tmp1[j] = pow(p->m02[i][s][j]/p->m02[i][s][idx],1.0 - p->zeta[i][s]);
		}
	      p->mu[i][s][idx] = 1.0/sum(tmp1,NC);
	      cnt=0;
	      for(j=0; j<NC; j++)
		{
		  cnt=cnt+1;
		  if(j != idx)
		    {
		      if(cnt<NC)
			{
			  p->mu[i][s][j] = p->mu[i][s][idx]*tmp1[j];
			}
		      else
			{
			  p->mu[i][s][j] = 1.0 - sum(p->mu[i][s],NC-1);
			}
		    }
	    
		}
	      tmp = pow( DOT_PROD_EX(p->m02[i][s],p->mu[i][s],NC,p->zeta[i][s]), 1.0/p->zeta[i][s] );
	      p->M[i][s] = p->m0[i][s]/tmp;
	    }
	}
    }

  for(i=0; i<NC; i++)
    {
      for(s=0; s<NS-1; s++)
	{
	  idx=i;
	  for(j=0; j<NC; j++)
	    {
	      tmp1[j] = pow(p->q02[i][s][j]/p->q02[i][s][idx],1.0 - p->sig[i][s]);
	    }
	  p->theta[i][s][idx] = 1.0/sum(tmp1,NC);
	  cnt=0;
	  for(j=0; j<NC; j++)
	    {
	      cnt=cnt+1;
	      if(j != idx)
		{
		  if(cnt<NC)
		    {
		      p->theta[i][s][j] = p->theta[i][s][idx]*tmp1[j];
		    }
		  else
		    {
		      p->theta[i][s][j] = 1.0 - sum(p->theta[i][s],NC-1);
		    }
		}
	    
	    }
	  tmp = pow( DOT_PROD_EX(p->q02[i][s],p->theta[i][s],NC,p->sig[i][s]), 1.0/p->sig[i][s] );
	  p->H[i][s] = p->q0[i][s]/tmp;
	}
    }

  for(i=0; i<NC; i++)
    {
      for(s=0; s<NS; s++)
	{
	  p->eps[i][1][s] = p->i0[i][s] / p->ii0[i];
	}
      p->G[i] = p->ii0[i]/ ( pow(p->i0[i][0],p->eps[i][1][0]) * 
			     pow(p->i0[i][1],p->eps[i][1][1]) *
			     pow(p->i0[i][2],p->eps[i][1][2]));
    }

  for(i=0; i<NC; i++)
    {
      for(s=0; s<NS-1; s++)
	{
	  if(noio_flag==0)
	    {
	      tmp = p->m0[i][s] - prod_m(p->m02[i][s],p->M[i][s],p->mu[i][s],p->zeta[i][s]);
	      if(fabs(tmp)>TINY)
		{
		  fprintf(logfile,KRED "Intermediate Armington production function for country/sector %d/%d, error = %f" RESET,i,s,tmp);
		  return 1;
		}
	    }

	  tmp = p->q0[i][s] - prod_q(p->q02[i][s],p->H[i][s],p->theta[i][s],p->sig[i][s]);
	  if(fabs(tmp)>TINY)
	    {
	      fprintf(logfile,KRED "Final Armington production function for country/sector %d/%d, error = %f" RESET,i,s,tmp);
	      return 1;
	    }

	  for(j=0; j<NC; j++)
	    {
	      if(noio_flag==0)
		{
		  tmp = 1.0 - p->mu[i][s][j] * pow(p->M[i][s],p->zeta[i][s]) * 
		    pow(p->m0[i][s]/p->m02[i][s][j],1.0-p->zeta[i][s]);

		  if(fabs(tmp)>TINY)
		    {
		      fprintf(logfile,KRED "Intermediate Armington FOC for country/sectorcountry %d/%d/%d, error = %f" RESET, i,s,j,tmp);
		      return 1;
		    }
		}

	      tmp = 1.0 - p->theta[i][s][j] * pow(p->H[i][s],p->sig[i][s]) * 
		pow(p->q0[i][s]/p->q02[i][s][j],1.0-p->sig[i][s]);

	      if(fabs(tmp)>TINY)
		{
		  fprintf(logfile,KRED "Final Armington FOC for country/sector/country %d/%d/%d, error = %f" RESET, i,s,j,tmp);
		  return 1;
		}

	    }

	  for(j=0; j<NC; j++)
	    {
	      for(jj=0; jj<NC; jj++)
		{
		  if(noio_flag==0)
		    {
		      tmp = 1.0 - (p->mu[i][s][j] / p->mu[i][s][jj]) * 
			pow(p->m02[i][s][jj]/p->m02[i][s][j],1.0-p->zeta[i][s]);
		      if(fabs(tmp)>TINY)
			{
			  fprintf(logfile,KRED "Intermediate Armington FOC v2 for country/sector/country/country %d/%d/%d/%d, error = %f" RESET,
				  i,s,j,j,tmp);
			  return 1;
			}
		    }

		  tmp = 1.0 - (p->theta[i][s][j] / p->theta[i][s][jj]) * 
		    pow(p->q02[i][s][jj]/p->q02[i][s][j],1.0-p->sig[i][s]);
		  if(fabs(tmp)>TINY)
		    {
		      fprintf(logfile,KRED "Final Armington FOC v2 for country/sector/country/country %d/%d/%d/%d, error = %f" RESET,
			      i,s,j,jj,tmp);
		      return 1;
		    }
		}
	    }

	  tmp = p->ii0[i] - prod_inv(p->i0[i],p->eps[i][1],p->G[i]);
	  if(fabs(tmp)>TINY)
	    {
	      fprintf(logfile,KRED "ii0 != prod_inv for country %d, error = %f" RESET,i,tmp);
	      return 1;
	    }

	  for(s=0; s<NS; s++)
	    {
	      tmp = 1.0 - p->eps[i][1][s]*(p->ii0[i]/p->i0[i][s]);
	      if(fabs(tmp)>TINY)
		{
		  fprintf(logfile,KRED "Investment FOC for country/sector %d/%d, error = %f" RESET,i,s,tmp);
		  return 1;
		}
	    }
	}
    }

  return 0;
  
}

uint calibrate_hh_params(params * p)
{
  uint i, s, t;
  double tmp;
  double tmp1[NS-1];
  
  for(i=0; i<NC; i++)
    {
      tmp = pow(p->c0[i][0]/p->c0[i][1],1.0 - p->rho);
      p->eps[i][0][0] = tmp/(1.0+tmp);
      p->eps[i][0][1] = 1.0 - p->eps[i][0][0];
      p->lbar0[i] = 3.0 * p->ll0[i];

      if(fixl==1)
	{
	  p->phi[i]=1.0;
	}
      else
	{
	  tmp = (p->eps[i][0][0] * pow(p->c0[i][0], p->rho) + p->eps[i][0][1]*pow(p->c0[i][1], p->rho)) /
	    (p->lbar0[i] - p->ll0[i]) / p->eps[i][0][0] / pow(p->c0[i][0], p->rho-1.0);
	  p->phi[i] = tmp/(1.0+tmp);
	}

      tmp = muc(p->c0[i],p->ll0[i],p->lbar0[i],1.0,1.0,p->eps[i][0],p->rho, p->phi[i], p->psi, 0) /
	muc(p->c0[i],p->ll0[i],p->lbar0[i],1.0,1.0,p->eps[i][0],p->rho, p->phi[i], p->psi, 1) - 1.0;
      if(fabs(tmp)>TINY)
	{
	  fprintf(logfile,KRED "HH intratemp FOC 1 for country %d, error = %f" RESET,i,tmp);
	  return 1;
	}

      if(fixl==0)
	{
	  tmp = muc(p->c0[i],p->ll0[i],p->lbar0[i],1.0,1.0,p->eps[i][0],p->rho, p->phi[i], p->psi, 0) /
	    mul(p->c0[i],p->ll0[i],p->lbar0[i],1.0,1.0,p->eps[i][0],p->rho,p->phi[i],p->psi) - 1.0;
	  if(fabs(tmp)>TINY)
	    {
	      fprintf(logfile,KRED "HH intratemp FOC 2 for country %d, error = %f" RESET,i,tmp);
	      return 1;
	    }
	}

      for(s=0; s<NS-1; s++)
	{
	  tmp1[s] = p->c0[i][s]*p->ga_bgp;
	}
      // note: beta = (1.0+rss)/gbgp^(phi*psi-1) > 1
      // but, as long as beta*gbgp^(phi*psi) < 1, we can calculate welfare just fine
      p->beta[i] = muc(p->c0[i],p->ll0[i],p->lbar0[i],1.0,1.0,p->eps[i][0],p->rho, p->phi[i], p->psi, 0) /
	muc(tmp1,p->ll0[i],p->lbar0[i],1.0,1.0,p->eps[i][0],p->rho, p->phi[i], p->psi, 0) / (1.0 + p->rss);
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
  file = fopen("output/params_2c.txt","wb");

  if(file)
    {
      fprintf(file,"Scalar parameters:\n");
      fprintf(file,"rho: %0.4f\n",p->rho);
      fprintf(file,"delta: %0.4f\n",p->delta);
      fprintf(file,"rss: %0.4f\n",p->rss);
      fprintf(file,"psi: %0.4f\n",p->psi);
      
      fprintf(file,"\nVECTOR PARAMETERS (1 x NC):\n");

      fprintf(file,"G:");
      fprintf_vec(file,p->G,NC);

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

      fprintf(file,"\nMATRIX PARAMETERS (NC x NS or NX x NS-1):\n\n");
      
      fprintf(file,"sig:\n");
      fprintf_mat_2(file,p->sig,NC);

      fprintf(file,"H:\n");
      fprintf_mat_2(file,p->H,NC);

      fprintf(file,"zeta:\n");
      fprintf_mat_2(file,p->zeta,NC);

      fprintf(file,"M:\n");
      fprintf_mat_2(file,p->M,NC);

      fprintf(file,"lam_va:\n");
      fprintf_mat_1(file,p->lam_va,NC);

      fprintf(file,"alpha:\n");
      fprintf_mat_1(file,p->alpha,NC);

      fprintf(file,"A:\n");
      fprintf_mat_1(file,p->A,NC);

      fprintf(file,"\n3D PARAMETERS:\n\n");

      fprintf(file,"eps (NC x 2 x NS):\n");
      fprintf_3d_1(file,p->eps,NC);

      fprintf(file,"theta (NC x NS-1 x NC):\n");
      fprintf_3d_2(file,p->theta,NC);

      fprintf(file,"mu (NC x NS-1 x NC):\n");
      fprintf_3d_2(file,p->mu,NC);

      fprintf(file,"lam (NC x NS x NS-1):\n");
      fprintf_3d_3(file,p->lam,NC);

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
