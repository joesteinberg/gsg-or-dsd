#ifndef __EQM_C__
#define __EQM_C__

#include "eqm.h"

// wage (NC), sectoral capital (NC*NS), sectoral labor (NC*NS), gross output prices (NC*NS),
// sectoral consumption (NC*NS)
const uint nbgp = NC+3*NC*NS+NC*(NS-1);

void set_neqm()
{
  // key eqm vars used in solver: wage (NC), sectoral investment (NC*NS, but not in last period),
  // sectoral labor (NC*NS), gross output prices (NC*NS), sectoral consumption (NC*2), 
  // bonds (NC-1, but not in first period), bond price (1, but not in last period)

  // counterfactual or fixed-wedge equilibria
  neqm = (NT+1)*NC + 2*(NT+1)*NC*NS + (NT+1)*NC*(NS-1) + NT*NC*NS + NT*(NC-1) + NT;
  if(f_adj_cost)
    {
      neqm = neqm + NT*NC*(NS-1)*NC + NT*NC*(NS-1);
    }
  if(m_adj_cost)
    {
      neqm = neqm + NT*NC*(NS-1)*NC + NT*NC*(NS-1);
    }
  if(scenario>=1)
    {
      neqm = neqm + NMATCH*2*NC + (NMATCH-2)*2*NC;
    }
}

void init_vars(eqm * e)
{
  SET_ALL_V(e->b_t,(NT+1)*NC,0.0);
  SET_ALL_V(e->cc_t,(NT+1)*NC,0.0);
  SET_ALL_V(e->ii_t,(NT+1)*NC,0.0);
  SET_ALL_V(e->ll_t,(NT+1)*NC,0.0);
  SET_ALL_V(e->kk_t,(NT+1)*NC,0.0);
  SET_ALL_V(e->cpi_t,(NT+1)*NC,0.0);
  SET_ALL_V(e->pi_t,(NT+1)*NC,0.0);
  SET_ALL_V(e->w_t,(NT+1)*NC,0.0);
  SET_ALL_V(e->rk_t,(NT+1)*NC,0.0);
  SET_ALL_V(e->ngdp_t,(NT+1)*NC,0.0);
  SET_ALL_V(e->rgdp_t,(NT+1)*NC,0.0);
  SET_ALL_V(e->iy_t,(NT+1)*NC,0.0);
  SET_ALL_V(e->cy_t,(NT+1)*NC,0.0);
  SET_ALL_V(e->tby_t,(NT+1)*NC,0.0);
  SET_ALL_V(e->reer_t,(NT+1)*NC,0.0);
  SET_ALL_V(e->lp_agg_t,(NT+1)*NC,0.0);
  SET_ALL_V(e->tau_a_ts,(NT+1)*NC,1.0);
  SET_ALL_V(e->tau_i_ts,(NT+1)*NC,1.0);
  SET_ALL_V(e->tau_s_ts,(NT+1)*NC,1.0);
  SET_ALL_V(e->tau_t_ts,(NT+1)*NC,1.0);

  SET_ALL_V(e->nx_t,(NT+1)*NC*NC,0.0);
  SET_ALL_V(e->ex_t,(NT+1)*NC*NC,0.0);
  SET_ALL_V(e->im_t,(NT+1)*NC*NC,0.0);
  SET_ALL_V(e->rer_t,(NT+1)*NC*NC,0.0);

  SET_ALL_V(e->y_t,(NT+1)*NC*NS,0.0);
  SET_ALL_V(e->py_t,(NT+1)*NC*NS,0.0);
  SET_ALL_V(e->va_t,(NT+1)*NC*NS,0.0);
  SET_ALL_V(e->md_t,(NT+1)*NC*NS*(NS-1),0.0);
  SET_ALL_V(e->k_t,(NT+1)*NC*NS,0.0);
  SET_ALL_V(e->l_t,(NT+1)*NC*NS,0.0);
  SET_ALL_V(e->is_t,(NT+1)*NC*NS,0.0);
  SET_ALL_V(e->lp_t,(NT+1)*NC*NS,0.0);

  SET_ALL_V(e->nxsf_t,(NT+1)*NC*(NS-1)*NC,0.0);
  SET_ALL_V(e->nxsm_t,(NT+1)*NC*(NS-1)*NC,0.0);

  SET_ALL_V(e->pm_t,(NT+1)*NC*(NS-1),0.0);
  SET_ALL_V(e->m_t,(NT+1)*NC*(NS-1),0.0);
  SET_ALL_V(e->m2_t,(NT+1)*NC*(NS-1)*NC,0.0);

  SET_ALL_V(e->p_t,(NT+1)*NC*(NS-1),0.0);
  SET_ALL_V(e->q_t,(NT+1)*NC*(NS-1),0.0);
  SET_ALL_V(e->q2_t,(NT+1)*NC*(NS-1)*NC,0.0);

  SET_ALL_V(e->c_t,(NT+1)*NC*(NS-1),0.0);
  SET_ALL_V(e->i_t,(NT+1)*NC*NS,0.0);
}

void copy_vars(eqm * e1, const eqm * e0)
{
  memcpy((double *)(e1->pb_t),(const double *)(e0->pb_t),sizeof(double)*(NT+1));

  memcpy((double *)(e1->b_t),(const double *)(e0->b_t),sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->cc_t),(const double *)(e0->cc_t),sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->ii_t),(const double *)(e0->ii_t),sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->ll_t),(const double *)(e0->ll_t),sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->kk_t),(const double *)(e0->kk_t),sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->cpi_t),
	 (const double *)(e0->cpi_t),
	 sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->pi_t),(const double *)(e0->pi_t),sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->w_t),(const double *)(e0->w_t),sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->rk_t),(const double *)(e0->rk_t),sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->ngdp_t),
	 (const double *)(e0->ngdp_t),
	 sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->rgdp_t),
	 (const double *)(e0->rgdp_t),
	 sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->iy_t),(const double *)(e0->iy_t),sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->cy_t),(const double *)(e0->cy_t),sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->tby_t),(const double *)(e0->tby_t),sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->reer_t),(const double *)(e0->reer_t),sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->tau_a_ts),(const double *)(e0->tau_a_ts),sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->tau_i_ts),(const double *)(e0->tau_i_ts),sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->tau_s_ts),(const double *)(e0->tau_s_ts),sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->lp_agg_t),
	 (const double *)(e0->lp_agg_t),
	 sizeof(double)*(NT+1)*NC);

  memcpy((double *)(e1->nx_t),
	 (const double *)(e0->nx_t),
	 sizeof(double)*(NT+1)*NC*NC);
  memcpy((double *)(e1->rer_t),
	 (const double *)(e0->rer_t),
	 sizeof(double)*(NT+1)*NC*NC);

  memcpy((double *)(e1->y_t),
	 (const double *)(e0->y_t),
	 sizeof(double)*(NT+1)*NC*NS);
  memcpy((double *)(e1->py_t),
	 (const double *)(e0->py_t),
	 sizeof(double)*(NT+1)*NC*NS);
  memcpy((double *)(e1->va_t),
	 (const double *)(e0->va_t),
	 sizeof(double)*(NT+1)*NC*NS);
  memcpy((double *)(e1->k_t),
	 (const double *)(e0->k_t),
	 sizeof(double)*(NT+1)*NC*NS);
  memcpy((double *)(e1->l_t),
	 (const double *)(e0->l_t),
	 sizeof(double)*(NT+1)*NC*NS);
  memcpy((double *)(e1->is_t),
	 (const double *)(e0->is_t),
	 sizeof(double)*(NT+1)*NC*NS);
  memcpy((double *)(e1->md_t),
	 (const double *)(e0->md_t),
	 sizeof(double)*(NT+1)*NC*NS*(NS-1));
  memcpy((double *)(e1->lp_t),
	 (const double *)(e0->lp_t),
	 sizeof(double)*(NT+1)*NC*NS);

  memcpy((double *)(e1->nxsf_t),
	 (const double *)(e0->nxsf_t),
	 sizeof(double)*(NT+1)*NC*(NS-1)*NC);
  memcpy((double *)(e1->nxsm_t),
	 (const double *)(e0->nxsm_t),
	 sizeof(double)*(NT+1)*NC*(NS-1)*NC);

  memcpy((double *)(e1->pm_t),
	 (const double *)(e0->pm_t),
	 sizeof(double)*(NT+1)*NC*(NS-1));
  memcpy((double *)(e1->m_t),
	 (const double *)(e0->m_t),
	 sizeof(double)*(NT+1)*NC*(NS-1));
  memcpy((double *)(e1->m2_t),
	 (const double *)(e0->m2_t),
	 sizeof(double)*(NT+1)*NC*(NS-1)*NC);

  memcpy((double *)(e1->p_t),
	 (const double *)(e0->p_t),
	 sizeof(double)*(NT+1)*NC*(NS-1));
  memcpy((double *)(e1->q_t),
	 (const double *)(e0->q_t),
	 sizeof(double)*(NT+1)*NC*(NS-1));
  memcpy((double *)(e1->q2_t),
	 (const double *)(e0->q2_t),
	 sizeof(double)*(NT+1)*NC*(NS-1)*NC);

  memcpy((double *)(e1->c_t),
	 (const double *)(e0->c_t),
	 sizeof(double)*(NT+1)*NC*(NS-1));
  memcpy((double *)(e1->i_t),
	 (const double *)(e0->i_t),
	 sizeof(double)*(NT+1)*NC*NS);
}

uint stack_bgp_vars(double * myx, const eqm * e)
{
  uint nx = 0;
  uint t = NT;
  
  COPY_SUBVECTOR_LOG(myx+nx,e->w_t[t],NC);
  nx=nx+NC;

  COPY_SUBVECTOR_LOG(myx+nx,e->k_t[t],NC*NS);
  nx=nx+NC*NS;

  COPY_SUBVECTOR_LOG(myx+nx,e->l_t[t],NC*NS);
  nx=nx+NC*NS;

  COPY_SUBVECTOR_LOG(myx+nx,e->py_t[t],NC*NS);
  nx=nx+NC*NS;

  COPY_SUBVECTOR_LOG(myx+nx,e->c_t[t],NC*(NS-1));
  nx=nx+NC*(NS-1);

  if(nx != nbgp)
    {
      fprintf(logfile,KRED "Error stacking bgp vars! nx = %d, nbgp = %d\n" RESET,nx,nbgp);
      return 1;
    }

    return 0;
}

uint unstack_bgp_vars(eqm * e, const double * myx)
{
  uint nx = 0;
  uint t = NT;

  copy_subvector_exp( (double *)(e->w_t[t]), myx+nx, NC);
  nx=nx+NC;

  COPY_SUBVECTOR_EXP(e->k_t[t],myx+nx,NC*NS);
  nx=nx+NC*NS;

  COPY_SUBVECTOR_EXP(e->l_t[t],myx+nx,NC*NS);
  nx=nx+NC*NS;

  COPY_SUBVECTOR_EXP(e->py_t[t],myx+nx,NC*NS);
  nx=nx+NC*NS;

  COPY_SUBVECTOR_EXP(e->c_t[t],myx+nx,NC*(NS-1));
  nx=nx+NC*(NS-1);

  if(nx != nbgp)
    {
      fprintf(logfile,KRED "Error stacking bgp vars! nx = %d, nbgp = %d\n" RESET,nx,nbgp);
      return 1;
    }

    return 0;
}

uint stack_eqm_vars(double * myx, const eqm * e)
{
  uint nx = 0;
  uint t0 = 0;
  uint nn = NT+1;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->w_t[t0]),(nn)*NC);
  nx = nx + (nn)*NC;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->is_t[t0]),(nn-1)*NC*NS);
  nx = nx + (nn-1)*NC*NS;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->l_t[t0]),(nn)*NC*NS);
  nx = nx + (nn)*NC*NS;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->py_t[t0]),(nn)*NC*NS);
  nx = nx + (nn)*NC*NS;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->c_t[t0]),(nn)*NC*(NS-1));
  nx = nx + (nn)*NC*(NS-1);

  uint i = 0;
  uint t = 0;
  for(t=t0+1; t<(NT+1); t++)
    {
      for(i=0; i<(NC-1); i++)
	{
	  *(myx+nx) = e->b_t[t][i];
	  nx = nx+1;
	}
    }

  COPY_SUBVECTOR_LOG(myx+nx,&(e->pb_t[t0]),(nn-1));
  nx = nx + (nn-1);

  if(f_adj_cost)
    {
      COPY_SUBVECTOR_LOG(myx+nx,&(e->q2_t[t0]),(nn-1)*NC*(NS-1)*NC);
      nx = nx + (nn-1)*NC*(NS-1)*NC;

      COPY_SUBVECTOR_LOG(myx+nx,&(e->p_t[t0]),(nn-1)*NC*(NS-1));
      nx = nx + (nn-1)*NC*(NS-1);
    }

  if(m_adj_cost)
    {
      if(noio_flag==0)
	{
	  COPY_SUBVECTOR_LOG(myx+nx,&(e->m2_t[t0]),(nn-1)*NC*(NS-1)*NC);
	  nx = nx + (nn-1)*NC*(NS-1)*NC;

	  COPY_SUBVECTOR_LOG(myx+nx,&(e->pm_t[t0]),(nn-1)*NC*(NS-1));
	  nx = nx + (nn-1)*NC*(NS-1);
	}
      else
	{
	  COPY_SUBVECTOR(myx+nx,&(e->m2_t[t0]),(nn-1)*NC*(NS-1)*NC);
	  nx = nx + (nn-1)*NC*(NS-1)*NC;

	  COPY_SUBVECTOR(myx+nx,&(e->pm_t[t0]),(nn-1)*NC*(NS-1));
	  nx = nx + (nn-1)*NC*(NS-1);
	}
    }

  if(scenario>=1)
    {
      for(t=t0; t<NMATCH; t++)
	{
	  for(i=0; i<NC; i++)
	    {
	      *(myx+nx) = log(e->tau_s_ts[t][i]);
	      nx=nx+1;

	      *(myx+nx) = log(e->tau_i_ts[t][i]);
	      nx=nx+1;

	      if(t<NMATCH-2)
		{
		  *(myx+nx) = log(e->tau_a_ts[t+1][i]);
		  nx=nx+1;

		  *(myx+nx) = log(e->tau_t_ts[t+1][i]);
		  nx=nx+1;
		}
	    }
	}
    }

  if(nx != neqm)
    {
      fprintf(logfile,KRED "Error stacking eqm vars! nx = %d, neqm = %d\n" RESET,nx,neqm);
      return 1;
    }

    return 0;
}

uint unstack_eqm_vars(eqm * e, const double * myx)
{
  uint nx = 0;
  uint t0 = 0;
  uint nn = NT+1;

  COPY_SUBVECTOR_EXP(&(e->w_t[t0]),myx+nx,(nn)*NC);
  nx = nx + (nn)*NC;
  
  COPY_SUBVECTOR_EXP(&(e->is_t[t0]),myx+nx,(nn-1)*NC*NS);
  nx = nx + (nn-1)*NC*NS;

  COPY_SUBVECTOR_EXP(&(e->l_t[t0]),myx+nx,(nn)*NC*NS);
  nx = nx + (nn)*NC*NS;

  COPY_SUBVECTOR_EXP(&(e->py_t[t0]),myx+nx,(nn)*NC*NS);
  nx = nx + (nn)*NC*NS;

  COPY_SUBVECTOR_EXP(&(e->c_t[t0]),myx+nx,(nn)*NC*(NS-1));
  nx = nx + (nn)*NC*(NS-1);

  uint i = 0;
  uint t = 0;
  for(t=t0+1; t<(NT+1); t++)
    {
      for(i=0; i<(NC-1); i++)
	{
	  e->b_t[t][i] = *(myx+nx);
	  nx = nx+1;
	}
    }

  COPY_SUBVECTOR_EXP(&(e->pb_t[t0]),myx+nx,(nn-1));
  nx = nx + (nn-1);

  if(f_adj_cost)
    {
      COPY_SUBVECTOR_EXP(&(e->q2_t[t0]),myx+nx,(nn-1)*NC*(NS-1)*NC);
      nx = nx + (nn-1)*NC*(NS-1)*NC;

      COPY_SUBVECTOR_EXP(&(e->p_t[t0]),myx+nx,(nn-1)*NC*(NS-1));
      nx = nx + (nn-1)*NC*(NS-1);
    }

  if(m_adj_cost)
    {
      if(noio_flag==0)
	{
	  COPY_SUBVECTOR_EXP(&(e->m2_t[t0]),myx+nx,(nn-1)*NC*(NS-1)*NC);
	  nx = nx + (nn-1)*NC*(NS-1)*NC;
	  
	  COPY_SUBVECTOR_EXP(&(e->pm_t[t0]),myx+nx,(nn-1)*NC*(NS-1));
	  nx = nx + (nn-1)*NC*(NS-1);
	}
      else
	{
	  COPY_SUBVECTOR(&(e->m2_t[t0]),myx+nx,(nn-1)*NC*(NS-1)*NC);
	  nx = nx + (nn-1)*NC*(NS-1)*NC;
	  
	  COPY_SUBVECTOR(&(e->pm_t[t0]),myx+nx,(nn-1)*NC*(NS-1));
	  nx = nx + (nn-1)*NC*(NS-1);

	}
    }

  if(scenario>=1)
    {
      for(t=t0; t<NMATCH; t++)
	{
	  for(i=0; i<NC; i++)
	    {
	      e->tau_s_ts[t][i] = exp(*(myx+nx));
	      nx=nx+1;

	      e->tau_i_ts[t][i] = exp(*(myx+nx));
	      nx=nx+1;

	      if(t<NMATCH-2)
		{
		  e->tau_a_ts[t+1][i] = exp(*(myx+nx));
		  nx=nx+1;

		  e->tau_t_ts[t+1][i] = exp(*(myx+nx));
		  nx=nx+1;
		}
	    }
	}
    }

  if(nx != neqm)
    {
      fprintf(logfile,KRED "Error unstacking eqm vars! nx = %d, neqm = %d\n" RESET ,nx,neqm);
      return 1;
    }

  return 0;
}

uint set_initial_bgp_guess()
{
  uint i, t, s;
  eqm * e = &(eee0[0]);
  params * p = &(ppp0[0]);
  t=NT;
  
  for(i=0; i<NC; i++)
    {
      e->w_t[t][i] = 1.0;
      for(s=0; s<NS; s++)
	{
	  if(s<NS-1)
	    {
	      e->c_t[t][i][s] = p->c0[i][s] * p->a_ts[t][i] * p->lbar_ts[t][i]/p->lbar0[i];
	    }
	  e->l_t[t][i][s] = p->l0[i][s] * p->lbar_ts[t][i]/p->lbar0[i];
	  e->k_t[t][i][s] = p->k0[i][s] * p->a_ts[t][i] * p->lbar_ts[t][i]/p->lbar0[i];
	  e->py_t[t][i][s] = 1.0;
	}
    }

  if(stack_bgp_vars(solver_x->data,e))
    {
      fprintf(logfile,KRED "Failed to create guess for balanced growth path!\n" RESET);
      return 1;
    }
  else
    {
      return 0;
    }
}

uint set_initial_eqm_guess()
{
  uint i,s,t,j;
  double bb[NC];

  eqm * e = &(eee0[0]);
  params * p = &(ppp0[0]);

  bb[0] = p->b0[0] * pow(p->ga_bgp,NT);
  bb[1] = p->b0[1] * pow(p->ga_bgp,NT);
  if(NC==3)
    {
      bb[2] = p->b0[2] * pow(p->ga_bgp,NT);
    }

  free_solver_mem();
  solver_n = nbgp;
  alloc_solver_mem();
  if(solve_bgp(bb))
    {
      fprintf(logfile, KRED "Error solving for balanced growth path!\n");
      return 1;
    }
  free_solver_mem();
  solver_n = neqm;
  alloc_solver_mem();

  // first construct bond guess... a little awkward to logspace this because we have to deal with
  // absolute values
  double tmpb0[NC]; //= {fabs(p->b0[0]),fabs(p->b0[1]),fabs(p->b0[2]),fabs(p->b0[3])};
  double tmpb1[NC]; //= {fabs(bb[0]),fabs(bb[1]),fabs(bb[2]),fabs(bb[3])};
  for(i=0; i<NC; i++)
    {
      tmpb0[i] = fabs(p->b0[i]);
      tmpb1[i] = fabs(bb[i]);
    }

  LOGSPACE_2D(tmpb0,tmpb1,NT+1,NC,e->b_t);
  for(i=0; i<(NC-1); i++)
    {
      if(fabs(p->b0[i])<1.0e-6)
	{
	  for(t=0; t<(NT+1); t++)
	    {
	      e->b_t[t][i] = 0.0;
	    }
	}
      else
	{
	  if(p->b0[i] < -TINY)
	    {
	      for(t=0; t<(NT+1); t++)
		{
		  e->b_t[t][i] = -e->b_t[t][i];
		}
	    }
	}
    }
  set_all_v(e->pb_t,NT+1,e->pb_t[NT]);

  // now construct guesses for prices real variables
  double tmpp[NC];
  double tmpp2[NC][NS];
  double tmpp3[NC][NS-1];
  for(i=0; i<NC; i++)
    {
      tmpp[i]=1.0;
      for(s=0; s<NS; s++)
	{
	  tmpp2[i][s] = 1.0;
	  if(s<NS-1)
	    {
	      tmpp3[i][s] = 1.0;
	    }
	}
    }

  LOGSPACE_2D(p->k0,e->k_t[NT],NT+1,NC*NS,e->k_t);
  LOGSPACE_2D(p->l0,e->l_t[NT],NT+1,NC*NS,e->l_t);
  LOGSPACE_2D(p->c0,e->c_t[NT],NT+1,NC*(NS-1),e->c_t);
  LOGSPACE_2D(p->q0,e->q_t[NT],NT+1,NC*(NS-1),e->q_t);
  LOGSPACE_2D(p->m0,e->m_t[NT],NT+1,NC*(NS-1),e->m_t);
  LINSPACE_2D(tmpp,e->w_t[NT],NT+1,NC,e->w_t);
  LINSPACE_2D(tmpp2,e->py_t[NT],NT+1,NC*NS,e->py_t);
  LINSPACE_2D(tmpp3,e->p_t[NT],NT+1,NC*(NS-1),e->p_t);
  LINSPACE_2D(tmpp3,e->pm_t[NT],NT+1,NC*(NS-1),e->pm_t);

  for(i=0; i<NC; i++)
    {
      for(s=0; s<NS; s++)
	{
	  for(t=0; t<NT; t++)
	    {
	      e->is_t[t][i][s] = e->k_t[t+1][i][s] - (1.0-p->delta)*e->k_t[t][i][s];
	    }
	}
    }

  for(t=0; t<(NT+1); t++)
    {
      for(i=0; i<NC; i++)
	{
	  for(s=0; s<NS-1; s++)
	    {
	      for(j=0; j<NC; j++)
		{
		  e->q2_t[t][i][s][j] = (p->q02[i][s][j]/p->q0[i][s]) * e->q_t[t][i][s];
		  if(noio_flag)
		    {
		      e->m2_t[t][i][s][j] = 0.0;
		    }
		  else
		    {
		      e->m2_t[t][i][s][j] = (p->m02[i][s][j]/p->m0[i][s]) * e->m_t[t][i][s];
		    }
		}
	    }
      
	}
    }

  if(stack_eqm_vars(solver_x->data,e))
    {
      fprintf(logfile,KRED "Failed to create guess for balanced growth path!\n" RESET);
      return 1;
    }
  else
    {      
      return 0;
    }
}

uint write_eqm_vars(const eqm * e, char * fname, uint i)
{
  char * fname2 = concat("output/",fname);

  char * fname3;
  if(!tau_t_flag && !tau_s_only_flag && !tau_s_only_flag2 && !noio_flag)
    {
      fname3 = concat(fname2,"_notrdw.csv");
    }
  else if(no_demo_flag)
    {
      fname3 = concat(fname2,"_nodemo.csv");
    }
  else if(sym_g_flag)
    {
      fname3 = concat(fname2,"_symg.csv");
    }
  else if(no_rw_tau_i_flag)
    {
      fname3 = concat(fname2,"_notaui.csv");
    }
  else if(no_rw_tau_i_flag2)
    {
      fname3 = concat(fname2,"_notaui2.csv");
    }
  else if(tau_s_only_flag)
    {
      fname3 = concat(fname2,"_taus.csv");
    }
  else if(tau_s_only_flag2)
    {
      fname3 = concat(fname2,"_taus2.csv");
    }
  else if(noio_flag==1)
    {
      fname3 = concat(fname2,"_noio.csv");
    }
  else if(high_etaK_flag==1)
    {
      fname3 = concat(fname2,"_etaKh.csv");
    }
  else if(low_etaK_flag==1)
    {
      fname3 = concat(fname2,"_etaKl.csv");
    }
  else if(high_arm_flag==1)
    {
      fname3 = concat(fname2,"_armh.csv");
    }
  else if(low_arm_flag==1)
    {
      fname3 = concat(fname2,"_arml.csv");
    }
  else if(high_rhow_flag==1)
    {
      fname3 = concat(fname2,"_rhowh.csv");
    }
  else if(low_rhow_flag==1)
    {
      fname3 = concat(fname2,"_rhowl.csv");
    }
  else if(low_r_flag==1)
    {
      fname3 = concat(fname2,"_lowr.csv");
    }
  else
    {
      fname3 = concat(fname2,".csv");
    }

  FILE * file = fopen(fname3,"wb");
  free(fname2);

  if(file)
    {
      uint s,j,t;
      fprintf(file,"period,rgdp,ngdp,iy,cy,tby,nfay,reer,rir,pi,tau_a,tau_s,tau_i,tau_t");
      for(s=0;s<NS;s++)
	{
	  fprintf(file,",y%d,va%d,i%d,l%d",s,s,s,s);
	}
      for(j=0; j<NC; j++)
	{
	  if(j!=i)
	    {
	      //fprintf(file,",rer%d,tb%d",j,j);
	      for(s=0;s<NS-1;s++)
		{
		  fprintf(file,",tbsm%d-%d,tbsf%d-%d",s,j,s,j);
		}
	    }
	}
      fprintf(file,"\n");

      for(t=0;t<(NT+1);t++)
	{
	  fprintf(file,"%d,",t);
	  fprintf(file,"%0.16f,",e->rgdp_t[t][i]);
	  fprintf(file,"%0.16f,",e->ngdp_t[t][i]);
	  fprintf(file,"%0.16f,",e->iy_t[t][i]);
	  fprintf(file,"%0.16f,",e->cy_t[t][i]);
	  fprintf(file,"%0.16f,",e->tby_t[t][i]);
	  fprintf(file,"%0.16f,",e->nfay_t[t][i]);
	  if(i==0 && tau_t_flag && scenario==1 && (t==16 || t==17))
	    {
	      fprintf(file,"%0.16f,",ppp0[0].RER_ts[t]);
	    }
	  else
	    {
	      fprintf(file,"%0.16f,",e->reer_t[t][i]);
	    }
	  fprintf(file,"%0.16f,",e->rir_t[t][i]);
	  fprintf(file,"%0.16f,",100.0*e->p_t[t][i][1]);
	  fprintf(file,"%0.16f,",e->tau_a_ts[t][i]);
	  fprintf(file,"%0.16f,",e->tau_s_ts[t][i]);
	  fprintf(file,"%0.16f,",e->tau_i_ts[t][i]);
	  fprintf(file,"%0.16f",e->tau_t_ts[t][i]);
	  for(s=0; s<NS; s++)
	    {
	      fprintf(file,",%0.16f",e->y_t[t][i][s]);
	      fprintf(file,",%0.16f",e->va_t[t][i][s]);
	      fprintf(file,",%0.16f",e->is_t[t][i][s]);
	      fprintf(file,",%0.16f",e->l_t[t][i][s]);
	    }
	  for(j=0; j<NC; j++)
	    {
	      if(j!=i)
		{
		  //fprintf(file,",%0.16f",e->rer_t[t][i][j]);
		  //fprintf(file,",%0.16f",100.0*e->nx_t[t][i][j]/e->ngdp_t[t][i]);
		  for(s=0; s<NS-1; s++)
		    {
		      fprintf(file,",%0.16f",100.0*e->nxsm_t[t][i][s][j]/e->ngdp_t[t][i]);
		      fprintf(file,",%0.16f",100.0*e->nxsf_t[t][i][s][j]/e->ngdp_t[t][i]);
		    }
		}
	    }
	  fprintf(file,"\n");
	}
      
      fclose(file);
      return 0;
    }
  else
    {
      fprintf(logfile,KRED "Error opening file to write equilibrium vars!\n" RESET);
      return 1;
    }
}

void set_vars(eqm * e, const params * p, uint t, uint bgp)
{
  uint i,s,r,j;

  e->b_t[t][NC-1] = -sum(e->b_t[t],NC-1);

  SET_ALL_V(e->ngdp_t[t],NC,0.0);
  SET_ALL_V(e->rgdp_t[t],NC,0.0);
  SET_ALL_V(e->nx_t[t],NC*NC,0.0);
  SET_ALL_V(e->ex_t[t],NC*NC,0.0);
  SET_ALL_V(e->im_t[t],NC*NC,0.0);
  SET_ALL_V(e->nxsm_t[t],NC*(NS-1)*NC,0.0);
  SET_ALL_V(e->nxsf_t[t],NC*(NS-1)*NC,0.0);

  for(i=0; i<NC; i++)
    {
      for(s=0; s<NS; s++)
	{
  
	  e->va_t[t][i][s] = e->tau_a_ts[t][i] * 
	    prod_va(e->k_t[t][i][s],e->l_t[t][i][s],p->A[i][s],p->alpha[i][s],p->a_ts[t][i]);
	  e->y_t[t][i][s] = e->va_t[t][i][s]/p->lam_va[i][s];
	  e->ngdp_t[t][i] = e->ngdp_t[t][i] + e->py_t[t][i][s]*e->y_t[t][i][s];
	  e->rgdp_t[t][i] = e->rgdp_t[t][i] + e->y_t[t][i][s];
	  e->lp_t[t][i][s] = e->va_t[t][i][s]/e->l_t[t][i][s];

	  for(r=0; r<NS-1; r++)
	    {
	      e->md_t[t][i][s][r] = e->y_t[t][i][s]*p->lam[i][s][r];
	    }
	}

      e->ll_t[t][i] = e->l_t[t][i][0] + e->l_t[t][i][1] + e->l_t[t][i][2];

      // adjustment costs are in units of labor, so we need to incremement labor demand by them
      // ??? do we need to do something to GDP as well with this?
      if(f_adj_cost && t<NT)
	{
	  for(s=0; s<NS-1; s++)
	    {
	      for(j=0; j<NC; j++)
		{
		  if(j!=i)
		    {
		      double tmp = 0.0;
		      if(t>0)
			{
			  tmp = e->q2_t[t][i][s][j]/e->q2_t[t-1][i][s][j] - 1.0;
			}
		      else
			{
			  tmp = e->q2_t[t][i][s][j]/p->q02[i][s][j] - 1.0;
			}
		      e->ll_t[t][i] = e->ll_t[t][i] + 
			(p->etaF/2.0) * tmp * tmp;
		    }
		}
	    }
	}

      if(m_adj_cost && t<NT && noio_flag==0)
	{
	  for(s=0; s<NS-1; s++)
	    {
	      for(j=0; j<NC; j++)
		{
		  if(j!=i)
		    {
		      double tmp = 0.0;
		      if(t>0)
			{
			  tmp = e->m2_t[t][i][s][j]/e->m2_t[t-1][i][s][j] - 1.0;
			}
		      else
			{
			  tmp = e->m2_t[t][i][s][j]/p->m02[i][s][j] - 1.0;
			}
		      e->ll_t[t][i] = e->ll_t[t][i] + 
			(p->etaM/2.0) * tmp * tmp;
		    }
		}
	    }
	}

      if(t<NT)
	{
	  if(t==(NT-1) || k_adj_cost==0)
	    {
	      for(s=0; s<NS; s++)
		{
		  e->k_t[t+1][i][s] = (1.0-p->delta) * e->k_t[t][i][s] + e->is_t[t][i][s];
		}
	    }
	  else
	    {
	      for(s=0; s<NS; s++)
		{
		  e->k_t[t+1][i][s] = (1.0-p->delta) * e->k_t[t][i][s] + 
		    phiK(e->is_t[t][i][s]/e->k_t[t][i][s],p->delta,p->ga_bgp,p->etaK) * e->k_t[t][i][s];
		}
	    }
	}
      else
	{
	  for(s=0; s<NS; s++)
	    {
	      e->is_t[t][i][s] = (p->ga_bgp-1.0+p->delta) * e->k_t[t][i][s];
	    }
	}
      
      e->kk_t[t][i] = 0.0;
      e->ii_t[t][i] = 0.0;
      for(s=0; s<NS; s++)
	{
	  e->kk_t[t][i] = e->kk_t[t][i] + e->k_t[t][i][s];
	  e->ii_t[t][i] = e->ii_t[t][i] + e->is_t[t][i][s];
	}

      if(!m_adj_cost || t==NT)
	{
	  for(s=0; s<NS-1; s++)
	    {
	      e->pm_t[t][i][s] = 0.0;
	      if(noio_flag==0)
		{
		  for(j=0; j<NC; j++)
		    {
		      double tmp = 1.0;
		      if(j!=i)
			{
			  tmp = e->tau_t_ts[t][i];
			}
		      e->pm_t[t][i][s] = e->pm_t[t][i][s] + 
			pow(p->mu[i][s][j],1.0/(1.0-p->zeta[i][s])) * 
			pow(tmp*e->py_t[t][j][s],p->zeta[i][s]/(p->zeta[i][s]-1.0));		      
		    }
		  e->pm_t[t][i][s] = (1.0/p->M[i][s]) * 
		    pow(e->pm_t[t][i][s],(p->zeta[i][s]-1.0)/p->zeta[i][s]);
		}
	    }
	}

      if(!f_adj_cost || t==NT)
	{
	  for(s=0; s<NS-1; s++)
	    {
	      e->p_t[t][i][s] = 0.0;
	      for(j=0; j<NC; j++)
		{
		  double tmp = 1.0;
		  if(j!=i)
		    {
		      tmp = e->tau_t_ts[t][i];
		    }
		  e->p_t[t][i][s] = e->p_t[t][i][s] + 
		    pow(p->theta[i][s][j],1.0/(1.0-p->sig[i][s])) * 
		    pow(tmp*e->py_t[t][j][s],p->sig[i][s]/(p->sig[i][s]-1.0));
		}
	      e->p_t[t][i][s] = (1.0/p->H[i][s]) * 
		pow(e->p_t[t][i][s],(p->sig[i][s]-1.0)/p->sig[i][s]);
	    }
	}

      // households' stochastic discount factor for dynamic firm's problem with adjustment costs
      if(t>0)
	{
	  double mutp = muc(
		e->c_t[t][i],
		e->ll_t[t][i],
		p->lbar_ts[t][i],
		p->pope_ts[t][i],
		p->popw_ts[t][i],
		p->eps[i][0],
		p->rho,
		p->phi[i],
		p->psi,
		0);
	  double mut = muc(
		  e->c_t[t-1][i],
		  e->ll_t[t-1][i],
		  p->lbar_ts[t][i],
		  p->pope_ts[t][i],
		  p->popw_ts[t][i],
		  p->eps[i][0],
		  p->rho,
		  p->phi[i],
		  p->psi,
		  0);

	  e->Q_t[t-1][i] = p->beta[i] * (mutp / e->p_t[t][i][0]) / (mut / e->p_t[t-1][i][0]);
	}
      if(t==NT)
	{
	  e->Q_t[t][i] = p->beta[i];
	}

      e->pi_t[t][i] = 1.0/p->G[i];
      for(s=0;s<NS-1; s++)
	{
	  e->pi_t[t][i] = e->pi_t[t][i] * pow(e->p_t[t][i][s]/p->eps[i][1][s],p->eps[i][1][s]);
	}
      e->pi_t[t][i] = e->pi_t[t][i] * pow(e->py_t[t][i][NS-1]/p->eps[i][1][NS-1],p->eps[i][1][NS-1]);
      for(s=0; s<NS-1; s++)
	{
	  e->i_t[t][i][s] = e->pi_t[t][i] * p->eps[i][1][s] * e->ii_t[t][i]/e->p_t[t][i][s];
	}
      e->i_t[t][i][NS-1] = e->pi_t[t][i] * p->eps[i][1][NS-1] * e->ii_t[t][i]/e->py_t[t][i][NS-1];

      for(s=0; s<NS-1; s++)
	{
	  e->q_t[t][i][s] = e->c_t[t][i][s] + e->i_t[t][i][s];
	  e->m_t[t][i][s] = 0.0;
	  for(r=0; r<NS; r++)
	    {
	      e->m_t[t][i][s] = e->m_t[t][i][s] + e->md_t[t][i][r][s];
	    }

	  if(!m_adj_cost || t==NT)
	    {
	      for(j=0; j<NC; j++)
		{
		  double tmp = 1.0;
		  if(j!=i)
		    {
		      tmp = e->tau_t_ts[t][i];
		    }
		  if(noio_flag)
		    {
		      e->m2_t[t][i][s][j] = 0.0;
		    }
		  else
		    {
		      e->m2_t[t][i][s][j] = e->m_t[t][i][s] * 
		      pow(tmp*e->py_t[t][j][s],1.0/(p->zeta[i][s]-1.0)) * 
		      pow(e->pm_t[t][i][s]*p->mu[i][s][j]*pow(p->M[i][s],p->zeta[i][s]),1.0/(1.0-p->zeta[i][s]));
		    }
		}
	    }

	  if(!f_adj_cost || t==NT)
	    {
	      for(j=0; j<NC; j++)
		{
		  double tmp = 1.0;
		  if(j!=i)
		    {
		      tmp = e->tau_t_ts[t][i];
		    }
		  e->q2_t[t][i][s][j] = e->q_t[t][i][s] * 
		    pow(tmp*e->py_t[t][j][s],1.0/(p->sig[i][s]-1.0)) * 
		    pow(e->p_t[t][i][s]*p->theta[i][s][j]*pow(p->H[i][s],p->sig[i][s]),1.0/(1.0-p->sig[i][s]));
		}
	    }
	}

      e->cpi_t[t][i] = 0.0;
      e->cc_t[t][i] = 0.0;
      for(s=0; s<NS-1; s++)
	{
	  e->cpi_t[t][i] = e->cpi_t[t][i] + e->p_t[t][i][s]*p->c0[i][s];
	  e->cc_t[t][i] = e->cc_t[t][i]+e->c_t[t][i][s];
	}
      e->cpi_t[t][i] = e->cpi_t[t][i]/SUM(p->c0[i],NS-1);
      
      if(t == 0)
	{
	  e->rk_t[t][i] = p->r0[i] + p->delta;
	}
      else if(t==NT)
	{
	  if(i==0)
	    {
	      e->pb_t[t] = e->cpi_t[t][i]/(1.0+p->rss);
	    }
	  if(bgp)
	    {
	      e->rk_t[t][i] = e->pi_t[t][i]*e->cpi_t[t][0]/e->pb_t[t] - (1.0-p->delta)*e->pi_t[t][i];
	    }
	  else
	    {
	      e->rk_t[t][i] = e->pi_t[t-1][i]*e->cpi_t[t][0]/e->pb_t[t-1] - (1.0-p->delta)*e->pi_t[t][i];
	    }
	}
      else
	{
	  e->rk_t[t][i] = e->pi_t[t-1][i]*e->cpi_t[t][0]/e->pb_t[t-1] - (1.0-p->delta)*e->pi_t[t][i];
	}
    }

  for(i=0; i<NC; i++)
    {
      for(j=0; j<NC; j++)
	{
	  for(s=0; s<NS-1; s++)
	    {
	      e->ngdp_t[t][i] = e->ngdp_t[t][i] - e->py_t[t][j][s]*e->m2_t[t][i][s][j];
	      e->rgdp_t[t][i] = e->rgdp_t[t][i] - e->m2_t[t][i][s][j];
	    }

	  if(j!=i)
	    {
	      e->rer_t[t][i][j] = 100*e->cpi_t[t][j]/e->cpi_t[t][i];
	      for(s=0; s<NS-1; s++)
		{
		  double m = e->py_t[t][i][s]*e->m2_t[t][j][s][i];
		  double f = e->py_t[t][i][s]*e->q2_t[t][j][s][i];
	
		  e->nxsm_t[t][i][s][j] = m;
		  e->nxsf_t[t][i][s][j] = f;
		  e->ex_t[t][i][j] = e->ex_t[t][i][j] + m + f;

		  m = e->py_t[t][j][s]*e->m2_t[t][i][s][j];
		  f = e->py_t[t][j][s]*e->q2_t[t][i][s][j];

		  e->nxsm_t[t][i][s][j] = e->nxsm_t[t][i][s][j] - m;
		  e->nxsf_t[t][i][s][j] = e->nxsf_t[t][i][s][j] - f;
		  e->im_t[t][i][j] = e->im_t[t][i][j] + m + f;

		  e->nx_t[t][i][j] = e->nx_t[t][i][j] + e->nxsm_t[t][i][s][j] + e->nxsf_t[t][i][s][j];
		}
	    }
	  else
	    {
	      e->rer_t[t][i][j] = 1.0;
	    }
	}

      e->iy_t[t][i] = 100.0*e->pi_t[t][i]*e->ii_t[t][i]/e->ngdp_t[t][i];
      e->cy_t[t][i] = 100.0*(e->p_t[t][i][0]*e->c_t[t][i][0] + e->p_t[t][i][1]*e->c_t[t][i][1])/e->ngdp_t[t][i];
      e->tby_t[t][i] = 100.0*sum(e->nx_t[t][i],NC)/e->ngdp_t[t][i];
      e->nfay_t[t][i] = 100.0*e->cpi_t[t][0]*e->b_t[t][i]/e->ngdp_t[t][i];
      e->lp_agg_t[t][i] = e->rgdp_t[t][i] / e->ll_t[t][i];

      if(t<NT)
	{
	  //e->rir_t[t][i] = 100.0*(e->cpi_t[t+1][0] * e->cpi_t[t][i] / e->pb_t[t] / e->cpi_t[t+1][i] - 1.0);
	  e->rir_t[t][i] = 100.0*(1.0 / e->pb_t[t] - 1.0);
	}
      else
	{
	  e->rir_t[t][i] = p->rss;
	}
      if(NC==2)
	{
	  e->reer_t[t][i] = e->rer_t[t][i][1-i];
	}
      else
	{
	  uint j1, j2;
	  if(i==0)
	    {
	      j1=1;
	      j2=2;
	    }
	  else if(i==1)
	    {
	      j1=0;
	      j2=2;
	    }
	  else
	    {
	      j1=0;
	      j2=1;
	    }
	  e->reer_t[t][i] = ((e->ex_t[t][i][j1]+e->im_t[t][i][j1])*e->rer_t[t][i][j1]+
			     (e->ex_t[t][i][j2]+e->im_t[t][i][j2])*e->rer_t[t][i][j2])/
	    (e->ex_t[t][i][j1]+e->im_t[t][i][j1]+e->ex_t[t][i][j2]+e->im_t[t][i][j2]);
	}
    }  

  
  if(tau_t_flag && noio_flag==0)
    {
      if(t==NMATCH-1 && (scenario==1 || (scenario==2 && fix_tau_t[0]==1)))
	{
	  e->tau_t_ts[t][1] = 1.15;
	}
      if(t==NMATCH && (scenario==1 || (scenario==2 && fix_tau_t[0]==1)))
	{
	  e->tau_t_ts[t][1]  = 1.08;
	}
      if(t==NMATCH+1 && (scenario==1 || (scenario==2 && fix_tau_t[0]==1)))
	{
	  e->tau_t_ts[t][1]  = 1.04;
	}
    }
  
  if(t>=NMATCH && t<NT+1)
    {
      for(i=0; i<NC; i++)
	{
	  e->tau_s_ts[t][i] = p->wedge_speed * e->tau_s_ts[t-1][i] + (1.0-p->wedge_speed)*1.0;
	  e->tau_i_ts[t][i] = p->wedge_speed * e->tau_i_ts[t-1][i] + (1.0-p->wedge_speed)*1.0;
	  //	  e->tau_t_ts[t][i] = p->wedge_speed * e->tau_t_ts[t-1][i] + (1.0-p->wedge_speed)*1.0;
	}
    }
}

uint eval_bgp_conds(const double * myx, double * myf, uint tn)
{
  uint i=0,s=0,t=NT,nx=0;
  eqm * e = &(eee0[tn]);
  params * p = &(ppp0[tn]);

  for(i=0; i<NC-1; i++)
    {
      e->b_t[t][i] = bbgp[i];
    }
  unstack_bgp_vars(e,myx);
  set_vars(e,p,t,1);
  nx=0;

  myf[nx] = price_norm(e,t);
  nx=nx+1;

  for(i=0; i<(NC-1); i++)
    {
      myf[nx] = bop(p,e,t,i);
      nx = nx+1;
    }

  for(i=0; i<NC; i++)
    {
      for(s=0; s<NS; s++)
	{
	  myf[nx] = mpk_rk(p,e,t,i,s);
	  nx=nx+1;

	  myf[nx] = mpl_w(p,e,t,i,s);
	  nx=nx+1;

	  myf[nx] = mkt_clear_y(p,e,t,i,s);
	  nx=nx+1;
	}

      myf[nx] = mucs_mucr(p,e,t,i,0,1);
      nx=nx+1;

      myf[nx] = muc_mul(p,e,t,i);
      nx=nx+1;
    }

  if(nx != nbgp)
    {
      fprintf(logfile,KRED "Error evaluating bgp eqns! nx = %d, nbgp = %d\n" RESET,nx,nbgp);
      return 1;
    }
  
  return 0;
}

uint solve_bgp(double bb[NC-1])
{
  uint i;
  for(i=0; i<NC-1; i++)
    {
      bbgp[i] = bb[i];
    }
  solver_n = nbgp;
  alloc_solver_mem();
  set_initial_bgp_guess();
  gsl_multiroot_function_fdf f = {&bgp_func_f,&bgp_func_df,&bgp_func_fdf,nbgp,NULL};
  par=1;
  uint status = find_root_deriv_mkl(&f);

  free_solver_mem();
  return status;
}

uint eval_eqm_conds(const double * myx, double * myf, uint tn)
{
  eqm * e = &(eee0[tn]);
  params * p = &(ppp0[tn]);
  uint i=0,s=0,t=NT,nx=0;
  uint t0 = 0;

  if(scenario == 1)
    {
      e = &(eee1[tn]);
      p = &(ppp1[tn]);
    }
  else if(scenario == 2)
    {
      e = &(eee2[tn]);
      p = &(ppp2[tn]);
    }


  unstack_eqm_vars(e,myx);

  for(i=0; i<NC; i++)
    {
      e->b_t[0][i] = p->b0[i];
    }

  for(i=0; i<NC; i++)
    {
      for(s=0; s<NS; s++)
	{
	  e->k_t[0][i][s] = p->k0[i][s];
	}
    }

  for(t=t0; t<(NT+1); t++)
    {
      set_vars(e,p,t,0);
    }

  nx=0;
  for(t=t0; t<(NT+1); t++)
    {
      myf[nx] = price_norm(e,t);
      nx=nx+1;

      for(i=0; i<(NC-1); i++)
	{
	  myf[nx] = bop(p,e,t,i);
	  nx = nx+1;
	}

      for(i=0; i<NC; i++)
	{
	  for(s=0; s<NS; s++)
	    {
	      if(t<NT)
		{
		  myf[nx] = mpk_rk(p,e,t+1,i,s);
		  nx=nx+1;
		}

	      myf[nx] = mpl_w(p,e,t,i,s);
	      nx=nx+1;

	      myf[nx] = mkt_clear_y(p,e,t,i,s);
	      nx=nx+1;
	    }

	  myf[nx] = mucs_mucr(p,e,t,i,0,1);
	  nx=nx+1;

	  myf[nx] = muc_mul(p,e,t,i);
	  nx=nx+1;

	  if(t<NT)
	    {
	      myf[nx] = euler(p,e,t,i);
	      nx = nx+1;
	    }

	  if(m_adj_cost && t<NT)
	    {
	      for(s=0; s<NS-1; s++)
		{
		  myf[nx] = prod_m_chk(p,e,t,i,s);
		  nx = nx+1;

		  uint j;
		  for(j=0; j<NC; j++)
		    {
		      myf[nx] = foc_m2(p,e,t,i,s,j);
		      nx = nx+1;
		    }
		}
	    }

	  if(f_adj_cost && t<NT)
	    {
	      for(s=0; s<NS-1; s++)
		{
		  myf[nx] = prod_q_chk(p,e,t,i,s);
		  nx = nx+1;

		  uint j;
		  for(j=0; j<NC; j++)
		    {
		      myf[nx] = foc_q2(p,e,t,i,s,j);
		      nx = nx+1;
		    }
		}
	    }
	}
    }
  if(scenario==1)
    {
      for(t=0; t<NMATCH; t++)
	{
	  if(tau_s_only_flag || tau_s_only_flag2)
	    {
	      myf[nx] = e->tau_i_ts[t][0] - 1.0;
	      nx=nx+1;
	      myf[nx] = e->tau_i_ts[t][1] - 1.0;
	      nx = nx+1;
	      
	      if(t<NMATCH-2)
		{
		  myf[nx] = e->tau_t_ts[t+1][0] - 1.0;
		  nx=nx+1;
		  myf[nx] = e->tau_t_ts[t+1][1] - 1.0;
		  nx=nx+1;
		  myf[nx] = e->tau_a_ts[t+1][0] - 1.0;
		  nx=nx+1;
		  myf[nx] = e->tau_a_ts[t+1][1] - 1.0;
		  nx = nx+1;
		}

	      myf[nx] = e->tby_t[t][0] - p->TB_frac_ts[t][0];
	      //myf[nx] = e->tau_s_ts[t][0] - 1.0;
	      nx = nx+1;

	      if(tau_s_only_flag)
		{
		  myf[nx] = e->rir_t[t][0] - p->RIR_ts[t];
		}
	      else if(tau_s_only_flag2)
		{
		  myf[nx] = e->iy_t[t][0] - p->I_frac_ts[t][0];
		}
	      nx = nx+1;
	    }
	  else
	    {
	      // investment wedges all match investment rates
	      for(i=0; i<NC; i++)
		{
		  if(i==1 && no_rw_tau_i_flag2)
		    {
		      myf[nx] = e->tau_i_ts[t][i] - 1.0;
		    }
		  else
		    {
		      myf[nx] = e->iy_t[t][i] - p->I_frac_ts[t][i];
		    }
		  nx=nx+1;
		}

	      // first NC-1 saving wedges match saving rates
	      for(i=0; i<NC-1; i++)
		{
		  myf[nx] = e->cy_t[t][i] - p->C_frac_ts[t][i];
		  nx=nx+1;
		}

	      // last saving wedge matches interest rate
	      if(no_rw_tau_i_flag)
		{
		  myf[nx] = e->tau_i_ts[t][1] - 1.0; 
		}
	      else
		{
		  myf[nx] = e->rir_t[t][0] - p->RIR_ts[t];
		}
	      nx=nx+1;

	      if(t<NMATCH-2)
		{
		  // tfp wedge
		  for(i=0; i<NC; i++)
		    {
		      myf[nx] = e->tau_a_ts[t+1][i] - 1.0;
		      //myf[nx] = e->rgdp_t[t+1][i]/p->popt_ts[t+1][i]/e->rgdp_t[0][i] - p->RGDP_ts[t+1][i];
		      nx = nx+1;
		    }
	      
		  // trade wedge in ROW matches RER
		  myf[nx] = e->tau_t_ts[t+1][0] - 1.0;
		  nx = nx+1;
	      
		  if(tau_t_flag)
		    {
		      myf[nx] = e->reer_t[t+1][0] - p->RER_ts[t+1];
		    }
		  else
		    {
		      myf[nx] = e->tau_t_ts[t+1][1] - 1.0;
		    }
		  nx = nx+1;
		}
	    }
	}
    }
  else if(scenario==2)
    {
       for(t=0; t<NMATCH; t++)
	{
	  // investment wedges all match investment rates
	  for(i=0; i<NC; i++)
	    {
	      if(fix_tau_i[i])
		{
		  myf[nx] = e->tau_i_ts[t][i] - eee1[0].tau_i_ts[t][i];
		}
	      else
		{
		  //myf[nx] = e->tau_i_ts[t][i] - eee1[0].tau_i_ts[0][i];
		  myf[nx] = e->tau_i_ts[t][i] - 1.0;
		}
	      nx=nx+1;
	    }

	  // first NC-1 saving wedges match saving rates
	  for(i=0; i<NC; i++)
	    {
	      if(fix_tau_s[i])
		{
		  myf[nx] = e->tau_s_ts[t][i] - eee1[0].tau_s_ts[t][i];
		}
	      else
		{
		  myf[nx] = e->tau_s_ts[t][i] - 1.0;
		}

	      nx=nx+1;
	    }

	  // first NC-1 saving wedges match saving rates
	  for(i=0; i<NC; i++)
	    {
	      if(t<NMATCH-2)
		{
		  if(fix_tau_a[i])
		    {
		      myf[nx] = e->tau_a_ts[t+1][i] - eee1[0].tau_a_ts[t+1][i];
		    }
		  else
		    {
		      myf[nx] = e->tau_a_ts[t+1][i] - 1.0;
		    }
		  nx = nx+1;
		  
		  if(fix_tau_t[i])
		    {
		      myf[nx] = e->tau_t_ts[t+1][i] - eee1[0].tau_t_ts[t+1][i];
		    }
		  else
		    {
		      myf[nx] = e->tau_t_ts[t+1][i] - 1.0;
		    }
		  nx=nx+1;
		}
	    }
	}     
    }

  if(nx != neqm)
    {
      fprintf(logfile,KRED "Error evaluating eqm eqns! nx = %d, neqm = %d\n" RESET,nx,neqm);
      return 1;
    }
  
  return 0;
}

uint solve_eqm()
{
  char * sname;
  eqm * e;

  if(scenario==0)
    {
      if(NC==2)
	{
	  sname = "output/seed0_2c.bin";
	}
      else
	{
	  sname = "output/seed0_3c.bin";
	}

      if(read_seed==1)
	{

	  free_solver_mem();
	  solver_n = neqm;
	  alloc_solver_mem();

	  if(read_vec_bin(solver_x->data, neqm, sname))
	    {
	      fprintf(logfile,KRED "Error loading equilibrium guess from seed file!\n" RESET);
	      free_solver_mem();
	      return 1;
	    }
	}
      else
	{
	  if(set_initial_eqm_guess())
	    {
	      fprintf(logfile,KRED "Error constructing equilibrium guess!\n" RESET);
	      free_solver_mem();
	      return 1;
	    }
	}
    }
  // otherwise we should use the solution from the previous exercise as the initial guess
  else if(scenario==1)
    {     
      free_solver_mem();
      solver_n = neqm;
      alloc_solver_mem();

      if(NC==2)
	{
	  sname = "output/seed1_2c.bin"; 
	}
      else
	{
	  sname = "output/seed1_3c.bin"; 
	}

      if(read_seed==1)
	{
	  if(read_vec_bin(solver_x->data, neqm, sname))
	    {
	      fprintf(logfile,KRED "Error loading equilibrium guess from seed file!\n" RESET);
	      free_solver_mem();
	      return 1;
	    }	  
	}
      else
	{
	  uint it;
	  for(it=0; it<NTH; it++)
	    {
	      copy_vars( &(eee1[it]), &(eee0[0])  );
	    }
	  e = &(eee0[0]);
	  if(stack_eqm_vars(solver_x->data,e))
	    {
	      fprintf(logfile,KRED "Failed to stack variables from previous exercise!\n" RESET);
	      free_solver_mem();
	      return 1;
	    }
	}
    }
  else
    {
      free_solver_mem();
      solver_n = neqm;
      alloc_solver_mem();

      if(scenario==2)
	{
	  uint it;
	  for(it=0; it<NTH; it++)
	    {
	      copy_vars( &(eee2[it]), &(eee1[0])  );
	    }
	}
      e = &(eee1[0]);
      if(stack_eqm_vars(solver_x->data,e))
	{
	  fprintf(logfile,KRED "Failed to stack variables from previous exercise!\n" RESET);
	  free_solver_mem();
	  return 1;
	}
    }

  uint status = 0;
  if(eval_eqm_once_flag)
    {
      status = eqm_func_f(solver_x,NULL,f0[0]);
      write_vec_txt(f0[0]->data,solver_n,"output/F.txt");
      if(status)
	fprintf(logfile,KRED "Error evaluating equilibrium function!\n" RESET);
    }
  else
    {
      gsl_multiroot_function_fdf f = {&eqm_func_f,&eqm_func_df,&eqm_func_fdf,neqm,NULL};

      par=1;
      status = find_root_deriv_mkl(&f);
      if(status)
	fprintf(logfile,KRED "Error solving for equilibrium!\n" RESET);

      if(scenario<=1 && write_seed==1 && !status)
	{
	  write_vec_bin(solver_x->data, neqm, sname);
	}

    }

  free_solver_mem();

  return status;
}

///////////////////////////////

int bgp_func_f(const gsl_vector * x, void * data, gsl_vector * f)
{
   //fcnt = fcnt + 1;
  uint tn;
  if(data==NULL)
    {
      tn = 0;
    }
  else
    {
      tn = *((uint *)data);
    }
  if(eval_bgp_conds(x->data,f->data,tn))
    {
      return GSL_EBADFUNC;
    }
  else
    {
      return GSL_SUCCESS;
    }
}

int bgp_func_df(const gsl_vector * x, void * data, gsl_matrix * J)
{
  if(jacobian(&bgp_func_f, x, J, 1))
    {
      return GSL_EFAILED;
    }
  else
    {
      return GSL_SUCCESS;
    }
}

int bgp_func_fdf(const gsl_vector * x, void * data, gsl_vector * f, gsl_matrix * J)
{
  if(bgp_func_f(x,NULL,f))
    {
      return GSL_EFAILED;
    }
  else
    {
      gsl_vector_memcpy(f0[0],f);
      if(jacobian(&bgp_func_f, x, J, 0))
	{
	  return GSL_EFAILED;
	}
      else
	{
	  return GSL_SUCCESS;
	}
    }
}

int eqm_func_f(const gsl_vector * x, void * data, gsl_vector * f)
{
  //fcnt = fcnt + 1;
  uint tn;
  if(data==NULL)
    {
      tn = 0;
    }
  else
    {
      tn = *((uint *)data);
    }
  if(eval_eqm_conds(x->data,f->data,tn))
    {
      return GSL_EBADFUNC;
    }
  else
    {
      return GSL_SUCCESS;
    }
}

int eqm_func_df(const gsl_vector * x, void * data, gsl_matrix * J)
{
  if(jacobian(&eqm_func_f, x, J, 1))
    {
      return GSL_EFAILED;
    }
  else
    {
      return GSL_SUCCESS;
    }
}

int eqm_func_fdf(const gsl_vector * x, void * data, gsl_vector * f, gsl_matrix * J)
{
  if(eqm_func_f(x,NULL,f))
    {
      return GSL_EFAILED;
    }
  else
    {
      gsl_vector_memcpy(f0[0],f);
      if(jacobian(&eqm_func_f, x, J, 0))
	{
	  return GSL_EFAILED;
	}
      else
	{
	  return GSL_SUCCESS;
	}
     }
}

#endif
