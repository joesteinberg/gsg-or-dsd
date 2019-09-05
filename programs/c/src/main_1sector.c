#ifndef __MAIN_C__
#define __MAIN_C__

#include "globals_1sector.h"
#include "calibrate_1sector.h"
#include "eqm_1sector.h"

uint parse_args(int argc, char **argv)
{
  int opt = 0;
  uint cnt=0;

  while((opt = getopt(argc, argv, "ntdg")) != -1){
    switch(opt){
    case 'n':
      noio_flag=1;
      cnt++;
      break;
    case 't':
      tau_t_flag=0;
      cnt++;
      break;
    case 'd':
      no_demo_flag=1;
      cnt++;
      break;
    case 'g':
      sym_g_flag=1;
      cnt++;
      break;
    case '?':
      fprintf(logfile,"\nIncorrect command line option: %c. Possible options: -t -d -g -n\n",opt);
      fprintf(logfile,"\t-t: do not include trade wedges to match RER\n");
      fprintf(logfile,"\t-d: no demographic change\n");
      fprintf(logfile,"\t-g: constant and symmetric productivity growth\n");
      fprintf(logfile,"\t-n: no input-output linkages\n");
      fprintf(logfile,"\tOptions are mutually exclusive (choose at most one!)\n");
      return 1;
      break;
    }
  }  
  
  if(cnt>1)
    {
      fprintf(logfile,"\nIncorrect command line option: %c. Possible options: -t -d -g -n\n",opt);
      fprintf(logfile,"\t-t: do not include trade wedges to match RER\n");
      fprintf(logfile,"\t-d: no demographic change\n");
      fprintf(logfile,"\t-g: constant and symmetric productivity growth\n");
      fprintf(logfile,"\t-n: no input-output linkages\n");
      fprintf(logfile,"\tOptions are mutually exclusive (choose at most one!)\n");
      return 1;
    }

  return 0;
}

int quant_exercise()
{
  // -----------------------------------------------------------------------------------------------------------
  // set up variable and parameter structures
  uint it;
  for(it=0; it<NTH; it++)
    {
      init_vars(&(eee0[it]));
      if(calibrate(&(ppp0[it])))
	{
	  fprintf(logfile, KRED "\nProgram failed!\n" RESET);
	  return 1;
	}
      if(tau_t_flag && sym_g_flag==0 && no_demo_flag==0 && it==0 && noio_flag==0 && write_params(&ppp0[it]))
	{
	  fprintf(logfile, KRED "\nProgram failed!\n" RESET);
	  return 1;
	}

      init_vars(&(eee1[it]));
      if(calibrate(&(ppp1[it])))
	{
	  fprintf(logfile, KRED "\nProgram failed!\n" RESET);
	  return 1;
	}

      init_vars(&(eee2[it]));
      if(calibrate(&(ppp2[it])))
	{
	  fprintf(logfile, KRED "\nProgram failed!\n" RESET);
	  return 1;
	}
    }

  // -------------------------------------------------------------------------------------------------------
  // no-saving glut counterfactual
  
  fprintf(logfile, KNRM "\n----------------------------------------------------------------------\n" RESET);

  scenario = 0;
  set_neqm();

  fprintf(logfile,KBLU "\nStep 0: Solving for counterfactual with no wedges...\n" RESET);
  if(solve_eqm())
    {
      fprintf(logfile, KRED "\nProgram failed!\n" RESET);
      return 1;
    }

  write_eqm_vars(&(eee0[0]),"vars0_2c_usa_1sector",0);
  write_eqm_vars(&(eee0[0]),"vars0_2c_row_1sector",1);

  // -------------------------------------------------------------------------------------------------------
  // wedge accounting
  fprintf(logfile, KNRM "\n----------------------------------------------------------------------\n" RESET);

  scenario = 1;
  set_neqm();

  fprintf(logfile,KBLU "\nStep 1: Solving for wedges...\n" RESET);

  if(solve_eqm())
    {
      fprintf(logfile, KRED "\nProgram failed!\n" RESET);
      return 1;
    }
  write_eqm_vars(&(eee1[0]),"vars1_2c_usa_1sector",0);
  write_eqm_vars(&(eee1[0]),"vars1_2c_row_1sector",1);

  if(!cal)
    {
      // -------------------------------------------------------------------------------------------------------
      // wedge decomps
      scenario=2;
      set_neqm();
      uint cnt = 0;

      uint i;
      for(i=0; i<NC; i++)
	{
	  fix_tau_i[i]=0;
	  fix_tau_s[i]=0;
	  fix_tau_a[i]=0;
	  fix_tau_t[i]=0;
	}

      for(i=0; i<NC; i++)
	{
	  fprintf(logfile, KNRM "\n----------------------------------------------------------------------\n" RESET);
	  fprintf(logfile,KBLU "\nStep 2.%d: Fixing country %d's saving wedge...\n" RESET,cnt,i);
	  fix_tau_s[i] = 1;
	  if(solve_eqm())
	    {
	      fprintf(logfile, KRED "\nProgram failed!\n" RESET);
	      return 1;
	    }

	  char *s1 = malloc(23*sizeof(char));
	  char *s2 = malloc(23*sizeof(char));
	  snprintf(s1,23,"vars2_%d_2c_usa_1sector",cnt);
	  snprintf(s2,23,"vars2_%d_2c_row_1sector",cnt);
	  write_eqm_vars(&(eee2[0]),s1,0);
	  write_eqm_vars(&(eee2[0]),s2,1);
	  free(s1);
	  free(s2);

	  fix_tau_s[i]=0;
	  cnt++;

	  fprintf(logfile, KNRM "----------------------------------------------------------------------\n" RESET);
	  fprintf(logfile,KBLU "\nStep 2.%d: Fixing country %d's investment wedge...\n" RESET,cnt,i);
	  fix_tau_i[i] = 1;
	  if(solve_eqm())
	    {
	      fprintf(logfile, KRED "\nProgram failed!\n" RESET);
	      return 1;
	    }

	  s1 = malloc(23*sizeof(char));
	  s2 = malloc(23*sizeof(char));
	  snprintf(s1,23,"vars2_%d_2c_usa_1sector",cnt);
	  snprintf(s2,23,"vars2_%d_2c_row_1sector",cnt);
	  write_eqm_vars(&(eee2[0]),s1,0);
	  write_eqm_vars(&(eee2[0]),s2,1);
	  free(s1);
	  free(s2);

	  fix_tau_i[i]=0;
	  cnt++;
	  
	  if(i==1 && tau_t_flag)
	    {
	      fprintf(logfile, KNRM "----------------------------------------------------------------------\n" RESET);
	      fprintf(logfile,KBLU "\nStep 2.%d: Fixing country %d's trade wedge...\n" RESET,cnt,i);
	      fix_tau_t[i] = 1;
	      if(solve_eqm())
		{
		  fprintf(logfile, KRED "\nProgram failed!\n" RESET);
		  return 1;
		}

	      char *s1 = malloc(23*sizeof(char));
	      char *s2 = malloc(23*sizeof(char));
	      snprintf(s1,23,"vars2_%d_2c_usa_1sector",cnt);
	      snprintf(s2,23,"vars2_%d_2c_row_1sector",cnt);
	      write_eqm_vars(&(eee2[0]),s1,0);
	      write_eqm_vars(&(eee2[0]),s2,1);
	      free(s1);
	      free(s2);

	      fix_tau_a[i]=0;
	      cnt++;
	    }
	}      
    }

  return 0;
}

int main(int argc, char * argv[])
{
  slow_step=0;
  noio_flag=0;
  par = 0;
  solver_verbose=1;
  f_adj_cost=0;
  m_adj_cost=0;
  k_adj_cost=1;
  fixl=0;
  eval_eqm_once_flag=0;
  read_seed=1;
  write_seed=0;
  logfile = stdout;
  homotopy=1;
  hfrac=0.6;
  cal=1;
  tau_t_flag = 1;
  no_demo_flag=0;
  sym_g_flag=0;

  fprintf(logfile, KGRN "\nOn the Source of US Trade Deficits:\nGlobal Saving Glut or Domestic Saving Drought?\n" RESET);
  fprintf(logfile, KGRN "\nJoseph Steinberg, University of Toronto" RESET);
  fprintf(logfile, KNRM "\n" RESET);

  fprintf(logfile, KNRM "\n----------------------------------------------------------------------\n" RESET);
  fprintf(logfile, KBLU "\nSetting up environment...\n\n" RESET);

  if(parse_args(argc,argv))
    {
      fprintf(logfile, KRED "\nProgram failed!\n" RESET);
      return 1;
    };

  // ---------------------------------------------------------------------------------
  // set up parallel environment
#ifdef _OPENMP
  omp_set_num_threads(NTH);
  mkl_set_num_threads(NTH);
  uint nt = omp_get_max_threads();
  fprintf(logfile, KBLU "\n\tParallel processing using %d OMP threads\n" RESET,nt);
#pragma omp parallel num_threads(nt)
  {
    int it = omp_get_thread_num();
    fprintf(logfile,KBLU "\t\tHello from thread %d out of %d\n" RESET,
	    it, nt);
  }
  fprintf(logfile,"\n");
#endif

  fprintf(logfile, KNRM "\n\n\n////////////////////////////////////////////////////////////////////////\n" RESET);
  fprintf(logfile,KGRN "One-sector version\n" RESET);
  
  if(tau_t_flag==0)
    {
      fprintf(logfile,KGRN "Sensitivity analysis with notrade wedges\n\n\n" RESET);
    }
  else if(no_demo_flag==1)
    {
      fprintf(logfile,KGRN "Sensitivity analysis with no demographic change\n\n\n" RESET);
    }
  else if(sym_g_flag==1)
    {
      fprintf(logfile,KGRN "Sensitivity analysis with constant, symmetric productivity growth\n\n\n" RESET);
    }
  else if(noio_flag==1)
    {
      fprintf(logfile,KGRN "Sensitivity analysis with no input-output linkages\n\n\n" RESET);
    }
  else
    {
      fprintf(logfile,KGRN "Baseline quantitative exercise\n\n\n" RESET);
    }

  quant_exercise();

  // ---------------------------------------------------------------------------------
  fprintf(logfile, KNRM "\n----------------------------------------------------------------------\n" RESET);
  fprintf(logfile, KGRN "\nProgram complete!\n\n" RESET);

  return 0;
}

#endif
