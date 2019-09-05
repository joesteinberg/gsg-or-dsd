///////////////////////////////////////////////////////////////////////////////////////////////

"On the Source of U.S. Trade Deficits: Global Saving Glut or Domestic Saving Drought?"

Joseph B. Steinberg

University of Toronto

RED Manuscript RED-16-198R1

This document describes the data and programs needed to replicate my empirical and quantitative 

analysis. The document is organized into three sections:


0. Datasets

1. Python scripts for data processing 
2. C program for quantitative analysis
3. Python scripts for 



The python scripts in section 1 should be run before the C program in section 2, which should in
turn be run before the python scripts in section 3.

All of my analyses were performed on a Dell Workstation 
running Ubuntu Linux 12.04. I use Python
version 2 with the SciPy, NumPy, and Matplotlib packages. The C program requires the GSL library
and the Intel MKL, and takes about 10 minutes to run in a 12-core, 32-gigabyte Dell Workstation
running Ubuntu Linux 12.04. There are no random numbers generated in any of the parts of the 
analyses.

///////////////////////////////////////////////////////////////////////////////////////////////

0. Datasets

0.1. wiot_full.dta: full World Input Output Database dataset (WIOD).

0.2. EWN19702011.xlsx: Lane and Milesi-Feretti dataset (EWN).

0.3. ted.csv: Total Economy Database dataset (TED).

0.4. pwt90.dta: Penn World Tables 9.0 (PWT).

0.5. WPP2015_INT_F2A_Annual_Population_Indicators: UN World Population Prospects total
population data (UNWPP).

0.6. WPP2015_INT_F2A_Annual_Population_Indicators_DependencyRatios: UN World Population
Prospects dependency ratios (UNWPP).

0.7. kaopen.dta: Chinn-Ito capital account openness dataset (KAOPEN).

0.8. FinStructure_November_2013.csv: Beck et al. dataset (BECK).

0.9. wpp_codes.csv: mapping between country isocodes and numeric country codes in UNWPP data.

0.10. forpython.csv: csv file containing US real exchange rate (REER from IMF IFS) and real
interest rate (10-year nominal T-bond yield less CPI-U inflation from FRED).

///////////////////////////////////////////////////////////////////////////////////////////////

1. Data processing scripts located "programs/python" (run in the following order)

1.1. wiod_dta_to_pik.py: Converts the WIOD Stata dataset into a Pandas DataFrame and pickles it.

1.2. wiod_preproc.py: aggregates WIOD data according to scheme described in paper, writes output
files used by C program.

1.3. wiod_weights.py: uses the disaggregated WIOD data to compute the time-varying weights used
to construct the RoW average.

1.4. iomat.py: uses the aggregated WIOD data to construct the IO tables used in the calibration.
Also creates Table 1.

1.5. lp.py: uses the WIOD weights and the TED data to compute productivity time series.

1.6. demo.py: uses the WIOD weights and the UNWPP data to compute demographic time series.

1.7. k.py: uses the WIOD weights and PWT data to compute initial capital stocks.

1.8. kaopen_findev.py: uses the WIOD weights, the KAOPEN data, and the Beck data to compute
time series for the rest of the world's capital account openness and domestic financial
development.

///////////////////////////////////////////////////////////////////////////////////////////////

2. C program located in "programs/c"

The C program that performs the quantitative analysis has a number of source files, all located
in the "src" subfolder.

2.1. main.c: main function
2.2. globals.h: macros, flags, simple utilities
2.3. calibrate.h/calibrate.c: calibration routines and input data processing
2.4. eqm.h/eqm.c: equilibrium conditions and driver for equilibrium solver
2.5. solver.h/solver.c: nonlinear system solver with parallelized jacobian evaluation
2.6. gnewton.h/gnewton.c: augmented newton's method iterator

There is also a makefile in the main C program folder. I use gcc with the Intel MKL on
Ubuntu Linux. The GSL library is also required. To compile, simply type "make". To run the
baeline model, type "./bin/gsg_dsd". There are a number of command-line options that trigger
sensitivity analyses. Type "./bin/gsg_dsd -h" to see these options. The program runs in about
10 minutes on a Dell Workstation with 12 cores and 32 gigabytes of RAM. Note that it takes about
10 times as long without parallelization.

The program uses a number of input files created by the scripts above. All of these files
are located in the folder "programs/python/output."

The program writes a number of output files that are used by the python scripts in the next
section. They are located in the "output" subfolder.

///////////////////////////////////////////////////////////////////////////////////////////////

3. Scripts to make tables and figures located "programs/python"

3.1. 2period_analytical_results.py: creates Figure 3.
3.2. plots.py: creates all other figures.
3.3. tables.py: creates Table 4.
3.4. tables_sens.py: creates Table 5 and Table 1 in the appendix.