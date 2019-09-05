import numpy as np
import pandas as pd

def colname(s,w):
    return '\\multicolumn{1}{p{'+str(w)+'cm}}{\centering '+s+'}'

def tex_header(file,loc,caption,label,ncols,w):
    file.write('\\begin{table}['+loc+']\n')
    file.write('\\footnotesize\n')
    #file.write('\\renewcommand{\\arraystretch}{1.2}\n')
    file.write('\\begin{center}\n')
    file.write('\\begin{threeparttable}')
    file.write('\\caption{'+caption+'}\n')
    file.write('\\label{tab:'+label+'}\n')
    file.write('\\begin{tabular}{l')
    for i in range(ncols):
        #file.write('C{'+str(w)+'cm}')
        file.write('c')
    file.write('}\n')
    file.write('\\toprule\n')
    return

def tex_footer(file,note):
    file.write('\\bottomrule\n')
    file.write('\\end{tabular}\n')
    file.write('\\begin{tablenotes}\n')
    file.write('\\item Notes: '+note)
    file.write('\\end{tablenotes}\n')
    file.write('\\end{threeparttable}\n')
    file.write('\\end{center}\n')
    file.write('\\normalsize\n')
    file.write('\\end{table}\n')
    return

def rmse(data,model):
    return np.sqrt(np.mean( (data-model)*(data-model) ))

def fex(counter,data,model):
    d = (data-counter)
    m = (model-counter)
    return m.sum()/d.sum()

# ------------------------------------------------------------------------------------
# load the data and model results
suff=''
capsuff = ''

# data
aggs_us = pd.read_pickle('output/key_us_data.pik').reset_index(drop=True)

# baseline (where we solve for wedges)
baseline_us = pd.read_csv('../c/output/vars1_2c_usa'+suff+'.csv')

# counterfactual without wedges
counter_us = pd.read_csv('../c/output/vars0_2c_usa'+suff+'.csv')

# fix US saving wedge
fix0_us = pd.read_csv('../c/output/vars2_0_2c_usa'+suff+'.csv')

# fix US investment wedge
fix1_us = pd.read_csv('../c/output/vars2_1_2c_usa'+suff+'.csv')

# fix RW saving wedge
fix2_us = pd.read_csv('../c/output/vars2_2_2c_usa'+suff+'.csv')

# fix RW inv wedge
fix3_us = pd.read_csv('../c/output/vars2_3_2c_usa'+suff+'.csv')

# fix RW trd wedge
fix4_us = pd.read_csv('../c/output/vars2_4_2c_usa'+suff+'.csv')

# ------------------------------------------------------------------------------------
# load the data and model results
with open('output/tab-results-1'+suff+'.tex','wb') as file:
    
    #file.write('\\begin{landscape}\n')

    # header
    tex_header(file,
               'p',
               'Wedge accounting results for U.S. quantities and prices'+capsuff,
               'tb-results-1'+suff,
               7,
               1.75)

    #top-level colnames
    file.write('& & & \\multicolumn{2}{c}{United States} & \\multicolumn{3}{c}{Rest of the world}\\\\\n')
    file.write('\\cmidrule(rl){4-5}\\cmidrule(rl){6-8}\n')
    file.write('Measure & Data & No wedges & Saving wedge & Inv. wedge & Saving wedge & Inv. wedge & Trade wedge\\\\\n')

    file.write('\\midrule\n')

    # trade balance
    file.write('\\multicolumn{8}{l}{\\textit{(a) Trade balance (percent GDP)}}\\\\\n')
    file.write(r'Minimum')
    tmp1 = baseline_us.tby[0:17].min()
    tmp2 = counter_us.tby[0:17].min()
    tmp3 = fix0_us.tby[0:17].min()
    tmp4 = fix1_us.tby[0:17].min()
    tmp5 = fix2_us.tby[0:17].min()
    tmp6 = '%0.2f' % fix3_us.tby[0:17].min()
    tmp7 = '%0.2f' % fix4_us.tby[0:17].min()
    file.write('& %0.2f & %0.2f & %0.2f & %0.2f & %0.2f & %s & %s' % (tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7))
    file.write('\\\\\n')

    file.write(r'Average')
    tmp1 = np.mean(baseline_us.tby[0:17])
    tmp2 = np.mean(counter_us.tby[0:17])
    tmp3 = np.mean(fix0_us.tby[0:17])
    tmp4 = np.mean(fix1_us.tby[0:17])
    tmp5 = np.mean(fix2_us.tby[0:17])
    tmp6 = '%0.2f' % np.mean(fix3_us.tby[0:17])
    tmp7 = '%0.2f' % np.mean(fix4_us.tby[0:17])
    file.write('& %0.2f & %0.2f & %0.2f & %0.2f & %0.2f & %s & %s' % (tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7))
    file.write('\\\\\n')

    file.write(r'RMSE from data')
    tmp2 = rmse(baseline_us.tby[0:17],counter_us.tby[0:17])
    tmp3 = rmse(baseline_us.tby[0:17],fix0_us.tby[0:17])
    tmp4 = rmse(baseline_us.tby[0:17],fix1_us.tby[0:17])
    tmp5 = rmse(baseline_us.tby[0:17],fix2_us.tby[0:17])
    tmp6 = '%0.2f' % rmse(baseline_us.tby[0:17],fix3_us.tby[0:17])
    tmp7 = '%0.2f' % rmse(baseline_us.tby[0:17],fix4_us.tby[0:17])
    file.write('& 0.00 & %0.2f & %0.2f & %0.2f & %0.2f & %s & %s' % (tmp2,tmp3,tmp4,tmp5,tmp6,tmp7))
    file.write('\\\\\n')

    file.write(r'Fraction CED explained')
    tmp3 = fex(counter_us.tby[0:17],baseline_us.tby[0:17],fix0_us.tby[0:17])
    tmp4 = fex(counter_us.tby[0:17],baseline_us.tby[0:17],fix1_us.tby[0:17])
    tmp5 = fex(counter_us.tby[0:17],baseline_us.tby[0:17],fix2_us.tby[0:17])
    tmp6 = '%0.2f' % fex(counter_us.tby[0:17],baseline_us.tby[0:17],fix3_us.tby[0:17])
    tmp7 = '%0.2f' % fex(counter_us.tby[0:17],baseline_us.tby[0:17],fix4_us.tby[0:17])
    file.write('& 1.00 & 0.00 & %0.2f & %0.2f & %0.2f & %s & %s' % (tmp3,tmp4,tmp5,tmp6,tmp7))
    file.write('\\\\\n')
    file.write('\\\\\n')



    file.write('\\multicolumn{8}{l}{\\textit{(b) Real exchange rate (1995 data = 100)}}\\\\\n')
    file.write(r'Minimum')
    tmp1 = baseline_us.reer[0:17].min()
    tmp2 = counter_us.reer[0:17].min()
    tmp3 = fix0_us.reer[0:17].min()
    tmp4 = fix1_us.reer[0:17].min()
    tmp5 = fix2_us.reer[0:17].min()
    tmp6 = '%0.2f' % fix3_us.reer[0:17].min()
    tmp7 = '%0.2f'% fix4_us.reer[0:17].min()
    file.write('& %0.2f & %0.2f & %0.2f & %0.2f & %0.2f & %s & %s' % (tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7))
    file.write('\\\\\n')

    file.write(r'Average')
    tmp1 = np.mean(baseline_us.reer[0:17])
    tmp2 = np.mean(counter_us.reer[0:17])
    tmp3 = np.mean(fix0_us.reer[0:17])
    tmp4 = np.mean(fix1_us.reer[0:17])
    tmp5 = np.mean(fix2_us.reer[0:17])
    tmp6 = '%0.2f' % np.mean(fix3_us.reer[0:17])
    tmp7 = '%0.2f' % np.mean(fix4_us.reer[0:17])
    file.write('& %0.2f & %0.2f & %0.2f & %0.2f & %0.2f & %s & %s' % (tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7))
    file.write('\\\\\n')

    file.write(r'RMSE from data')
    tmp2 = rmse(baseline_us.reer[0:17],counter_us.reer[0:17])
    tmp3 = rmse(baseline_us.reer[0:17],fix0_us.reer[0:17])
    tmp4 = rmse(baseline_us.reer[0:17],fix1_us.reer[0:17])
    tmp5 = rmse(baseline_us.reer[0:17],fix2_us.reer[0:17])
    tmp6 = '%0.2f' % rmse(baseline_us.reer[0:17],fix3_us.reer[0:17])
    tmp7 = '%0.2f' % rmse(baseline_us.reer[0:17],fix4_us.reer[0:17])
    file.write('& 0.00 & %0.2f & %0.2f & %0.2f & %0.2f & %s & %s' % (tmp2,tmp3,tmp4,tmp5,tmp6,tmp7))
    file.write('\\\\\n')

    file.write(r'Fraction CEA explained')
    tmp3 = fex(counter_us.reer[0:17],baseline_us.reer[0:17],fix0_us.reer[0:17])
    tmp4 = fex(counter_us.reer[0:17],baseline_us.reer[0:17],fix1_us.reer[0:17])
    tmp5 = fex(counter_us.reer[0:17],baseline_us.reer[0:17],fix2_us.reer[0:17])
    tmp6 = '%0.2f' % fex(counter_us.reer[0:17],baseline_us.reer[0:17],fix3_us.reer[0:17])
    tmp7 = '%0.2f' % fex(counter_us.reer[0:17],baseline_us.reer[0:17],fix4_us.reer[0:17])
    file.write('& 1.00 & 0.00 & %0.2f & %0.2f & %0.2f & %s & %s' % (tmp3,tmp4,tmp5,tmp6,tmp7))
    file.write('\\\\\n')

    file.write('\\\\\n')




    # investment
    file.write('\\multicolumn{8}{l}{\\textit{(c) Investment rate (percent GDP)}}\\\\\n')
    file.write(r'RMSE from data')
    tmp2 = rmse(baseline_us.iy[0:17],counter_us.iy[0:17])
    tmp3 = rmse(baseline_us.iy[0:17],fix0_us.iy[0:17])
    tmp4 = rmse(baseline_us.iy[0:17],fix1_us.iy[0:17])
    tmp5 = rmse(baseline_us.iy[0:17],fix2_us.iy[0:17])
    tmp6 = '%0.2f' % rmse(baseline_us.iy[0:17],fix3_us.iy[0:17])
    tmp7 = '%0.2f' % rmse(baseline_us.iy[0:17],fix4_us.iy[0:17])
    file.write('& 0.00 & %0.2f & %0.2f & %0.2f & %0.2f & %s & %s' % (tmp2,tmp3,tmp4,tmp5,tmp6,tmp7))
    file.write('\\\\\n')

    file.write(r'Fraction CEI explained (1995--2006)')
    tmp3 = fex(counter_us.iy[0:12],baseline_us.iy[0:12],fix0_us.iy[0:12])
    tmp4 = fex(counter_us.iy[0:12],baseline_us.iy[0:12],fix1_us.iy[0:12])
    tmp5 = fex(counter_us.iy[0:12],baseline_us.iy[0:12],fix2_us.iy[0:12])
    tmp6 = '%0.2f' % fex(counter_us.iy[0:12],baseline_us.iy[0:12],fix3_us.iy[0:12])
    tmp7 = '%0.2f' % fex(counter_us.iy[0:12],baseline_us.iy[0:12],fix4_us.iy[0:12])
    file.write('& 1.00 & 0.00 & %0.2f & %0.2f & %0.2f & %s & %s' % (tmp3,tmp4,tmp5,tmp6,tmp7))
    file.write('\\\\\n')

    file.write(r'Fraction CEI explained (2007--2011)')
    tmp3 = fex(counter_us.iy[12:17],baseline_us.iy[12:17],fix0_us.iy[12:17])
    tmp4 = fex(counter_us.iy[12:17],baseline_us.iy[12:17],fix1_us.iy[12:17])
    tmp5 = fex(counter_us.iy[12:17],baseline_us.iy[12:17],fix2_us.iy[12:17])
    tmp6 = '%0.2f' % fex(counter_us.iy[12:17],baseline_us.iy[12:17],fix3_us.iy[12:17])
    tmp7 = '%0.2f' % fex(counter_us.iy[12:17],baseline_us.iy[12:17],fix4_us.iy[12:17])
    file.write('& 1.00 & 0.00 & %0.2f & %0.2f & %0.2f & %s & %s' % (tmp3,tmp4,tmp5,tmp6,tmp7))
    file.write('\\\\\n')

    file.write('\\\\\n')

    # real interest rate
    file.write('\\multicolumn{8}{l}{\\textit{(d) Real interest rate (percent per year)}}\\\\\n')
    file.write(r'Minimum')
    tmp1 = baseline_us.rir[0:17].min()
    tmp2 = counter_us.rir[0:17].min()
    tmp3 = fix0_us.rir[0:17].min()
    tmp4 = fix1_us.rir[0:17].min()
    tmp5 = fix2_us.rir[0:17].min()
    tmp6 = '%0.2f' % fix3_us.rir[0:17].min()
    tmp7 = '%0.2f' % fix4_us.rir[0:17].min()
    file.write('& %0.2f & %0.2f & %0.2f & %0.2f & %0.2f & %s & %s' % (tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7))
    file.write('\\\\\n')

    file.write(r'Average')
    tmp1 = np.mean(baseline_us.rir[0:17])
    tmp2 = np.mean(counter_us.rir[0:17])
    tmp3 = np.mean(fix0_us.rir[0:17])
    tmp4 = np.mean(fix1_us.rir[0:17])
    tmp5 = np.mean(fix2_us.rir[0:17])
    tmp6 = '%0.2f' % np.mean(fix3_us.rir[0:17])
    tmp7 = '%0.2f' % np.mean(fix4_us.rir[0:17])
    file.write('& %0.2f & %0.2f & %0.2f & %0.2f & %0.2f & %s & %s' % (tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7))
    file.write('\\\\\n')

    file.write(r'RMSE from data')
    tmp2 = rmse(baseline_us.rir[0:17],counter_us.rir[0:17])
    tmp3 = rmse(baseline_us.rir[0:17],fix0_us.rir[0:17])
    tmp4 = rmse(baseline_us.rir[0:16],fix1_us.rir[0:17])
    tmp5 = rmse(baseline_us.rir[0:17],fix2_us.rir[0:17])
    tmp6 = '%0.2f' % rmse(baseline_us.rir[0:17],fix3_us.rir[0:17])
    tmp7 = '%0.2f' % rmse(baseline_us.rir[0:17],fix4_us.rir[0:17])
    file.write('& 0.00 & %0.2f & %0.2f & %0.2f & %0.2f & %s & %s' % (tmp2,tmp3,tmp4,tmp5,tmp6,tmp7))
    file.write('\\\\\n')

    # file.write(r'Correlation with data')
    # tmp2 = np.corrcoef(baseline_us.rir[0:17],counter_us.rir[0:17])[0,1]
    # tmp3 = np.corrcoef(baseline_us.rir[0:17],fix0_us.rir[0:17])[0,1]
    # tmp4 = np.corrcoef(baseline_us.rir[0:16],fix1_us.rir[0:16])[0,1]
    # tmp5 = np.corrcoef(baseline_us.rir[0:17],fix2_us.rir[0:17])[0,1]
    # tmp6 = '%0.2f' % np.corrcoef(baseline_us.rir[0:17],fix3_us.rir[0:17])[0,1]
    # tmp7 = '%0.2f' % np.corrcoef(baseline_us.rir[0:17],fix4_us.rir[0:17])[0,1]
    # file.write('& 1.00 & %0.2f & %0.2f & %0.2f & %0.2f & %s & %s' % (tmp2,tmp3,tmp4,tmp5,tmp6,tmp7))
    # file.write('\\\\\n')

    file.write(r'Fraction decline explained')
    tmp1 = baseline_us.rir[16]-baseline_us.rir[0]
    tmp2 = (counter_us.rir[16]-counter_us.rir[0])/tmp1
    tmp3 = (fix0_us.rir[16]-fix0_us.rir[0])/tmp1
    tmp4 = (fix1_us.rir[15]-fix1_us.rir[0])/tmp1
    tmp5 = (fix2_us.rir[16]-fix2_us.rir[0])/tmp1
    tmp6 = '%0.2f' % ((fix3_us.rir[16]-fix3_us.rir[0])/tmp1)
    tmp7 = '%0.2f' % ((fix4_us.rir[16]-fix4_us.rir[0])/tmp1)
    tmp1=tmp1/tmp1
    file.write('& %0.2f & %0.2f & %0.2f & %0.2f & %0.2f & %s & %s' % (tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7))
    file.write('\\\\\n')

    note='The second column reports counterfactual model outcomes when all wedges are set to one. Columns 3--7 report counterfactual model outcomes with one wedge set to its calibrated value in each period and all other wedges held constant. Fraction of CED explained calculated as (cumulative difference between trade balance in model and no wedge counterfactual) divided by (cumulative difference between trade balance in data and no-wedge counterfactual). Fraction of cumulative excess RER appreciation (CEA) and cumulative excess investment (CEI) computed analogously.'
    tex_footer(file,note)
    
    #file.write('\\end{landscape}\n')
    
