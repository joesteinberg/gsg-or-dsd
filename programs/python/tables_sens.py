import numpy as np
import pandas as pd

def colname(s,w):
    return '\\multicolumn{1}{p{'+str(w)+'cm}}{\centering '+s+'}'

def tex_header(file,loc,caption,label,ncols,w):
    file.write('\\begin{table}['+loc+']\n')
    file.write('\\footnotesize\n')
    file.write('\\begin{center}\n')
    file.write('\\caption{'+caption+'}\n')
    file.write('\\label{tab:'+label+'}\n')
    file.write('\\begin{threeparttable}')
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

# data
aggs_us = pd.read_pickle('output/key_us_data.pik').reset_index(drop=True)

# ------------------------------------------------------------------------------------
# load the data and model results

# data
aggs_us = pd.read_pickle('output/key_us_data.pik').reset_index(drop=True)

# baseline (where we solve for wedges)
baseline_us = pd.read_csv('../c/output/vars1_2c_usa.csv')
baseline_us_symg = pd.read_csv('../c/output/vars1_2c_usa_symg.csv')
baseline_us_nodemo = pd.read_csv('../c/output/vars1_2c_usa_nodemo.csv')
baseline_us_notrdw = pd.read_csv('../c/output/vars1_2c_usa_notrdw.csv')
baseline_us_noio = pd.read_csv('../c/output/vars1_2c_usa_noio.csv')
baseline_us_etaKh = pd.read_csv('../c/output/vars1_2c_usa_etaKh.csv')
baseline_us_etaKl = pd.read_csv('../c/output/vars1_2c_usa_etaKl.csv')
baseline_us_armh = pd.read_csv('../c/output/vars1_2c_usa_armh.csv')
baseline_us_arml = pd.read_csv('../c/output/vars1_2c_usa_arml.csv')
baseline_us_rhowh = pd.read_csv('../c/output/vars1_2c_usa_rhowh.csv')
baseline_us_rhowl = pd.read_csv('../c/output/vars1_2c_usa_rhowl.csv')
baseline_us_1sector = pd.read_csv('../c/output/vars1_2c_usa_1sector.csv')
baseline_us_lowr = pd.read_csv('../c/output/vars1_2c_usa_lowr.csv')

# counterfactual without wedges
counter_us = pd.read_csv('../c/output/vars0_2c_usa.csv')
counter_us_symg = pd.read_csv('../c/output/vars0_2c_usa_symg.csv')
counter_us_nodemo = pd.read_csv('../c/output/vars0_2c_usa_nodemo.csv')
counter_us_notrdw = pd.read_csv('../c/output/vars0_2c_usa_notrdw.csv')
counter_us_noio = pd.read_csv('../c/output/vars0_2c_usa_noio.csv')
counter_us_etaKh = pd.read_csv('../c/output/vars0_2c_usa_etaKh.csv')
counter_us_etaKl = pd.read_csv('../c/output/vars0_2c_usa_etaKl.csv')
counter_us_armh = pd.read_csv('../c/output/vars0_2c_usa_armh.csv')
counter_us_arml = pd.read_csv('../c/output/vars0_2c_usa_arml.csv')
counter_us_rhowh = pd.read_csv('../c/output/vars0_2c_usa_rhowh.csv')
counter_us_rhowl = pd.read_csv('../c/output/vars0_2c_usa_rhowl.csv')
counter_us_1sector = pd.read_csv('../c/output/vars0_2c_usa_1sector.csv')
counter_us_lowr = pd.read_csv('../c/output/vars0_2c_usa_lowr.csv')

# fix US saving wedge
fix0_us = pd.read_csv('../c/output/vars2_0_2c_usa.csv')
fix0_us_symg = pd.read_csv('../c/output/vars2_0_2c_usa_symg.csv')
fix0_us_nodemo = pd.read_csv('../c/output/vars2_0_2c_usa_nodemo.csv')
fix0_us_notrdw = pd.read_csv('../c/output/vars2_0_2c_usa_notrdw.csv')
fix0_us_noio = pd.read_csv('../c/output/vars2_0_2c_usa_noio.csv')
fix0_us_etaKh = pd.read_csv('../c/output/vars2_0_2c_usa_etaKh.csv')
fix0_us_etaKl = pd.read_csv('../c/output/vars2_0_2c_usa_etaKl.csv')
fix0_us_armh = pd.read_csv('../c/output/vars2_0_2c_usa_armh.csv')
fix0_us_arml = pd.read_csv('../c/output/vars2_0_2c_usa_arml.csv')
fix0_us_rhowh = pd.read_csv('../c/output/vars2_0_2c_usa_rhowh.csv')
fix0_us_rhowl = pd.read_csv('../c/output/vars2_0_2c_usa_rhowl.csv')
fix0_us_1sector = pd.read_csv('../c/output/vars2_0_2c_usa_1sector.csv')
fix0_us_lowr = pd.read_csv('../c/output/vars2_0_2c_usa_lowr.csv')

# fix US inv wedge
fix1_us = pd.read_csv('../c/output/vars2_1_2c_usa.csv')
fix1_us_symg = pd.read_csv('../c/output/vars2_1_2c_usa_symg.csv')
fix1_us_nodemo = pd.read_csv('../c/output/vars2_1_2c_usa_nodemo.csv')
fix1_us_notrdw = pd.read_csv('../c/output/vars2_1_2c_usa_notrdw.csv')
fix1_us_noio = pd.read_csv('../c/output/vars2_1_2c_usa_noio.csv')
fix1_us_etaKh = pd.read_csv('../c/output/vars2_1_2c_usa_etaKh.csv')
fix1_us_etaKl = pd.read_csv('../c/output/vars2_1_2c_usa_etaKl.csv')
fix1_us_armh = pd.read_csv('../c/output/vars2_1_2c_usa_armh.csv')
fix1_us_arml = pd.read_csv('../c/output/vars2_1_2c_usa_arml.csv')
fix1_us_rhowh = pd.read_csv('../c/output/vars2_1_2c_usa_rhowh.csv')
fix1_us_rhowl = pd.read_csv('../c/output/vars2_1_2c_usa_rhowl.csv')
fix1_us_1sector = pd.read_csv('../c/output/vars2_1_2c_usa_1sector.csv')
fix1_us_lowr = pd.read_csv('../c/output/vars2_1_2c_usa_lowr.csv')

# fix RW saving wedge
fix2_us = pd.read_csv('../c/output/vars2_2_2c_usa.csv')
fix2_us_symg = pd.read_csv('../c/output/vars2_2_2c_usa_symg.csv')
fix2_us_nodemo = pd.read_csv('../c/output/vars2_2_2c_usa_nodemo.csv')
fix2_us_notrdw = pd.read_csv('../c/output/vars2_2_2c_usa_notrdw.csv')
fix2_us_noio = pd.read_csv('../c/output/vars2_2_2c_usa_noio.csv')
fix2_us_etaKh = pd.read_csv('../c/output/vars2_2_2c_usa_etaKh.csv')
fix2_us_etaKl = pd.read_csv('../c/output/vars2_2_2c_usa_etaKl.csv')
fix2_us_armh = pd.read_csv('../c/output/vars2_2_2c_usa_armh.csv')
fix2_us_arml = pd.read_csv('../c/output/vars2_2_2c_usa_arml.csv')
fix2_us_rhowh = pd.read_csv('../c/output/vars2_2_2c_usa_rhowh.csv')
fix2_us_rhowl = pd.read_csv('../c/output/vars2_2_2c_usa_rhowl.csv')
fix2_us_1sector = pd.read_csv('../c/output/vars2_2_2c_usa_1sector.csv')
fix2_us_lowr = pd.read_csv('../c/output/vars2_2_2c_usa_lowr.csv')

# fix RW inv
fix3_us = pd.read_csv('../c/output/vars2_3_2c_usa.csv')
fix3_us_symg = pd.read_csv('../c/output/vars2_3_2c_usa_symg.csv')
fix3_us_nodemo = pd.read_csv('../c/output/vars2_3_2c_usa_nodemo.csv')
fix3_us_notrdw = pd.read_csv('../c/output/vars2_3_2c_usa_notrdw.csv')
fix3_us_noio = pd.read_csv('../c/output/vars2_3_2c_usa_noio.csv')
fix3_us_etaKh = pd.read_csv('../c/output/vars2_3_2c_usa_etaKh.csv')
fix3_us_etaKl = pd.read_csv('../c/output/vars2_3_2c_usa_etaKl.csv')
fix3_us_armh = pd.read_csv('../c/output/vars2_3_2c_usa_armh.csv')
fix3_us_arml = pd.read_csv('../c/output/vars2_3_2c_usa_arml.csv')
fix3_us_rhowh = pd.read_csv('../c/output/vars2_3_2c_usa_rhowh.csv')
fix3_us_rhowl = pd.read_csv('../c/output/vars2_3_2c_usa_rhowl.csv')
fix3_us_1sector = pd.read_csv('../c/output/vars2_3_2c_usa_1sector.csv')
fix3_us_lowr = pd.read_csv('../c/output/vars2_3_2c_usa_lowr.csv')

# fix RW trade wedge
fix4_us = pd.read_csv('../c/output/vars2_4_2c_usa.csv')
fix4_us_symg = pd.read_csv('../c/output/vars2_4_2c_usa_symg.csv')
fix4_us_nodemo = pd.read_csv('../c/output/vars2_4_2c_usa_nodemo.csv')
fix4_us_noio = pd.read_csv('../c/output/vars2_4_2c_usa_noio.csv')
fix4_us_etaKh = pd.read_csv('../c/output/vars2_4_2c_usa_etaKh.csv')
fix4_us_etaKl = pd.read_csv('../c/output/vars2_4_2c_usa_etaKl.csv')
fix4_us_armh = pd.read_csv('../c/output/vars2_4_2c_usa_armh.csv')
fix4_us_arml = pd.read_csv('../c/output/vars2_4_2c_usa_arml.csv')
fix4_us_rhowh = pd.read_csv('../c/output/vars2_4_2c_usa_rhowh.csv')
fix4_us_rhowl = pd.read_csv('../c/output/vars2_4_2c_usa_rhowl.csv')
fix4_us_1sector = pd.read_csv('../c/output/vars2_4_2c_usa_1sector.csv')
fix4_us_lowr = pd.read_csv('../c/output/vars2_4_2c_usa_lowr.csv')

# ------------------------------------------------------------------------------------
# load the data and model results
with open('output/tab-results-2.tex','wb') as file:
    
    #file.write('\\begin{landscape}\n')

    # header
    tex_header(file,
               'p',
               'Wedge accounting for U.S. trade balance in baseline model and sensitivity analyses',
               'tb-results-2',
               6,
               1.75)

    #top-level colnames
    file.write('& & \\multicolumn{2}{c}{United States} & \\multicolumn{3}{c}{Rest of the world}\\\\\n')
    file.write('\\cmidrule(rl){3-4}\\cmidrule(rl){5-7}\n')
    file.write('Model & No wedges & Saving wedge & Inv. wedge & Saving wedge & Inv. wedge & Trade wedge\\\\\n')

    file.write('\\midrule\n')

    # trade balance
    file.write('\\multicolumn{7}{l}{\\textit{(a) RMSE from data}}\\\\\n')

    file.write(r'Baseline')
    tmp2 = rmse(baseline_us.tby[0:17],counter_us.tby[0:17])
    tmp3 = rmse(baseline_us.tby[0:17],fix0_us.tby[0:17])
    tmp4 = rmse(baseline_us.tby[0:17],fix1_us.tby[0:17])
    tmp5 = rmse(baseline_us.tby[0:17],fix2_us.tby[0:17])
    tmp6 = rmse(baseline_us.tby[0:17],fix3_us.tby[0:17])
    tmp7 = rmse(baseline_us.tby[0:17],fix4_us.tby[0:17])
    file.write('& %0.2f & %0.2f & %0.2f & %0.2f & %0.2f & %0.2f' % (tmp2,tmp3,tmp4,tmp5,tmp6,tmp7))
    file.write('\\\\\n')
    file.write(r'No demographic change')
    tmp2 = rmse(baseline_us_nodemo.tby[0:17],counter_us_nodemo.tby[0:17])
    tmp3 = rmse(baseline_us_nodemo.tby[0:17],fix0_us_nodemo.tby[0:17])
    tmp4 = rmse(baseline_us_nodemo.tby[0:17],fix1_us_nodemo.tby[0:17])
    tmp5 = rmse(baseline_us_nodemo.tby[0:17],fix2_us_nodemo.tby[0:17])
    tmp6 = rmse(baseline_us_nodemo.tby[0:17],fix3_us_nodemo.tby[0:17])
    tmp7 = rmse(baseline_us_nodemo.tby[0:17],fix4_us_nodemo.tby[0:17])
    file.write('& %0.2f & %0.2f & %0.2f & %0.2f & %0.2f & %0.2f' % (tmp2,tmp3,tmp4,tmp5,tmp6,tmp7))
    file.write('\\\\\n')
    file.write(r'Symmetric growth')
    tmp2 = rmse(baseline_us_symg.tby[0:17],counter_us_symg.tby[0:17])
    tmp3 = rmse(baseline_us_symg.tby[0:17],fix0_us_symg.tby[0:17])
    tmp4 = rmse(baseline_us_symg.tby[0:17],fix1_us_symg.tby[0:17])
    tmp5 = rmse(baseline_us_symg.tby[0:17],fix2_us_symg.tby[0:17])
    tmp6 = rmse(baseline_us_symg.tby[0:17],fix3_us_symg.tby[0:17])
    tmp7 = rmse(baseline_us_symg.tby[0:17],fix4_us_symg.tby[0:17])
    file.write('& %0.2f & %0.2f & %0.2f & %0.2f & %0.2f & %0.2f' % (tmp2,tmp3,tmp4,tmp5,tmp6,tmp7))
    file.write('\\\\\n')
    file.write(r'No trade wedge')
    tmp2 = rmse(baseline_us_notrdw.tby[0:17],counter_us_notrdw.tby[0:17])
    tmp3 = rmse(baseline_us_notrdw.tby[0:17],fix0_us_notrdw.tby[0:17])
    tmp4 = rmse(baseline_us_notrdw.tby[0:17],fix1_us_notrdw.tby[0:17])
    tmp5 = rmse(baseline_us_notrdw.tby[0:17],fix2_us_notrdw.tby[0:17])
    tmp6 = rmse(baseline_us_notrdw.tby[0:17],fix3_us_notrdw.tby[0:17])
    file.write('& %0.2f & %0.2f & %0.2f & %0.2f & %0.2f & -' % (tmp2,tmp3,tmp4,tmp5,tmp6))
    file.write('\\\\\n')

    file.write('\\\\\n')

    file.write('\\multicolumn{7}{l}{\\textit{(b) Fraction of CED explained}}\\\\\n')
    file.write(r'Baseline')
    tmp3 = fex(counter_us.tby[0:17],baseline_us.tby[0:17],fix0_us.tby[0:17])
    tmp4 = fex(counter_us.tby[0:17],baseline_us.tby[0:17],fix1_us.tby[0:17])
    tmp5 = fex(counter_us.tby[0:17],baseline_us.tby[0:17],fix2_us.tby[0:17])
    tmp6 = fex(counter_us.tby[0:17],baseline_us.tby[0:17],fix3_us.tby[0:17])
    tmp7 = fex(counter_us.tby[0:17],baseline_us.tby[0:17],fix4_us.tby[0:17])
    file.write('& 0.00 & %0.2f & %0.2f & %0.2f & %0.2f & %0.2f' % (tmp3,tmp4,tmp5,tmp6,tmp7))
    file.write('\\\\\n')
    file.write(r'Symmetric growth')
    tmp3 = fex(counter_us_symg.tby[0:17],baseline_us_symg.tby[0:17],fix0_us_symg.tby[0:17])
    tmp4 = fex(counter_us_symg.tby[0:17],baseline_us_symg.tby[0:17],fix1_us_symg.tby[0:17])
    tmp5 = fex(counter_us_symg.tby[0:17],baseline_us_symg.tby[0:17],fix2_us_symg.tby[0:17])
    tmp6 = fex(counter_us_symg.tby[0:17],baseline_us_symg.tby[0:17],fix3_us_symg.tby[0:17])
    tmp7 = fex(counter_us_symg.tby[0:17],baseline_us_symg.tby[0:17],fix4_us_symg.tby[0:17])
    file.write('& 0.00 & %0.2f & %0.2f & %0.2f & %0.2f & %0.2f' % (tmp3,tmp4,tmp5,tmp6,tmp7))
    file.write('\\\\\n')
    file.write(r'No demographic change')
    tmp3 = fex(counter_us_nodemo.tby[0:17],baseline_us_nodemo.tby[0:17],fix0_us_nodemo.tby[0:17])
    tmp4 = fex(counter_us_nodemo.tby[0:17],baseline_us_nodemo.tby[0:17],fix1_us_nodemo.tby[0:17])
    tmp5 = fex(counter_us_nodemo.tby[0:17],baseline_us_nodemo.tby[0:17],fix2_us_nodemo.tby[0:17])
    tmp6 = fex(counter_us_nodemo.tby[0:17],baseline_us_nodemo.tby[0:17],fix3_us_nodemo.tby[0:17])
    tmp7 = fex(counter_us_nodemo.tby[0:17],baseline_us_nodemo.tby[0:17],fix4_us_nodemo.tby[0:17])
    file.write('& 0.00 & %0.2f & %0.2f & %0.2f & %0.2f & %0.2f' % (tmp3,tmp4,tmp5,tmp6,tmp7))
    file.write('\\\\\n')
    file.write(r'No trade wedge')
    tmp3 = fex(counter_us_notrdw.tby[0:17],baseline_us_notrdw.tby[0:17],fix0_us_notrdw.tby[0:17])
    tmp4 = fex(counter_us_notrdw.tby[0:17],baseline_us_notrdw.tby[0:17],fix1_us_notrdw.tby[0:17])
    tmp5 = fex(counter_us_notrdw.tby[0:17],baseline_us_notrdw.tby[0:17],fix2_us_notrdw.tby[0:17])
    tmp6 = fex(counter_us_notrdw.tby[0:17],baseline_us_notrdw.tby[0:17],fix3_us_notrdw.tby[0:17])
    file.write('& 0.00 & %0.2f & %0.2f & %0.2f & %0.2f & -' % (tmp3,tmp4,tmp5,tmp6))
    file.write('\\\\\n')
    note='The second column reports counterfactual model outcomes when all wedges are set to one. Columns 3--7 report counterfactual model outcomes with one wedge set to its calibrate value in each period and all other wedges held constant. Fraction of CED explained calculated as (cumulative difference between trade balance in model and no wedge counterfactual) divided by (cumulative difference between trade balance in data and no-wedge counterfactual).'

    tex_footer(file,note)
    
    #file.write('\\end{landscape}\n')














# ------------------------------------------------------------------------------------
# load the data and model results
with open('output/tab-results-extra.tex','wb') as file:
    
    #file.write('\\begin{landscape}\n')

    # header
    tex_header(file,
               'p',
               'Wedge accounting for U.S. trade balance in baseline model and additional sensitivity analyses',
               'tb-results-2',
               6,
               1.75)

    #top-level colnames
    file.write('& & \\multicolumn{2}{c}{United States} & \\multicolumn{3}{c}{Rest of the world}\\\\\n')
    file.write('\\cmidrule(rl){3-4}\\cmidrule(rl){5-7}\n')
    file.write('Model & No wedges & Saving wedge & Inv. wedge & Saving wedge & Inv. wedge & Trade wedge\\\\\n')

    file.write('\\midrule\n')

    # trade balance
    file.write('\\multicolumn{7}{l}{\\textit{(a) RMSE from data}}\\\\\n')

    file.write(r'Baseline')
    tmp2 = rmse(baseline_us.tby[0:17],counter_us.tby[0:17])
    tmp3 = rmse(baseline_us.tby[0:17],fix0_us.tby[0:17])
    tmp4 = rmse(baseline_us.tby[0:17],fix1_us.tby[0:17])
    tmp5 = rmse(baseline_us.tby[0:17],fix2_us.tby[0:17])
    tmp6 = rmse(baseline_us.tby[0:17],fix3_us.tby[0:17])
    tmp7 = rmse(baseline_us.tby[0:17],fix4_us.tby[0:17])
    file.write('& %0.2f & %0.2f & %0.2f & %0.2f & %0.2f & %0.2f' % (tmp2,tmp3,tmp4,tmp5,tmp6,tmp7))
    file.write('\\\\\n')
    file.write(r'No input-output linkages')
    tmp2 = rmse(baseline_us_noio.tby[0:17],counter_us_noio.tby[0:17])
    tmp3 = rmse(baseline_us_noio.tby[0:17],fix0_us_noio.tby[0:17])
    tmp4 = rmse(baseline_us_noio.tby[0:17],fix1_us_noio.tby[0:17])
    tmp5 = rmse(baseline_us_noio.tby[0:17],fix2_us_noio.tby[0:17])
    tmp6 = rmse(baseline_us_noio.tby[0:17],fix3_us_noio.tby[0:17])
    tmp7 = rmse(baseline_us_noio.tby[0:17],fix4_us_noio.tby[0:17])
    file.write('& %0.2f & %0.2f & %0.2f & %0.2f & %0.2f & %0.2f' % (tmp2,tmp3,tmp4,tmp5,tmp6,tmp7))
    file.write('\\\\\n')
    file.write(r'One sector')
    tmp2 = rmse(baseline_us_1sector.tby[0:17],counter_us_1sector.tby[0:17])
    tmp3 = rmse(baseline_us_1sector.tby[0:17],fix0_us_1sector.tby[0:17])
    tmp4 = rmse(baseline_us_1sector.tby[0:17],fix1_us_1sector.tby[0:17])
    tmp5 = rmse(baseline_us_1sector.tby[0:17],fix2_us_1sector.tby[0:17])
    tmp6 = rmse(baseline_us_1sector.tby[0:17],fix3_us_1sector.tby[0:17])
    tmp7 = rmse(baseline_us_1sector.tby[0:17],fix4_us_1sector.tby[0:17])
    file.write('& %0.2f & %0.2f & %0.2f & %0.2f & %0.2f & %0.2f' % (tmp2,tmp3,tmp4,tmp5,tmp6,tmp7))
    file.write('\\\\\n')
    file.write(r'Low capital adj. costs')
    tmp2 = rmse(baseline_us_etaKh.tby[0:17],counter_us_etaKh.tby[0:17])
    tmp3 = rmse(baseline_us_etaKh.tby[0:17],fix0_us_etaKh.tby[0:17])
    tmp4 = rmse(baseline_us_etaKh.tby[0:17],fix1_us_etaKh.tby[0:17])
    tmp5 = rmse(baseline_us_etaKh.tby[0:17],fix2_us_etaKh.tby[0:17])
    tmp6 = rmse(baseline_us_etaKh.tby[0:17],fix3_us_etaKh.tby[0:17])
    tmp7 = rmse(baseline_us_etaKh.tby[0:17],fix4_us_etaKh.tby[0:17])
    file.write('& %0.2f & %0.2f & %0.2f & %0.2f & %0.2f & %0.2f' % (tmp2,tmp3,tmp4,tmp5,tmp6,tmp7))
    file.write('\\\\\n')
    file.write(r'High capital adj. costs')
    tmp2 = rmse(baseline_us_etaKl.tby[0:17],counter_us_etaKl.tby[0:17])
    tmp3 = rmse(baseline_us_etaKl.tby[0:17],fix0_us_etaKl.tby[0:17])
    tmp4 = rmse(baseline_us_etaKl.tby[0:17],fix1_us_etaKl.tby[0:17])
    tmp5 = rmse(baseline_us_etaKl.tby[0:17],fix2_us_etaKl.tby[0:17])
    tmp6 = rmse(baseline_us_etaKl.tby[0:17],fix3_us_etaKl.tby[0:17])
    tmp7 = rmse(baseline_us_etaKl.tby[0:17],fix4_us_etaKl.tby[0:17])
    file.write('& %0.2f & %0.2f & %0.2f & %0.2f & %0.2f & %0.2f' % (tmp2,tmp3,tmp4,tmp5,tmp6,tmp7))
    file.write('\\\\\n')
    file.write(r'Low Armington elasticities')
    tmp2 = rmse(baseline_us_arml.tby[0:17],counter_us_arml.tby[0:17])
    tmp3 = rmse(baseline_us_arml.tby[0:17],fix0_us_arml.tby[0:17])
    tmp4 = rmse(baseline_us_arml.tby[0:17],fix1_us_arml.tby[0:17])
    tmp5 = rmse(baseline_us_arml.tby[0:17],fix2_us_arml.tby[0:17])
    tmp6 = rmse(baseline_us_arml.tby[0:17],fix3_us_arml.tby[0:17])
    tmp7 = rmse(baseline_us_arml.tby[0:17],fix4_us_arml.tby[0:17])
    file.write('& %0.2f & %0.2f & %0.2f & %0.2f & %0.2f & %0.2f' % (tmp2,tmp3,tmp4,tmp5,tmp6,tmp7))
    file.write('\\\\\n')
    file.write(r'High Armington elasticities')
    tmp2 = rmse(baseline_us_armh.tby[0:17],counter_us_armh.tby[0:17])
    tmp3 = rmse(baseline_us_armh.tby[0:17],fix0_us_armh.tby[0:17])
    tmp4 = rmse(baseline_us_armh.tby[0:17],fix1_us_armh.tby[0:17])
    tmp5 = rmse(baseline_us_armh.tby[0:17],fix2_us_armh.tby[0:17])
    tmp6 = rmse(baseline_us_armh.tby[0:17],fix3_us_armh.tby[0:17])
    tmp7 = rmse(baseline_us_armh.tby[0:17],fix4_us_armh.tby[0:17])
    file.write('& %0.2f & %0.2f & %0.2f & %0.2f & %0.2f & %0.2f' % (tmp2,tmp3,tmp4,tmp5,tmp6,tmp7))
    file.write('\\\\\n')
    file.write(r'Low wedge persistence')
    tmp2 = rmse(baseline_us_rhowl.tby[0:17],counter_us_rhowl.tby[0:17])
    tmp3 = rmse(baseline_us_rhowl.tby[0:17],fix0_us_rhowl.tby[0:17])
    tmp4 = rmse(baseline_us_rhowl.tby[0:17],fix1_us_rhowl.tby[0:17])
    tmp5 = rmse(baseline_us_rhowl.tby[0:17],fix2_us_rhowl.tby[0:17])
    tmp6 = rmse(baseline_us_rhowl.tby[0:17],fix3_us_rhowl.tby[0:17])
    tmp7 = rmse(baseline_us_rhowl.tby[0:17],fix4_us_rhowl.tby[0:17])
    file.write('& %0.2f & %0.2f & %0.2f & %0.2f & %0.2f & %0.2f' % (tmp2,tmp3,tmp4,tmp5,tmp6,tmp7))
    file.write('\\\\\n')
    file.write(r'High wedge persistence')
    tmp2 = rmse(baseline_us_rhowh.tby[0:17],counter_us_rhowh.tby[0:17])
    tmp3 = rmse(baseline_us_rhowh.tby[0:17],fix0_us_rhowh.tby[0:17])
    tmp4 = rmse(baseline_us_rhowh.tby[0:17],fix1_us_rhowh.tby[0:17])
    tmp5 = rmse(baseline_us_rhowh.tby[0:17],fix2_us_rhowh.tby[0:17])
    tmp6 = rmse(baseline_us_rhowh.tby[0:17],fix3_us_rhowh.tby[0:17])
    tmp7 = rmse(baseline_us_rhowh.tby[0:17],fix4_us_rhowh.tby[0:17])
    file.write('& %0.2f & %0.2f & %0.2f & %0.2f & %0.2f & %0.2f' % (tmp2,tmp3,tmp4,tmp5,tmp6,tmp7))
    file.write('\\\\\n')

    file.write(r'Lower long-run real interest rate')
    tmp2 = rmse(baseline_us_lowr.tby[0:17],counter_us_lowr.tby[0:17])
    tmp3 = rmse(baseline_us_lowr.tby[0:17],fix0_us_lowr.tby[0:17])
    tmp4 = rmse(baseline_us_lowr.tby[0:17],fix1_us_lowr.tby[0:17])
    tmp5 = rmse(baseline_us_lowr.tby[0:17],fix2_us_lowr.tby[0:17])
    tmp6 = rmse(baseline_us_lowr.tby[0:17],fix3_us_lowr.tby[0:17])
    tmp7 = rmse(baseline_us_lowr.tby[0:17],fix4_us_lowr.tby[0:17])
    file.write('& %0.2f & %0.2f & %0.2f & %0.2f & %0.2f & %0.2f' % (tmp2,tmp3,tmp4,tmp5,tmp6,tmp7))
    file.write('\\\\\n')

    file.write('\\\\\n')

    file.write('\\multicolumn{7}{l}{\\textit{(b) Fraction of CED explained}}\\\\\n')
    file.write(r'Baseline')
    tmp3 = fex(counter_us.tby[0:17],baseline_us.tby[0:17],fix0_us.tby[0:17])
    tmp4 = fex(counter_us.tby[0:17],baseline_us.tby[0:17],fix1_us.tby[0:17])
    tmp5 = fex(counter_us.tby[0:17],baseline_us.tby[0:17],fix2_us.tby[0:17])
    tmp6 = fex(counter_us.tby[0:17],baseline_us.tby[0:17],fix3_us.tby[0:17])
    tmp7 = fex(counter_us.tby[0:17],baseline_us.tby[0:17],fix4_us.tby[0:17])
    file.write('& 0.00 & %0.2f & %0.2f & %0.2f & %0.2f & %0.2f' % (tmp3,tmp4,tmp5,tmp6,tmp7))
    file.write('\\\\\n')
    file.write(r'No input-output linkages')
    tmp3 = fex(counter_us_noio.tby[0:17],baseline_us_noio.tby[0:17],fix0_us_noio.tby[0:17])
    tmp4 = fex(counter_us_noio.tby[0:17],baseline_us_noio.tby[0:17],fix1_us_noio.tby[0:17])
    tmp5 = fex(counter_us_noio.tby[0:17],baseline_us_noio.tby[0:17],fix2_us_noio.tby[0:17])
    tmp6 = fex(counter_us_noio.tby[0:17],baseline_us_noio.tby[0:17],fix3_us_noio.tby[0:17])
    tmp7 = fex(counter_us_noio.tby[0:17],baseline_us_noio.tby[0:17],fix4_us_noio.tby[0:17])
    file.write('& 0.00 & %0.2f & %0.2f & %0.2f & %0.2f & %0.2f' % (tmp3,tmp4,tmp5,tmp6,tmp7))
    file.write('\\\\\n')
    file.write(r'One sector')
    tmp3 = fex(counter_us_1sector.tby[0:17],baseline_us_1sector.tby[0:17],fix0_us_1sector.tby[0:17])
    tmp4 = fex(counter_us_1sector.tby[0:17],baseline_us_1sector.tby[0:17],fix1_us_1sector.tby[0:17])
    tmp5 = fex(counter_us_1sector.tby[0:17],baseline_us_1sector.tby[0:17],fix2_us_1sector.tby[0:17])
    tmp6 = fex(counter_us_1sector.tby[0:17],baseline_us_1sector.tby[0:17],fix3_us_1sector.tby[0:17])
    tmp7 = fex(counter_us_1sector.tby[0:17],baseline_us_1sector.tby[0:17],fix4_us_1sector.tby[0:17])
    file.write('& 0.00 & %0.2f & %0.2f & %0.2f & %0.2f & %0.2f' % (tmp3,tmp4,tmp5,tmp6,tmp7))
    file.write('\\\\\n')
    file.write(r'Low capital adj. costs')
    tmp3 = fex(counter_us_etaKh.tby[0:17],baseline_us_etaKh.tby[0:17],fix0_us.tby[0:17])
    tmp4 = fex(counter_us_etaKh.tby[0:17],baseline_us_etaKh.tby[0:17],fix1_us.tby[0:17])
    tmp5 = fex(counter_us_etaKh.tby[0:17],baseline_us_etaKh.tby[0:17],fix2_us.tby[0:17])
    tmp6 = fex(counter_us_etaKh.tby[0:17],baseline_us_etaKh.tby[0:17],fix3_us.tby[0:17])
    tmp7 = fex(counter_us_etaKh.tby[0:17],baseline_us_etaKh.tby[0:17],fix4_us.tby[0:17])
    file.write('& 0.00 & %0.2f & %0.2f & %0.2f & %0.2f & %0.2f' % (tmp3,tmp4,tmp5,tmp6,tmp7))
    file.write('\\\\\n')
    file.write(r'High capital adj. costs')
    tmp3 = fex(counter_us_etaKl.tby[0:17],baseline_us_etaKl.tby[0:17],fix0_us.tby[0:17])
    tmp4 = fex(counter_us_etaKl.tby[0:17],baseline_us_etaKl.tby[0:17],fix1_us.tby[0:17])
    tmp5 = fex(counter_us_etaKl.tby[0:17],baseline_us_etaKl.tby[0:17],fix2_us.tby[0:17])
    tmp6 = fex(counter_us_etaKl.tby[0:17],baseline_us_etaKl.tby[0:17],fix3_us.tby[0:17])
    tmp7 = fex(counter_us_etaKl.tby[0:17],baseline_us_etaKl.tby[0:17],fix4_us.tby[0:17])
    file.write('& 0.00 & %0.2f & %0.2f & %0.2f & %0.2f & %0.2f' % (tmp3,tmp4,tmp5,tmp6,tmp7))
    file.write('\\\\\n')
    file.write(r'Low Armington elasticities')
    tmp3 = fex(counter_us_arml.tby[0:17],baseline_us_arml.tby[0:17],fix0_us.tby[0:17])
    tmp4 = fex(counter_us_arml.tby[0:17],baseline_us_arml.tby[0:17],fix1_us.tby[0:17])
    tmp5 = fex(counter_us_arml.tby[0:17],baseline_us_arml.tby[0:17],fix2_us.tby[0:17])
    tmp6 = fex(counter_us_arml.tby[0:17],baseline_us_arml.tby[0:17],fix3_us.tby[0:17])
    tmp7 = fex(counter_us_arml.tby[0:17],baseline_us_arml.tby[0:17],fix4_us.tby[0:17])
    file.write('& 0.00 & %0.2f & %0.2f & %0.2f & %0.2f & %0.2f' % (tmp3,tmp4,tmp5,tmp6,tmp7))
    file.write('\\\\\n')
    file.write(r'High Armington elasticities')
    tmp3 = fex(counter_us_armh.tby[0:17],baseline_us_armh.tby[0:17],fix0_us.tby[0:17])
    tmp4 = fex(counter_us_armh.tby[0:17],baseline_us_armh.tby[0:17],fix1_us.tby[0:17])
    tmp5 = fex(counter_us_armh.tby[0:17],baseline_us_armh.tby[0:17],fix2_us.tby[0:17])
    tmp6 = fex(counter_us_armh.tby[0:17],baseline_us_armh.tby[0:17],fix3_us.tby[0:17])
    tmp7 = fex(counter_us_armh.tby[0:17],baseline_us_armh.tby[0:17],fix4_us.tby[0:17])
    file.write('& 0.00 & %0.2f & %0.2f & %0.2f & %0.2f & %0.2f' % (tmp3,tmp4,tmp5,tmp6,tmp7))
    file.write('\\\\\n')
    file.write(r'Low wedge persistence')
    tmp3 = fex(counter_us_rhowl.tby[0:17],baseline_us_rhowl.tby[0:17],fix0_us.tby[0:17])
    tmp4 = fex(counter_us_rhowl.tby[0:17],baseline_us_rhowl.tby[0:17],fix1_us.tby[0:17])
    tmp5 = fex(counter_us_rhowl.tby[0:17],baseline_us_rhowl.tby[0:17],fix2_us.tby[0:17])
    tmp6 = fex(counter_us_rhowl.tby[0:17],baseline_us_rhowl.tby[0:17],fix3_us.tby[0:17])
    tmp7 = fex(counter_us_rhowl.tby[0:17],baseline_us_rhowl.tby[0:17],fix4_us.tby[0:17])
    file.write('& 0.00 & %0.2f & %0.2f & %0.2f & %0.2f & %0.2f' % (tmp3,tmp4,tmp5,tmp6,tmp7))
    file.write('\\\\\n')
    file.write(r'High wedge persistence')
    tmp3 = fex(counter_us_rhowh.tby[0:17],baseline_us_rhowh.tby[0:17],fix0_us.tby[0:17])
    tmp4 = fex(counter_us_rhowh.tby[0:17],baseline_us_rhowh.tby[0:17],fix1_us.tby[0:17])
    tmp5 = fex(counter_us_rhowh.tby[0:17],baseline_us_rhowh.tby[0:17],fix2_us.tby[0:17])
    tmp6 = fex(counter_us_rhowh.tby[0:17],baseline_us_rhowh.tby[0:17],fix3_us.tby[0:17])
    tmp7 = fex(counter_us_rhowh.tby[0:17],baseline_us_rhowh.tby[0:17],fix4_us.tby[0:17])
    file.write('& 0.00 & %0.2f & %0.2f & %0.2f & %0.2f & %0.2f' % (tmp3,tmp4,tmp5,tmp6,tmp7))
    file.write('\\\\\n')

    file.write(r'Lower long-run interest rate')
    tmp3 = fex(counter_us_lowr.tby[0:17],baseline_us_lowr.tby[0:17],fix0_us.tby[0:17])
    tmp4 = fex(counter_us_lowr.tby[0:17],baseline_us_lowr.tby[0:17],fix1_us.tby[0:17])
    tmp5 = fex(counter_us_lowr.tby[0:17],baseline_us_lowr.tby[0:17],fix2_us.tby[0:17])
    tmp6 = fex(counter_us_lowr.tby[0:17],baseline_us_lowr.tby[0:17],fix3_us.tby[0:17])
    tmp7 = fex(counter_us_lowr.tby[0:17],baseline_us_lowr.tby[0:17],fix4_us.tby[0:17])
    file.write('& 0.00 & %0.2f & %0.2f & %0.2f & %0.2f & %0.2f' % (tmp3,tmp4,tmp5,tmp6,tmp7))
    file.write('\\\\\n')

    note='The second column reports counterfactual model outcomes when all wedges are set to one. Columns 3--7 report counterfactual model outcomes with one wedge set to its calibrate value in each period and all other wedges held constant. Fraction of CED explained calculated as (cumulative difference between trade balance in model and no wedge counterfactual) divided by (cumulative difference between trade balance in data and no-wedge counterfactual). The low (high) capital adjustment cost models set $\\varphi$ to 0.9 (0.7). The low (high) Armington elasticity models set $1/(1-\\zeta^1)$ to 2 (6) and $1/(1-\\sigma^1)$ to 1.5 (3). The low (high) wedge persistence models set $\\rho_\\tau$ to 0.5 (0.8). The lower long-run real interest rate model sets the long-run real interest rate to 0.5\%.'

    tex_footer(file,note)
    
    #file.write('\\end{landscape}\n')
    
