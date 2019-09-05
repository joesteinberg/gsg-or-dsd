import sys
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import MultipleLocator

mpl.rc('text', usetex=True)
mpl.rc('savefig',bbox='tight')
mpl.rc('savefig',format='pdf')

alpha=0.8
colors=['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#a65628']
dashes=[(16,8),(5,3)]
patches = [mpatches.Patch(color=c) for c in colors]


# ------------------------------------------------------------------------------------
# load the data and model results

suff = ''

# data
aggs = pd.read_pickle('output/key_us_data.pik').reset_index(drop=True)

# baseline (where we solve for wedges)
baseline = pd.read_csv('../c/output/vars1_2c_usa'+suff+'.csv')
baseline_row = pd.read_csv('../c/output/vars1_2c_row'+suff+'.csv')

# counterfactual without wedges
counter = pd.read_csv('../c/output/vars0_2c_usa'+suff+'.csv')

# fix US saving wedge
fix0 = pd.read_csv('../c/output/vars2_0_2c_usa'+suff+'.csv')

# fix US investment wedge
fix1 = pd.read_csv('../c/output/vars2_1_2c_usa'+suff+'.csv')

# fix RW saving wedge
fix2 = pd.read_csv('../c/output/vars2_2_2c_usa'+suff+'.csv')

# fix RW inv wedge
fix3 = pd.read_csv('../c/output/vars2_3_2c_usa'+suff+'.csv')

# fix RW trd wedge
fix4 = pd.read_csv('../c/output/vars2_4_2c_usa'+suff+'.csv')

fix = [fix0,fix2,fix1,fix3,fix4]

years = baseline.period + 1995
maxyr = 2025

# ------------------------------------------------------------------------------------
# load the alternative model results

suffs = ['_symg','_nodemo','_notrdw','_noio','_1sector']

base_sens = [baseline]
counter_sens = [counter]
fix_sens = [fix]

for suff in suffs:
    baseline_ = pd.read_csv('../c/output/vars1_2c_usa'+suff+'.csv')
    counter_ = pd.read_csv('../c/output/vars0_2c_usa'+suff+'.csv')
    fix0_ = pd.read_csv('../c/output/vars2_0_2c_usa'+suff+'.csv')
    fix1_ = pd.read_csv('../c/output/vars2_1_2c_usa'+suff+'.csv')
    fix2_ = pd.read_csv('../c/output/vars2_2_2c_usa'+suff+'.csv')
    fix3_ = pd.read_csv('../c/output/vars2_3_2c_usa'+suff+'.csv')
    if(suff!='_notrdw'):
        fix4_ = pd.read_csv('../c/output/vars2_4_2c_usa'+suff+'.csv')
    else:
        fix4_ = 0
    fix_ = [fix0_,fix1_,fix2_,fix3_,fix4_]
    base_sens.append(baseline_)
    counter_sens.append(counter_)
    fix_sens.append(fix_)


# ------------------------------------------------------------------------------------
# data plot: trade balance, real exchange rate, real interest rate, investment rate
mpl.rc('font',**{'family':'serif','serif':['Palatino'],'size':10})
mpl.rc('font',size=10)
mpl.rc('lines',linewidth=1.5)

fig,axes=plt.subplots(2,2,figsize=(6.5,6.5),sharex=False,sharey=False)

axes[0,0].plot(years,baseline['tby'],color=colors[0])
axes[0,0].set_title('(a) Trade balance (pct. GDP)',y=1.04,size=10)
axes[0,0].set_ylim(-6,0)

axes[0,1].plot(years,baseline['rir'],color=colors[1])
axes[0,1].set_title('(b) Real interest rate (pct. per year)',y=1.04,size=10)
axes[0,1].set_ylim(0.5,4.5)

axes[1,0].plot(years,baseline['reer'],color=colors[2])
axes[1,0].set_title('(c) Real exchange rate (1995 = 100)',y=1.04,size=10)
axes[1,0].set_ylim(75,110)

axes[1,1].plot(years,baseline['iy'],color=colors[3])
axes[1,1].set_title('(d) Investment (pct. GDP)',y=1.04,size=10)

ml=MultipleLocator(1)
row=0
col=0
for i in range(4):
    axes[row,col].set_xlim(1995,2011)
    axes[row,col].set_xticks([1998,2002,2006,2010])
    axes[row,col].xaxis.set_minor_locator(ml)
    col = col+1
    if col==2:
        col=0
        row=row+1

# save and close
axes[1,0].text(1995,71,"Notes: The data source for panels (a) and (d) is the World Input Output Database\n(WIOD). In panel (b), the real interest rate is calculated using the yield on 10-year\nTreasury bonds and the realized rate of CPI-U inflation. The source for panel (c) is\nthe International Monetary Fund's International Financial Statistics (IFS) database.",fontsize=10,va='top',ha='left')
fig.subplots_adjust(hspace=0.3,wspace=0.2)
plt.savefig('output/main_data.pdf',bbox='tight')
plt.clf()
plt.close()

# ------------------------------------------------------------------------------------
# model fit plot 1: trade balance disaggregation

mpl.rc('font',**{'family':'serif','serif':['Palatino'],'size':10})
mpl.rc('font',size=10)
mpl.rc('lines',linewidth=1.5)

# plot the data
fig=plt.subplots(1,1,figsize=(5,3.5))
plt.plot(aggs.year,aggs.TB_M_G_frac,color=colors[0],linestyle='-',marker='o',markeredgewidth=0,alpha=alpha,markersize=5)
plt.plot(aggs.year,aggs.TB_F_G_frac,color=colors[1],linestyle='-',marker='o',markeredgewidth=0,alpha=alpha,markersize=5)
plt.plot(aggs.year,aggs.TB_M_S_frac,color=colors[2],linestyle='-',marker='o',markeredgewidth=0,alpha=alpha,markersize=5)
plt.plot(aggs.year,aggs.TB_F_S_frac,color=colors[3],linestyle='-',marker='o',markeredgewidth=0,alpha=alpha,markersize=5)

# plot the model
plt.plot(years,baseline['tbsm0-1'],color=colors[0],linestyle='--',alpha=alpha)
plt.plot(years,baseline['tbsf0-1'],color=colors[1],linestyle='--',alpha=alpha)
plt.plot(years,baseline['tbsm1-1'],color=colors[2],linestyle='--',alpha=alpha)
plt.plot(years,baseline['tbsf1-1'],color=colors[3],linestyle='--',alpha=alpha)

# labels and footnote
plt.text(2013.75,-0.6,r'Intermediate goods',fontsize=10)
plt.text(2015,-1.65,r'Final goods',fontsize=10)
plt.text(2010,2.35,r'Intermediate services',fontsize=10)
plt.text(2010,0.5,r'Final services',fontsize=10)
plt.text(1995,-4.6,'Notes: Colors denote trade categories (blue = final goods, red = \nintermediate goods,purple = final services, green = intermediate\nservices). Data are represented by solid lines with round markers.\nModel results are represented by dashed lines.',fontsize=10,va='top',ha='left')
plt.ylabel('percent GDP')

# adjust axes
plt.xlim(1995,maxyr)
plt.ylim(-4,3)

# save and close
plt.savefig('output/tb_model_vs_data.pdf',bbox='tight')
plt.clf()
plt.close()

# ------------------------------------------------------------------------------------
# plot the wedges

mpl.rc('font',**{'family':'serif','serif':['Palatino'],'size':10})
mpl.rc('font',size=10)
mpl.rc('lines',linewidth=1.5)

cols_us = ['tau_s','tau_i']
cols_rw = ['tau_s','tau_i','tau_t']
fig, axes = plt.subplots(3,2,figsize=(6.5,9),sharey='row',sharex=False)

#baseline['tau_i'] = baseline.tau_i*(1.0-0.39)
baseline_row['tau_i'] = baseline_row.tau_i*0.4+0.6
plot_data=[baseline.tau_s,baseline_row.tau_s,baseline.tau_i,baseline_row.tau_i,baseline_row.tau_t]

row=0
col=0
for i in range(5):
    axes[row,col].plot(years,plot_data[i],
                 color=colors[i],
                 linestyle='-',
                 alpha=0.8)
    col = col+1
    if col==2:
        col=0
        row=row+1
    
axes[0,0].set_title('(a) US saving wedge',y=1.04,size=10)
axes[0,1].set_title('(b) RW saving wedge',y=1.04,size=10)
axes[1,0].set_title('(c) US inv. wedge',y=1.04,size=10)
axes[1,1].set_title('(d) RW inv. wedge',y=1.04,size=10)
axes[2,0].set_title('(e) Trade wedge',y=1.04,size=10)

axes[0,0].set_ylim(0.8,1.2)
axes[1,0].set_ylim(0.4,1.7)
axes[2,0].set_ylim(0.4,1.7)

ml=MultipleLocator(1)
row=0
col=0
for i in range(5):
    axes[row,col].set_xlim(1995,maxyr)
    axes[row,col].set_xticks([2000,2010,2020])

    axes[row,col].xaxis.set_minor_locator(ml)
    col = col+1
    if col==2:
        col=0
        row=row+1

axes[2,1].set_xticks([])
#axes[2,1].set_yticks([])

fig.subplots_adjust(hspace=0.25,wspace=0.1)

#axes[2,1].get_xaxis().set_visible(False)
#axes[2,1].get_yaxis().set_visible(False)
axes[2,1].axis('off')

#axes[2,0].text(1996,0.2,"Notes: Each p,fontsize=10,va='top',ha='left')

# save and close
plt.savefig('output/wedges.pdf',bbox='tight')
plt.clf()
plt.close()


# ------------------------------------------------------------------------------------
# baseline vs. symg and nodemo models

mpl.rc('font',**{'family':'serif','serif':['Palatino'],'size':10})
mpl.rc('font',size=10)
mpl.rc('lines',linewidth=1.5)

fig, axes = plt.subplots(1,2,sharey=True,figsize=(7,3.5))

axes[0].plot(years,baseline['tby'],color=colors[0],linestyle='-',alpha=0.8,zorder=0,marker='o',markeredgewidth=0,markersize=5)
axes[0].plot(years,counter['tby'],color=colors[1],linestyle='-',alpha=0.8,zorder=1)
axes[0].plot(years,counter_sens[1]['tby'],color=colors[2],linestyle='--',dashes=dashes[0],alpha=0.8,zorder=2)
axes[0].plot(years,counter_sens[2]['tby'],color=colors[3],linestyle='--',dashes=dashes[1],alpha=0.8,zorder=3)
axes[0].set_xlim(1995,maxyr)
axes[0].set_title('(a) No-wedge counterfactuals',y=1.04,size=10)
axes[0].text(2000,6.65,r'No demographics',fontsize=10)
axes[0].text(1997,4.2,r'Baseline',size=10)
axes[0].text(2000,-0.5,r'Sym. growth',size=10)
axes[0].text(2009,-5,r'Data + all wedges',size=10)


axes[1].plot(years,baseline['tby'],color=colors[0],linestyle='-',alpha=0.8,zorder=0,marker='o',markeredgewidth=0,markersize=5)
axes[1].plot(years,fix[2]['tby'],color=colors[1],linestyle='-',alpha=0.8,zorder=1)
axes[1].plot(years,fix_sens[1][2]['tby'],color=colors[2],linestyle='--',dashes=dashes[0],alpha=0.8,zorder=2)
axes[1].plot(years,fix_sens[2][2]['tby'],color=colors[3],linestyle='--',dashes=dashes[1],alpha=0.8,zorder=3)
axes[1].set_xlim(1995,maxyr)
axes[1].set_title('(b) RW saving wedge only',y=1.04,size=10)
axes[1].text(2011,3,r'No demographics',fontsize=10)
axes[1].annotate('Baseline',size=10,xytext=(1998,1),xy=(2005,-3),xycoords='data',textcoords='data',arrowprops=dict(arrowstyle='->'))
axes[1].text(2014,-2.5,r'Sym. growth',size=10)
axes[1].text(2009,-5,r'Data + all wedges',size=10)

axes[0].set_ylim(-6,8)
axes[0].set_ylabel('percent GDP',size=10)

axes[0].text(1995,-7.5,"Notes: Both panels plot the observed U.S. trade balance (and future projections from the\nmodel with all wedges) in red with round markers. Panel (a) plots the U.S. trade balance\nin the no-wedge counterfactual for three versions of the analysis: baseline (solid line), sy-\nmmetric productivity growth (long dashes), and without demographics (short dashes).\nPanel (b) plots the counterfactual outcomes implied by the rest of the world's saving we-\ndge using the same line style scheme.",fontsize=10,va='top',ha='left')

ml=MultipleLocator(1)
for i in [0,1]:
    axes[i].set_xlim(1995,maxyr)
    axes[i].set_xticks([2000,2010,2020])
    axes[i].xaxis.set_minor_locator(ml)

fig.subplots_adjust(hspace=0.1,wspace=0.1)

plt.savefig('output/symg_nodemo_sens.pdf',bbox='tight')
plt.clf()
plt.close()

# ------------------------------------------------------------------------------------
# baseline vs. no trade wedge model

mpl.rc('font',**{'family':'serif','serif':['Palatino'],'size':10})
mpl.rc('font',size=10)
mpl.rc('lines',linewidth=1.5)

fig, axes = plt.subplots(1,2,sharey=False,figsize=(7,3.5))

axes[0].plot(years,baseline['reer'],color=colors[0],linestyle='-',alpha=0.8,zorder=0,marker='o',markeredgewidth=0,markersize=5)
axes[0].plot(years,base_sens[3]['reer'],color=colors[4],linestyle='--',dashes=dashes[1],alpha=0.8,zorder=1)
axes[0].plot(years,counter['reer'],color=colors[1],linestyle='-',alpha=0.8,zorder=2)
axes[0].plot(years,fix[2]['reer'],color=colors[2],linestyle='-',alpha=0.8,zorder=3)
axes[0].plot(years,fix_sens[3][2]['reer'],color=colors[3],linestyle='--',dashes=dashes[1],alpha=0.8,zorder=4)
axes[0].set_xlim(1995,maxyr)
axes[0].set_title('(a) US real exchange rate',y=1.04,size=10)
axes[0].set_ylabel('1995 = 100')
axes[0].text(1997,122,r'No wedges',size=10)
axes[0].text(2000,76.5,r'Data',size=10)
axes[0].annotate('RW saving\nwedge\n(baseline)',size=10,xytext=(1999,105),xy=(2001,87),xycoords='data',textcoords='data',arrowprops=dict(arrowstyle='->'))
axes[0].annotate('All wedges\n(no trd. wedge)',size=10,
                 xy=(2008,83),xytext=(2010,73),
                 xycoords='data',textcoords='data',
                 arrowprops=dict(arrowstyle='->'))
axes[0].annotate('RW saving\nwedge (no\ntrd. wedge)',size=10,
                 xy=(2010,90),xytext=(2015,89),
                 xycoords='data',textcoords='data',
                 arrowprops=dict(arrowstyle='->'))

axes[1].plot(years,baseline['tby'],color=colors[0],linestyle='-',alpha=0.8,zorder=0,marker='o',markeredgewidth=0,markersize=5)
axes[1].plot(years,counter['tby'],color=colors[1],linestyle='-',alpha=0.8,zorder=1)
axes[1].plot(years,fix[2]['tby'],color=colors[2],linestyle='-',alpha=0.8,zorder=3)
axes[1].plot(years,fix_sens[3][2]['tby'],color=colors[3],linestyle='--',dashes=dashes[1],alpha=0.8,zorder=4)
axes[1].set_xlim(1995,maxyr)
axes[1].set_title('(b) US trade balance',y=1.04,size=10)
axes[1].set_ylabel('percent GDP')
axes[1].text(1997,4.2,r'No wedges',size=10)
axes[1].text(2009,-5,r'Data + all wedges',size=10)
axes[1].text(2002.5,-1,'RW saving\nwedge\n(baseline)',fontsize=10)
axes[1].annotate('RW saving\nwedge (no\ntrd. wedge)',size=10,xytext=(2014.8,-2),xy=(2010,-2.375),xycoords='data',textcoords='data',arrowprops=dict(arrowstyle='->'))

ml=MultipleLocator(1)
for i in [0,1]:
    axes[i].set_xlim(1995,maxyr)
    axes[i].set_xticks([2000,2010,2020])
    axes[i].xaxis.set_minor_locator(ml)

axes[0].text(1995,65,"Notes: Panel (a) plots the observed U.S. real exchange rate (solid red), the exchange rate\nwith no wedges (solid blue), the exchange rate in the no-trade-wedge version of the anal-\nysis (dashed yellow), the exchange rate implied by the rest of the world's saving wedge\nin the baseline (solid green) and no-trade-wedge (dashed purple) versions of the analysis.\nPanel (b) presents a similar plot for the U.S. trade balance; both versions of the model ma-\ntch the trade balance data exactly, however, so there is one less line.",fontsize=10,va='top',ha='left')

fig.subplots_adjust(hspace=0.5,wspace=0.3)

plt.savefig('output/notrdw_sens.pdf',bbox='tight')
plt.clf()
plt.close()


fig, axes = plt.subplots(1,1,figsize=(7,3.5))
axes.plot(years,baseline['rir'],color=colors[0],linestyle='-',alpha=0.8,zorder=0,marker='o',markeredgewidth=0,markersize=5)
axes.plot(years,counter['rir'],color=colors[1],linestyle='-',alpha=0.8,zorder=1)
axes.plot(years,fix[3]['rir'],color=colors[2],linestyle='-',alpha=0.8,zorder=3)
axes.plot(years,fix_sens[3][3]['rir'],color=colors[3],linestyle='--',dashes=dashes[1],alpha=0.8,zorder=4)
axes.set_xlim(1995,maxyr)
axes.set_title('(c) US real interest rate',y=1.04,size=10)
axes.set_ylabel('percent GDP')
#axes.text(1997,4.2,r'No wedges',size=10)
#axes.text(2009,-5,r'Data + all wedges',size=10)
#axes.text(2002.5,-1,'RW saving\nwedge\n(baseline)',fontsize=10)
#axes.annotate('RW saving\nwedge (no\ntrd. wedge)',size=10,xytext=(2014.8,-2),xy=(2010,-2.375),xycoords='data',textcoords='data',arrowprops=dict(arrowstyle='->'))

plt.savefig('output/notrdw_sens_r.pdf',bbox='tight')
plt.clf()
plt.close()

# ------------------------------------------------------------------------------------
# "panel" style plot with wedges in columns and variables in rows

mpl.rc('font',**{'family':'serif','serif':['Palatino'],'size':10})
mpl.rc('font',size=10)
mpl.rc('lines',linewidth=1.5)

cols = ['tby','reer','iy','rir']
#fig, axes = plt.subplots(4,5,sharey='row',sharex='col',figsize=(9,6.5))

for i in range(4):

    fig, axes = plt.subplots(3,2,sharey=True,sharex=False,figsize=(6.5,9))

    col = cols[i]

    r=0
    c=0

    for j in range(5):

        axes[r,c].plot(years,baseline[col],
                     color=colors[0],
                     linestyle='-',
                     alpha=0.8,
                     marker='o',markersize=5,markeredgewidth=0,
                     zorder=0)

        axes[r,c].plot(years,counter[col],
                     color=colors[1],
                     linestyle='-',
                     alpha=0.8,
                       #dashes=dashes[0],
                     #marker='s',markersize=2,markeredgewidth=0,
                     zorder=1)

        axes[r,c].plot(years,fix[j][col],
                     color=colors[2],
                     linestyle='--',
                     alpha=0.8,
                       dashes=dashes[0],
                     zorder=5)

        c = c+1
        if c==2:
            c=0
            r=r+1
    

    axes[0,0].set_title('(a) US saving wedge',y=1.04,size=10)
    axes[0,1].set_title('(b) RW saving wedge',y=1.04,size=10)
    axes[1,0].set_title('(c) US inv. wedge',y=1.04,size=10)
    axes[1,1].set_title('(d) RW inv. wedge',y=1.04,size=10)
    axes[2,0].set_title('(e) Trade wedge',y=1.04,size=10)

    txty=0
    if col=='tby':
        txty=12
        axes[0,0].set_ylabel('percent GDP')
        axes[1,0].set_ylabel('percent GDP')
        axes[2,0].set_ylabel('percent GDP')
    elif col=='reer':
        txty=140
        axes[0,0].set_ylabel('1995 = 100')
        axes[1,0].set_ylabel('1995 = 100')
        axes[2,0].set_ylabel('1995 = 100')
    elif col=='iy':
        txty=22
        axes[0,0].set_ylabel('percent GDP')
        axes[1,0].set_ylabel('percent GDP')
        axes[2,0].set_ylabel('percent GDP')
    elif col=='rir':
        txty=6
        axes[0,0].set_ylabel('percent per year')
        axes[1,0].set_ylabel('percent per year')
        axes[2,0].set_ylabel('percent per year')



    axes[2,0].text(2028,txty,"Notes: Each panel shows the effect of\none wedge. Solid red lines marked with\ncircles show the data (and projections\nfrom the model with all wedges). Blue\nlines show the no-wedge counterfactual.\nDashed green lines plot counterfactual\noutcomes with the relevant wedge set to\nits calibrated value in each period hold-\ning all other  wedges constant.",fontsize=10,va='top',ha='left')
        #n  .\n Long-dashed \n Short-\n\n


    if col=='tby':
        axes[0,0].set_ylim(-6,12)
    elif col=='reer':
        axes[0,0].set_ylim(70,140)
    elif col=='iy':
        axes[0,0].set_ylim(13,22)
    else:
        axes[0,0].set_ylim(0,6)

    r=0
    c=0
    for i in range(5):
        axes[r,c].set_xlim(1995,maxyr)
        axes[r,c].set_xticks([2000,2010,2020])
        axes[r,c].xaxis.set_minor_locator(ml)
        
        c = c+1
        if c==2:
            c=0
            r=r+1


    #axes[2,1].set_xticks([])
    axes[2,1].axis('off')


    fig.subplots_adjust(hspace=0.25,wspace=0.1)

    plt.savefig('output/wedge_acc_'+col+'.pdf',bbox='tight')
    plt.clf()
    plt.close('all')


# ------------------------------------------------------------------------------------
# demographic and productivity plots
mpl.rc('font',**{'family':'serif','serif':['Palatino'],'size':10})
mpl.rc('font',size=10)
mpl.rc('lines',linewidth=1.5)

demo = pd.read_pickle('output/demo_2c_gdp.pik')
lp = pd.read_pickle('output/lp_2c_gdp.pik')

for t in range(2010,2026):
    i=t-2011+16
    lp.loc[i,'year']=t
    lp.loc[i,'lp_usa']=lp.loc[i-1,'lp_usa']*1.020175
    lp.loc[i,'lp_row']=lp.loc[i-1,'lp_row']*1.0309967
    demo['pop_adeq_usa'] = (demo.pop_wa_usa*0.5608+0.5*(demo.pop_total_usa-demo.pop_wa_usa*0.5608))/1.5608
    demo['pop_adeq_row'] = (demo.pop_wa_row*0.5608+0.5*(demo.pop_total_row-demo.pop_wa_row*0.5608))/1.5608
    demo['pop_adeq_usa'] = demo.pop_adeq_usa/demo.pop_adeq_usa[0]
    demo['pop_adeq_row'] = demo.pop_adeq_row/demo.pop_adeq_row[0]

# make a graph
fig, axes = plt.subplots(1,2,sharey=False,figsize=(7,3.5))

axes[0].plot(demo.year,demo.pop_adeq_usa,color=colors[0],linestyle='-',alpha=0.8)
axes[0].plot(demo.year,demo.pop_wa_usa,color=colors[0],linestyle='--',dashes=dashes[1],alpha=0.8)
axes[0].plot(demo.year,demo.pop_adeq_row,color=colors[1],linestyle='-',alpha=0.8)
axes[0].plot(demo.year,demo.pop_wa_row,color=colors[1],linestyle='--',dashes=dashes[1],alpha=0.8)
axes[0].set_ylabel('1995 = 1')
axes[0].set_title('(a) Demographics',y=1.04,size=10)
axes[0].set_xlim(1995,maxyr)
axes[0].set_ylim(0.95,1.35)
axes[0].text(2002,1.01,'US (adult eq.)')
axes[0].text(2012,1.29,'RW (adult eq.)')
axes[0].text(1998,1.2,'RW\n(working age)')
axes[0].annotate('US (working age)',xytext=(2009,0.97),xy=(2014,1.165),xycoords='data',textcoords='data',arrowprops=dict(arrowstyle='->'))

axes[1].plot(lp.year,lp.lp_usa,color=colors[0],linestyle='-',alpha=0.8)
axes[1].plot(lp.year,lp.lp_row,color=colors[1],linestyle='--',dashes=dashes[1],alpha=0.8)
axes[1].set_title('(b) Labor productivity',y=1.04,size=10)
axes[1].set_xlim(1995,maxyr)
axes[1].set_ylim(0.8,2.6)
axes[1].text(2014,2.0,r'RW')
axes[1].text(2010,1.25,r'US')

axes[0].text(1995,0.91,"Notes: Panel (a) plots the demographic time series parameters. The U.S. demographic da-\nta is shown in blue and the rest of the world's data is shown in red. Solid lines plot adult-\nequivalent populations and dashed lines plot working-age populations. Panel (b) plots\nthe labor productivity parameters. U.S. labor productivity is a solid red line and the rest\nof the world's labor productivity is a dashed blue line.",fontsize=10,va='top',ha='left')

ml=MultipleLocator(1)
for i in [0,1]:
    axes[i].set_xlim(1995,maxyr)
    axes[i].set_xticks([2000,2010,2020])
    axes[i].xaxis.set_minor_locator(ml)

fig.subplots_adjust(hspace=0.2,wspace=0.15)

plt.savefig('output/demo_lp.pdf')
plt.clf()
plt.close('all')












# ------------------------------------------------------------------------------------
# 2 period model graph



mpl.rc('font',**{'family':'serif','serif':['Palatino'],'size':10})
mpl.rc('font',size=10)


y=1.1
A=1.563
alpha=0.3
beta_us=0.628
beta_rw=0.628

def i_us(r):
    return ((1.0+r)/(A*alpha))**(1.0/(alpha-1.0))

def s_us(r):
    k=i_us(r)
    W = y - k + A*k**alpha/(1.0+r)
    c1 = W/(1.0+beta_us)
    s = y-c1
    return s

def tb_us(r):
    return s_us(r) - i_us(r)

def tb_rw(r):
    W = y + y/(1.0+r)
    c1 = W/(1.0+beta_rw)
    return y-c1

def range_brace(x_min, x_max, mid=0.75, 
                beta1=50.0, beta2=100.0, height=1, 
                initial_divisions=11, resolution_factor=1.5):
    # determine x0 adaptively values using second derivitive
    # could be replaced with less snazzy:
    #   x0 = np.arange(0, 0.5, .001)
    x0 = np.array(())
    tmpx = np.linspace(0, 0.5, initial_divisions)
    tmp = beta1**2 * (np.exp(beta1*tmpx)) * (1-np.exp(beta1*tmpx)) / np.power((1+np.exp(beta1*tmpx)),3)
    tmp += beta2**2 * (np.exp(beta2*(tmpx-0.5))) * (1-np.exp(beta2*(tmpx-0.5))) / np.power((1+np.exp(beta2*(tmpx-0.5))),3)
    for i in range(0, len(tmpx)-1):
        t = int(np.ceil(resolution_factor*max(np.abs(tmp[i:i+2]))/float(initial_divisions)))
        x0 = np.append(x0, np.linspace(tmpx[i],tmpx[i+1],t))
    x0 = np.sort(np.unique(x0)) # sort and remove dups
    # half brace using sum of two logistic functions
    y0 = mid*2*((1/(1.+np.exp(-1*beta1*x0)))-0.5)
    y0 += (1-mid)*2*(1/(1.+np.exp(-1*beta2*(x0-0.5))))
    # concat and scale x
    x = np.concatenate((x0, 1-x0[::-1])) * float((x_max-x_min)) + x_min
    y = np.concatenate((y0, y0[::-1])) * float(height)
    return (x,y)

r_vec = np.arange(0.2,1.0,0.02)
s_us_vec = np.array([s_us(r) for r in r_vec])
i_us_vec = np.array([i_us(r) for r in r_vec])
tb_us_vec = np.array([tb_us(r) for r in r_vec])
tb_rw_vec = np.array([tb_rw(r) for r in r_vec])

beta_us = 0.507
s_us_vec2 = np.array([s_us(r) for r in r_vec])
i_us_vec2 = np.array([i_us(r) for r in r_vec])
tb_us_vec2 = np.array([tb_us(r) for r in r_vec])

beta_rw = 0.705
tb_rw_vec2 = np.array([tb_rw(r) for r in r_vec])

fig, axes = plt.subplots(1,2,figsize=(7,3.5))

axes[0].plot(tb_us_vec,r_vec,color=colors[0],marker='None',alpha=0.85)
axes[0].plot(-tb_rw_vec,r_vec,color=colors[1],marker=None,linestyle='--',alpha=0.85)
axes[0].plot(-tb_rw_vec2,r_vec,color=colors[3],linestyle='--',marker='+',alpha=0.85)
axes[0].scatter([0.0],[0.59256],s=20,color='black',marker='o',zorder=4)
axes[0].scatter([-0.03],[0.517],s=20,color='black',marker='o',zorder=5)
axes[0].set_xlabel('trade balance')
axes[0].set_ylabel('interest rate')
axes[0].set_title('(a) Global saving glut',y=1.04,size=10)
axes[0].axvline(0,color='black',alpha=0.3,linestyle='--',dashes=(3,1),linewidth=1)
axes[0].axvline(-0.03,color='black',alpha=0.15,linestyle='--',dashes=(3,1),linewidth=1)
axes[0].axhline(0.592,color='black',alpha=0.6,linestyle=':',linewidth=1)
axes[0].axhline(0.517,color='black',alpha=0.3,linestyle=':',linewidth=1)
axes[0].set_xlim(-0.1,0.1)
axes[0].set_ylim(0.4,0.9)
axes[0].get_xaxis().set_ticks([])
axes[0].get_yaxis().set_ticks([])
axes[0].text(0.003,0.588,r'A')
axes[0].text(-0.098,0.595,r'$R$',fontsize=14)
axes[0].text(-0.026,0.51,r'B')
axes[0].text(-0.098,0.52,r"$R'$",fontsize=14)
axes[0].text(-0.032,0.43,r'U.S. trade deficit',fontsize=10)
axes[0].text(0.05,0.8,r'$tb_{us}$')
axes[0].text(-0.048,0.8,r'$-tb_{rw}$')
axes[0].text(-0.088,0.75,r"$-tb'_{rw}$")

axes[0].annotate('',xy=(-0.075,0.585-0.06),xytext=(-0.075,0.585),
                 xycoords='data',textcoords='data',
                 arrowprops=dict(arrowstyle='->'))

axes[0].annotate('',xy=(-0.038-0.0335,0.725),xytext=(-0.038,0.725),
                 xycoords='data',textcoords='data',
                 arrowprops=dict(arrowstyle='->'))

x0,y0=range_brace(x_min=-0.029,x_max=-0.001,height=0.02)
y0=[y+0.405 for y in y0]
axes[0].plot(x0,y0,color='black',linewidth=1)

axes[1].plot(tb_us_vec,r_vec,color=colors[0],marker='None',alpha=0.85)
axes[1].plot(-tb_rw_vec,r_vec,color=colors[1],marker=None,linestyle='--',alpha=0.85)
axes[1].plot(tb_us_vec2,r_vec,color=colors[2],alpha=0.85,marker='x')
axes[1].scatter([0.0],[0.59256],s=20,color='black',marker='o',zorder=4)
axes[1].scatter([-0.03],[0.714],s=20,color='black',marker='o',zorder=5)
axes[1].set_xlabel('trade balance')
axes[1].set_title('(b) Domestic saving drought',y=1.04,size=10)
axes[1].axvline(0,color='black',alpha=0.3,linestyle='--',dashes=(3,1),linewidth=1)
axes[1].axvline(-0.03,color='black',alpha=0.15,linestyle='--',dashes=(3,1),linewidth=1)
axes[1].axhline(0.592,color='black',alpha=0.6,linestyle=':',linewidth=1)
axes[1].axhline(0.714,color='black',alpha=0.3,linestyle=':',linewidth=1)
axes[1].set_xlim(-0.1,0.1)
axes[1].set_ylim(0.4,0.9)
axes[1].get_xaxis().set_ticks([])
axes[1].get_yaxis().set_ticks([])
axes[1].text(0.003,0.588,r'A')
axes[1].text(-0.098,0.595,r'$R$',fontsize=14)
axes[1].text(-0.026,0.705,r'C')
axes[1].text(-0.098,0.715,r"$R'$",fontsize=14)
axes[1].text(-0.032,0.43,r'U.S. trade deficit',fontsize=10)
axes[1].text(0.05,0.8,r'$tb_{us}$')
axes[1].text(0.017,0.825,r"$tb'_{us}$")
axes[1].text(-0.048,0.8,r'$-tb_{rw}$')

axes[1].annotate('',xy=(-0.08,0.5965+0.11),xytext=(-0.08,0.5965),
                 xycoords='data',textcoords='data',
                 arrowprops=dict(arrowstyle='->'))

axes[1].annotate('',xy=(0.05375-0.059,0.77),xytext=(0.05375,0.77),
                 xycoords='data',textcoords='data',
                 arrowprops=dict(arrowstyle='->'))

x0,y0=range_brace(x_min=-0.029,x_max=-0.001,height=0.02)
y0=[y+0.405 for y in y0]
axes[1].plot(x0,y0,color='black',linewidth=1)

axes[0].text(-0.1,0.34,"Notes: Panel (a) illustrates the effects of an increase in the rest of the world's saving wedge\non the U.S. trade balance and real interest rate. The line $tb_{us}$ denotes the U.S. trade balance\ncurve. $-tb_{rw}$ denotes the rest of the world's initial trade deficit curve. $-tb'_{rw}$ denotes the\nrest of the world's trade deficit curve after its saving wedge rises. Point A denotes the\ninitial equilibrium and point B denotes the equilibrium after the wedge rises. Panel (b)\nillustrates the effects of a decrease in the U.S. saving wedge using similar notation. Point\nC denotes the equilibrim after the U.S. saving wedge falls.",fontsize=10,va='top',ha='left')



fig.subplots_adjust(hspace=0.2,wspace=0.1)

plt.savefig('output/2period.pdf')
plt.clf()
plt.close()


# ------------------------------------------------------------------------------------
# kaopen plot

mpl.rc('font',**{'family':'serif','serif':['Palatino'],'size':10})
mpl.rc('font',size=10)
mpl.rc('lines',linewidth=1.5)

kaopen=pd.read_csv('output/kaopen.txt')
findev=pd.read_csv('output/findev.txt')

fig, axes = plt.subplots(1,2,figsize=(7,3.5))

axes[0].plot(aggs.year,aggs.TB_frac,color=colors[0],linestyle='-')
axes[0].set_ylabel('US trade balance (pct. GDP)')
axes[0].set_title('(a) Capital account openness',y=1.04,size=10)
ax2=axes[0].twinx()
ax2.plot(kaopen.year,kaopen.ka_open,color=colors[1],linestyle='--')
ax2.set_ylabel('RW capital account openness (US = 1)')
ml=MultipleLocator(1)
axes[0].set_xlim(1995,2011)
axes[0].set_xticks([1998,2002,2006,2010])
axes[0].xaxis.set_minor_locator(ml)
ax2.set_xlim(1995,2011)
ax2.set_xticks([1998,2002,2006,2010])
ax2.xaxis.set_minor_locator(ml)
axes[0].text(2001,-5.45,'US trade balance')
axes[0].text(2000,-1.07,'RW capital\naccount openness')

axes[1].plot(aggs.year,aggs.TB_frac,color=colors[0],linestyle='-')
axes[1].set_ylabel('US trade balance (pct. GDP)')
axes[1].set_title('(b) Financial development',y=1.04,size=10)
ax2=axes[1].twinx()
ax2.plot(findev.year,findev.pcrdbofgdp,color=colors[2],linestyle='--')
ax2.set_ylabel('RW Financial development (US = 1)')
axes[1].set_xlim(1995,2011)
axes[1].set_xticks([1998,2002,2006,2010])
axes[1].xaxis.set_minor_locator(ml)
ax2.set_xlim(1995,2011)
ax2.set_xticks([1998,2002,2006,2010])
ax2.xaxis.set_minor_locator(ml)
axes[1].text(1996.5,-0.7,'US trade balance')
axes[1].text(2000.5,-2.4,'RW financial\ndevelopment')

axes[0].text(1995,-6.675,"Notes: Solid red lines in both panels plot the U.S. trade balance as a fraction of U.S. GDP\n(left axis). Dashed lines in panels (a) and (b), respectively, plot the rest of the world's\ncapital account openness and financial development relative to those of the United States\n(right axes). The former is based on the $KAOPEN$ variable from the Chinn and Ito (2006)\ndataset, while the latter is based on $PCRDBOFGDP$ from the Beck et al. (2000) dataset.",fontsize=10,va='top',ha='left')

fig.subplots_adjust(hspace=0.2,wspace=0.5)

# save and close
plt.savefig('output/kaopen_findev.pdf',bbox='tight')
plt.clf()
plt.close()




plt.close('all')
