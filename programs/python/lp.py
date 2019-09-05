##################################################################################
# demo.py
# this script calculates historical and projected demographic series for
# the model countries using UN World Population Projects data. the series
# for the rest of the world are weighted averages, with weights coming from
# the wiod_weights.py script

##################################################################################
# setup

# imports
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rc('font',**{'family':'serif','serif':['Palatino'],'size':20})
mpl.rc('font',size=20)
mpl.rc('text', usetex=True)
mpl.rc('lines',linewidth=1.5)
mpl.rc('savefig',bbox='tight')
mpl.rc('savefig',format='pdf')

# pick whether we are doing 2-country or 3-country version
flag = 0
regions=[]
suff = ''
weight='gdp'

if flag==0:
    regions = ['USA']
    suff='2c'
else:
    regions = ['USA','CHN']
    suff='3c'

colors=['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00']

##################################################################################
# load the data files

ted = pd.read_csv('../../data/ted.csv')

# shift rgdppc by one period
#ted['rgdppc_2015_usd'] = ted.groupby('country')['rgdppc_2015usd'].transform(lambda v: v.shift())

weights = pd.read_pickle('output/weights_'+suff+'.pik')

ted['lp_hr_2015usd'] = ted['lp_hr_2015usd'].astype(float)
ted['lp_person_2015usd'] = ted['lp_person_2015usd'].astype(float)
ted.loc[ted.country=='CHN1','country'] = 'CHN'

ted = pd.merge(left=ted,right=weights,how='left',on=['year','country'])
mask = ted.country.isin(regions)
mask = np.logical_or(mask,ted.gdp_weight.notnull())
ted = ted[mask]
ted = ted[ted.year>=1995]
ted = ted[ted.year<=2011]
ted['gdp_weight'] = ted['gdp_weight'].fillna(0.0)
ted['trd_weight'] = ted['trd_weight'].fillna(0.0)

##################################################################################
# calculate LP series

ted['lp'] = ted['lp_hr_2015usd']
tmp = ted.country[ted.lp.isnull()].unique()
ted.loc[ted.country.isin(tmp),'lp'] = ted['lp_person_2015usd'][ted.country.isin(tmp)]

# group by country and normalize
ted = ted.sort_values(['country','year'],ascending=[1,1]).reset_index()
ted['lp_n'] = ted.groupby('country')['lp'].transform(lambda v: v/v.iloc[0])
ted['rgdppc_n'] = ted.groupby('country')['rgdppc_2015usd'].transform(lambda v: v/v.iloc[0])

# pull out USA (and CHN if flag==1) data
usa = ted[['year','lp_n','rgdppc_n']][ted.country=='USA']
usa = usa.rename(columns={'lp_n':'lp_usa','rgdppc_n':'rgdppc_usa'})

chn=0
if flag==1:
    chn = ted[['year','lp_n','rgdppc_n']][ted.country=='CHN']
    chn = chn.rename(columns={'lp_n':'lp_chn','rgdppc_n':'rgdppc_chn'})

# weighted average of other countries
wavg = lambda x: np.average(x,weights=ted.loc[x.index,weight+'_weight'])
row = ted.groupby('year')['lp_n','rgdppc_n'].aggregate(wavg).reset_index()
row = row.rename(columns={'lp_n':'lp_row','rgdppc_n':'rgdppc_row'})

# save text files for program input
if flag==0: # no need to do this twice
    usa['lp_usa'].to_csv('output/lp_USA.txt',sep=' ',header=False,index=False)
    usa['rgdppc_usa'].to_csv('output/rgdppc_USA.txt',sep=' ',header=False,index=False)
    row['lp_row'].to_csv('output/lp_ROW_'+weight+'.txt',sep=' ',header=False,index=False)
    row['rgdppc_row'].to_csv('output/rgdppc_ROW_'+weight+'.txt',sep=' ',header=False,index=False)
else:
    chn['lp_chn'].to_csv('output/lp_CHN.txt',sep=' ',header=False,index=False)
    chn['rgdppc_chn'].to_csv('output/rgdppc_CHN.txt',sep=' ',header=False,index=False)
    row['lp_row'].to_csv('output/lp_ROW2_'+weight+'.txt',sep=' ',header=False,index=False)
    row['rgdppc_row'].to_csv('output/rgdppc_ROW2_'+weight+'.txt',sep=' ',header=False,index=False)

# save pickle files for later use in tables and plots
lp_agg = pd.merge(left=usa,right=row,how='left',on='year')
if(flag==1):
    lp_agg = pd.merge(left=lp_agg,right=chn,how='left',on='year')

lp_agg.to_pickle('output/lp_'+suff+'_'+weight+'.pik')

# display average annual growth rates (for use in projections)
first = lambda x: x.iloc[0]
last = lambda x: x.iloc[-1]
avg_chg = lambda x: 100*((x.iloc[-1]/x.iloc[0])**(1.0/(x.count()-1))-1)

cols = []
if flag==0:
    cols = ['lp_usa','rgdppc_usa','lp_row','rgdppc_row']
else:
    cols = ['lp_usa','rgdppc_usa','lp_chn','rgdppc_chn','lp_row','rgdppc_row']

avgs = lp_agg.groupby(lambda _ : True)[cols].aggregate(avg_chg).reset_index(drop=True)

print weight + '-weighted avg LP growth rates, 1995-2011'
print avgs
