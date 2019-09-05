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

colors=['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00']

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

##################################################################################
# load the data files

# load the WPP data
wpp1 = pd.read_csv('/home/joe/Datasets/world-population-prospects/WPP2015_INT_F2A_Annual_Population_Indicators.csv')
wpp2 = pd.read_csv('/home/joe/Datasets/world-population-prospects/WPP2015_INT_F2C_Annual_Population_Indicators_DependencyRatios.csv')
wpp1 = wpp1[['LocID','Time','PopTotal']]
wpp2 = wpp2[['LocID','Time','PopTotTDepR2064']]
wpp = pd.merge(left=wpp1,right=wpp2,how='left',on=['LocID','Time'])
wpp = wpp.rename(columns={'PopTotal':'pop_total','PopTotTDepR2064':'dep_ratio','Time':'year'})
wpp = wpp[wpp.year>=1995]

# merge on country code concordance
codes = pd.read_csv('/home/joe/Datasets/world-population-prospects/wpp-codes.csv')
wpp = pd.merge(left=wpp,right=codes,how='left',on='LocID')
wpp = wpp[wpp.country.notnull()]

# merge on weights for historical data
weights = pd.read_pickle('output/weights_'+suff+'.pik')
wpp_hist = pd.merge(left=wpp[wpp.year<=2011],right=weights,how='left',on=['country','year'])

# use 2011 weights in future projections
wpp_proj = pd.merge(left=wpp[wpp.year>2011],right=weights[weights.year==2011].drop('year',axis=1),how='left',on='country')

# append
wpp = wpp_hist.append(wpp_proj)

# get rid of extraneous data
mask = wpp.country.isin(regions)
mask = np.logical_or(mask,wpp.gdp_weight.notnull())
wpp = wpp[['country','year','gdp_weight','trd_weight','pop_total','dep_ratio']][mask]
wpp['gdp_weight'] = wpp['gdp_weight'].fillna(0.0)
wpp['trd_weight'] = wpp['trd_weight'].fillna(0.0)

##################################################################################
# compute demo series

# compute working age pop
wpp['pop_wa'] = wpp.pop_total - (wpp.dep_ratio/100.0)*wpp.pop_total/(1.0+wpp.dep_ratio/100.0)

# group by country and normalize
wpp = wpp.sort_values(['country','year'],ascending=[1,1]).reset_index()
wpp['pop_total_n'] = wpp.groupby('country')['pop_total'].transform(lambda v: v/v.iloc[0])
wpp['pop_wa_n'] = wpp.groupby('country')['pop_wa'].transform(lambda v: v/v.iloc[0])

# pull out USA (and CHN if flag==1) data
usa = wpp[['year','pop_total_n','pop_wa_n']][wpp.country=='USA']
usa = usa.rename(columns={'pop_total_n':'pop_total_usa','pop_wa_n':'pop_wa_usa'})

chn=0
if flag==1:
    chn = wpp[['year','pop_total_n','pop_wa_n']][wpp.country=='CHN']
    chn = chn.rename(columns={'pop_total_n':'pop_total_chn','pop_wa_n':'pop_wa_chn'})

# take average of other countries
wavg = lambda x: np.average(x,weights=wpp.loc[x.index,weight+'_weight'])
row = wpp.groupby('year')[['pop_total_n','pop_wa_n']].aggregate(wavg).reset_index()
row = row.rename(columns={'pop_total_n':'pop_total_row','pop_wa_n':'pop_wa_row'})

# save text files for program input
if flag==0: # no need to do this twice
    usa[['pop_total_usa','pop_wa_usa']].to_csv('output/demo_USA.txt',sep=' ',header=False,index=False)
    row[['pop_total_row','pop_wa_row']].to_csv('output/demo_ROW_'+weight+'.txt',sep=' ',header=False,index=False)
else:
    chn[['pop_total_chn','pop_wa_chn']].to_csv('output/demo_CHN.txt',sep=' ',header=False,index=False)
    row[['pop_total_row','pop_wa_row']].to_csv('output/demo_ROW2_'+weight+'.txt',sep=' ',header=False,index=False)

# save pickle files for later use in tables and plots
demo_agg = pd.merge(left=usa,right=row,how='left',on='year')
if(flag==1):
    demo_agg = pd.merge(left=demo_agg,right=chn,how='left',on='year')

demo_agg.to_pickle('output/demo_'+suff+'_'+weight+'.pik')

# display average annual growth rates
first = lambda x: x.iloc[0]
last = lambda x: x.iloc[-1]
avg_chg = lambda x: 100*((x.iloc[-1]/x.iloc[0])**(1.0/(x.count()-1))-1)

cols = []
if flag==0:
    cols = ['pop_total_usa','pop_wa_usa','pop_total_row','pop_wa_row']
else:
    cols = ['pop_total_usa','pop_wa_usa','pop_total_chn','pop_wa_chn','pop_total_row','pop_wa_row']

avgs = demo_agg[demo_agg.year<=2011].groupby(lambda _ : True)[cols].aggregate(avg_chg).reset_index(drop=True)

print weight + '-weighted avg pop growth rates, 1995-2011'
print avgs





