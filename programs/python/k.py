##################################################################################
# k.py
# this script calculates the initial capital stock and average depreciation rates
# for the model countries. values for the rest of the world are weighed averages,
# with weights coming from the wiod_weights.py script

##################################################################################
# setup

# imports
import numpy as np
import pandas as pd

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

pwt = pd.read_stata('/home/joe/Datasets/pwt/pwt90.dta').drop('country',axis=1).rename(columns={'countrycode':'country'})
weights = pd.read_pickle('output/weights_'+suff+'.pik')

# average the weights!
weights2 = weights[['country','gdp_weight','trd_weight']][weights.year==2011]

pwt = pd.merge(left=pwt,right=weights,how='left',on=['country','year'])
#pwt = pd.merge(left=pwt,right=weights2,how='left',on='country')
mask = pwt.country.isin(regions)
mask = np.logical_or(mask,pwt.gdp_weight.notnull())
pwt = pwt[mask]
pwt['gdp_weight'] = pwt['gdp_weight'].fillna(0.0)
pwt['trd_weight'] = pwt['trd_weight'].fillna(0.0)
pwt['k_y_ratio'] = 100*pwt.ck/pwt.cgdpo

##################################################################################
# initial capital stocks

p1995 = pwt[pwt.year==1995]

print 'USA initial K/Y: %0.4f' % p1995.k_y_ratio[p1995.country=='USA']

if flag==1:
    print 'CHN initial K/Y: %0.4f' % p1995.k_y_ratio[p1995.country=='CHN']

wavg = lambda x: np.average(x,weights=p1995.loc[x.index,weight+'_weight'])
print weight + '-weighted ROW initial K/Y: %0.4f' % wavg(p1995.k_y_ratio)

##################################################################################
# depreciation rates

pwt = pwt[pwt.year>=1995]
pwt = pwt[pwt.year<=2011]
dep = pwt.groupby('country')['delta'].mean()
