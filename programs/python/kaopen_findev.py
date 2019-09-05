##################################################################################
# kopen.py
# this script calculates the average capital account openness (Chinn and Ito, 2006)
# for the countries in the rest of the world. the average is GDP-weighted like the
# productivity and demographic series.
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
weight='trd'

if flag==0:
    regions = ['USA']
    suff='2c'
else:
    regions = ['USA','CHN']
    suff='3c'

colors=['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00']

##################################################################################
# kaopen

kaopen = pd.read_stata('../../data/kaopen_2015.dta').rename(columns={'ccode':'country'})
kaopen['year'] = kaopen.year.astype('int64')
kaopen=kaopen[kaopen.year>=1995]
kaopen=kaopen[kaopen.year<=2011]
kaopen=kaopen[kaopen.kaopen.notnull()]

# merge on weights
weights = pd.read_pickle('output/weights_'+suff+'.pik')
kaopen = pd.merge(left=kaopen,right=weights,how='left',on=['year','country'])
kaopen=kaopen[kaopen.gdp_weight.notnull()]

# weighted average
wavg = lambda x: np.average(x,weights=kaopen.loc[x.index,weight+'_weight'])
row = kaopen.groupby('year')['ka_open'].aggregate(wavg).reset_index()

# save to file
row.to_csv('output/kaopen.txt',sep=',',index=False)

##################################################################################
# beck et al

findev = pd.read_csv('../../data/FinStructure_November_2013.csv').drop('country',axis=1).rename(columns={'cncode':'country'})
findev=findev[findev.year>=1995]
findev=findev[findev.year<=2011]
findev=findev[findev.pcrdbofgdp.notnull()]

usa=findev[findev.country=='USA']
usa=usa[['year','pcrdbofgdp']].reset_index()

# merge on weights
weights = pd.read_pickle('output/weights_'+suff+'.pik')
findev = pd.merge(left=findev,right=weights,how='left',on=['year','country'])
findev = findev[findev.gdp_weight.notnull()]

# weighted average
wavg = lambda x: np.average(x,weights=findev.loc[x.index,weight+'_weight'])
row = findev.groupby('year')['pcrdbofgdp'].aggregate(wavg).reset_index()
row['pcrdbofgdp'] = row.pcrdbofgdp/usa.pcrdbofgdp

# save to file
row.to_csv('output/findev.txt',sep=',',index=False)
