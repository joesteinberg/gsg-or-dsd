##################################################################################
# wiod_weights.py
# this script uses WIOD data to compute weights that other parts of the paper
# use to aggregate or average data from other sources

##################################################################################
# setup

# imports
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

# load the database
wiod = pd.read_pickle('/home/joe/Datasets/WIOD/stata/wiot_stata_sep12/wiod.pik')

# pick whether we are doing 2-country or 3-country version
flag = 1
regions=[]
suff = ''

if flag==0:
    regions = ['USA']
    suff='2c'
else:
    regions = ['USA','CHN']
    suff='3c'

# sector aggregation
sectors={0:range(1,16),1:[16,17]+range(19,36),2:[18]}

# set sector and region codes
def which_sector(x):
    for i in sectors.keys():
        if x in sectors[i]:
            return i
    return -1

def which_region(c):
    if c in ['TOT','TXP','CIF','PUA','PUF','VA','ITM','GO']:
        return ''
    else:
        for i in regions:
            if c==i:
                return i
        return 'ROW'

def which_use(u):
    if u in range(1,36):
        return 'M'
    elif u in [37,38,39]:
        return 'C'
    elif u in [41,42]:
        return 'I'

wiod['row_sector']=wiod['row_item'].apply(which_sector)
wiod['col_sector']=wiod['col_item'].apply(which_sector)

wiod['row_region']=wiod['row_country'].apply(which_region)
wiod['col_region']=wiod['col_country'].apply(which_region)

wiod['col_use']=wiod['col_item'].apply(which_use)

##################################################################################
# compute weights

# weight 1: nominal GDP
mask = wiod.col_use=='M'
mask = np.logical_and(mask,wiod.col_region!='USA') # ... don't need the United States here
if flag==1:
    mask = np.logical_and(mask,wiod.col_region!='CHN')
mask = np.logical_and(mask,wiod.col_region.isin(regions+['ROW'])) # ... valid use region
mask = np.logical_and(mask,wiod.row_item.isin([99,61,62,63,64,65]))  # ...valid source sector
mask = np.logical_and(mask,wiod.row_region=='') # ... no assigned region
g = wiod[mask].groupby(['year','col_country','col_region']) # group by and aggregate
gdp = g['value'].sum().reset_index().rename(columns={'value':'gdp',
                                                     'col_country':'country',
                                                     'col_region':'region'})
g = gdp.groupby(['year','region'])
tot = g['gdp'].sum().reset_index().rename(columns={'gdp':'tot_gdp'})
merged_gdp = pd.merge(left=gdp,right=tot,how='left',on=['year','region'])
merged_gdp['gdp_weight'] = merged_gdp.gdp/merged_gdp.tot_gdp
merged_gdp.drop(['gdp','tot_gdp','region'],axis=1,inplace=True)

# weight 2: trade
# have to do exports and imports separately then add them together
mask = wiod.col_region =='USA'
mask = np.logical_and(mask,wiod.row_region !='USA')
if flag==1:
    mask = np.logical_and(mask,wiod.row_region!='CHN')
mask = np.logical_and(mask,wiod.col_use.isin(['M','C','I']))
mask = np.logical_and(mask,wiod.row_sector>=0)
g=wiod[mask].groupby(['year','row_country','row_region'])
im = g['value'].sum().reset_index().rename(columns={'value':'im',
                                                    'row_country':'country',
                                                    'row_region':'region'})

mask = wiod.row_region =='USA'
mask = np.logical_and(mask,wiod.col_region !='USA')
mask = np.logical_and(mask,wiod.col_use.isin(['M','C','I']))
mask = np.logical_and(mask,wiod.row_sector>=0)
g=wiod[mask].groupby(['year','col_country','col_region'])
ex = g['value'].sum().reset_index().rename(columns={'value':'ex',
                                                    'col_country':'country',
                                                    'col_region':'region'})

trd = pd.merge(left=im,right=ex,how='left',on=['year','country','region'])
trd['trd'] = trd.ex+trd.im
g = trd.groupby(['year','region'])
tot = g['trd'].sum().reset_index().rename(columns={'trd':'tot_trd'})
merged_trd = pd.merge(left=trd,right=tot,how='left',on=['year','region'])
merged_trd['trd_weight'] = merged_trd.trd/merged_trd.tot_trd
merged_trd.drop(['ex','im','trd','tot_trd','region'],axis=1,inplace=True)

weights = pd.merge(left=merged_gdp,right=merged_trd,how='left',on=['year','country'])
weights.to_pickle('output/weights_'+suff+'.pik')

