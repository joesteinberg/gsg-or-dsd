import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

##################################################################################
# setup

# load the database
wiod = pd.read_pickle('/home/joe/Datasets/WIOD/stata/wiot_stata_sep12/wiod.pik')

# pick whether we are doing 2-country or 3-country version
flag = 0
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
# aggregate data

# aggregate intermediate inputs by use region/sector and source region/sector...
# mask for...
mask = wiod.col_use=='M' # ...intermediate use
mask = np.logical_and(mask,wiod.col_region.isin(regions+['ROW'])) # ... valid use region
mask = np.logical_and(mask,wiod.row_sector>=0)  # ...valid source sector
mask = np.logical_and(mask,wiod.row_region.isin(regions+['ROW'])) # ... valid source region
# group by and aggregate
g = wiod[mask].groupby(['year','col_region','col_sector','row_region','row_sector'])
m = g['value'].sum().reset_index()
m.rename(columns={'value':'M'},inplace=True)

# aggregate consumption by use region and source region/sector...
# mask for...
mask = wiod.col_use=='C' # ...consumption
mask = np.logical_and(mask,wiod.col_region.isin(regions+['ROW'])) # ... valid use region
mask = np.logical_and(mask,wiod.row_sector>=0)  # ...valid source sector
mask = np.logical_and(mask,wiod.row_region.isin(regions+['ROW'])) # ... valid source region
# group by and aggregate
g = wiod[mask].groupby(['year','col_region','row_region','row_sector'])
c = g['value'].sum().reset_index()
c.rename(columns={'value':'C'},inplace=True)

# aggregate investment by use region and source region/sector...
# mask for...
mask = wiod.col_use=='I' # ...investment
mask = np.logical_and(mask,wiod.col_region.isin(regions+['ROW'])) # ... valid use region
mask = np.logical_and(mask,wiod.row_sector>=0)  # ...valid source sector
mask = np.logical_and(mask,wiod.row_region.isin(regions+['ROW'])) # ... valid source region
# group by and aggregate
g = wiod[mask].groupby(['year','col_region','row_region','row_sector'])
x = g['value'].sum().reset_index()
x.rename(columns={'value':'I'},inplace=True)

# aggregate VA by use region/sector
# mask for...
mask = wiod.col_use=='M'
mask = np.logical_and(mask,wiod.col_region.isin(regions+['ROW'])) # ... valid use region
mask = np.logical_and(mask,wiod.row_item.isin([99,61,62,63,64,65]))  # ...valid source sector
mask = np.logical_and(mask,wiod.row_region=='') # ... no assigned region
# group by and aggregate
g = wiod[mask].groupby(['year','col_region','col_sector'])
va = g['value'].sum().reset_index()
va.rename(columns={'value':'VA'},inplace=True)

# aggregate GO by use region/sector
# mask for...
mask = wiod.col_use=='M'
mask = np.logical_and(mask,wiod.col_region.isin(regions+['ROW'])) # ... valid use region
mask = np.logical_and(mask,wiod.row_item==69)  # ...valid source sector
mask = np.logical_and(mask,wiod.row_region=='') # ... no assigned region
# group by and aggregate
g = wiod[mask].groupby(['year','col_region','col_sector'])
go = g['value'].sum().reset_index()
go.rename(columns={'value':'GO'},inplace=True)    

##################################################################################
# merge and check consistency

# merge aggregations with same dimensionality
f = pd.merge(left=c,right=x,how='left',on=['year','col_region','row_region','row_sector'])
vg = pd.merge(left=va,right=go,how='left',on=['year','col_region','col_sector'])

# market clearing
msums = m.groupby(['year','row_region','row_sector'])['M'].sum().reset_index()
msums.rename(columns={'row_region':'region','row_sector':'sector'},inplace=True)

fsums = f.groupby(['year','row_region','row_sector'])[['C','I']].sum().reset_index()
fsums.rename(columns={'row_region':'region','row_sector':'sector'},inplace=True)

gsums = vg[['year','col_region','col_sector','GO']]
gsums = gsums.rename(columns={'col_region':'region','col_sector':'sector'})

sums = pd.merge(left=msums,right=fsums,how='left',on=['year','region','sector'])
sums = pd.merge(left=sums,right=gsums,how='left',on=['year','region','sector'])
sums['diff'] = (sums.GO - sums.M - sums.C - sums.I)/sums.GO
test = sum(sums['diff']>1e-6)

if test>0:
    print 'Market clearing failure!'

##################################################################################
# calculate trade vars

# intermediate trade
m_trd =  m.groupby(['year','col_region','row_region','row_sector'])['M'].sum().reset_index()
m_trd=m_trd[m_trd.col_region != m_trd.row_region]
im_m = m_trd.rename(columns={'col_region':'region','row_sector':'sector','row_region':'partner','M':'im_M'})
ex_m = m_trd.rename(columns={'row_region':'region','row_sector':'sector','col_region':'partner','M':'ex_M'})
m_trd2 = pd.merge(left=ex_m,right=im_m,how='left',on=['year','region','partner','sector'])

# final trade
f_trd = f[f.col_region != f.row_region]
im_f = f_trd.rename(columns={'col_region':'region','row_sector':'sector','row_region':'partner','C':'im_C','I':'im_I'})
ex_f = f_trd.rename(columns={'row_region':'region','row_sector':'sector','col_region':'partner','C':'ex_C','I':'ex_I'})
f_trd2 = pd.merge(left=ex_f,right=im_f,how='left',on=['year','region','partner','sector'])

# merge and calculate totals + balances
trd = pd.merge(left=m_trd2,right=f_trd2,how='left',on=['year','region','partner','sector'])

for d in ['ex','im']:
    trd[d+'_F'] = trd[d+'_C']+trd[d+'_I']
    trd[d] = trd[d+'_M']+trd[d+'_F']

for u  in ['_M','_C','_I','_F','']:
    trd['tb'+u] = trd['ex'+u] - trd['im'+u]

# aggregate by sector and append
cols = []
for d in ['ex','im','tb']:
    for u in ['','_M','_F','_C','_I']:
        cols.append(d+u)

g = trd.groupby(['year','region','partner'])
sums = g[cols].sum().reset_index()
sums['sector'] = 'T'
trd = trd.append(sums)

# aggregate by country and append
g = trd.groupby(['year','region','sector'])
sums = g[cols].sum().reset_index()
sums['partner'] = 'TOT'
trd = trd.append(sums)
trd = trd.sort_values(['year','region','partner','sector']).reset_index(drop=True)

# merge on gdp
gdp = vg.groupby(['year','col_region'])['VA'].sum().reset_index()
gdp.rename(columns={'col_region':'region','VA':'GDP'},inplace=True)

# merge on GDP and calculate fractions
trd = pd.merge(left=trd,right=gdp,how='left',on=['year','region'])

if(flag==0):
    trd = trd[trd.partner=='TOT']

##################################################################################
# calculate aggregates (real GDP + composition of NGDP, trade balances)

# aggregates
gdp = vg.groupby(['year','col_region'])['VA'].sum().reset_index()
cons_inv = f.groupby(['year','col_region'])[['C','I']].sum().reset_index()
aggs = pd.merge(left=gdp,right=cons_inv,how='left',on=['year','col_region'])
aggs.rename(columns={'col_region':'region'},inplace=True)
aggs['C_frac'] = 100*aggs.C/aggs.VA
aggs['I_frac'] = 100*aggs.I/aggs.VA

tb_T = trd[['year','region','tb','tb_M','tb_F']][np.logical_and(trd.sector=='T',trd.partner=='TOT')]
tb_G = trd[['year','region','tb','tb_M','tb_F']][np.logical_and(trd.sector==0,trd.partner=='TOT')]
tb_S = trd[['year','region','tb','tb_M','tb_F']][np.logical_and(trd.sector==1,trd.partner=='TOT')]
tb_G.rename(columns={'tb':'tb_G','tb_M':'tb_M_G','tb_F':'tb_F_G'},inplace=True)
tb_S.rename(columns={'tb':'tb_S','tb_M':'tb_M_S','tb_F':'tb_F_S'},inplace=True)
aggs = pd.merge(left=aggs,right=tb_T,how='left',on=['year','region'])
aggs = pd.merge(left=aggs,right=tb_G,how='left',on=['year','region'])
aggs = pd.merge(left=aggs,right=tb_S,how='left',on=['year','region'])
aggs['TB_frac'] = 100*aggs.tb/aggs.VA
aggs['TB_M_frac'] = 100*aggs.tb_M/aggs.VA
aggs['TB_F_frac'] = 100*aggs.tb_F/aggs.VA
aggs['TB_G_frac'] = 100*aggs.tb_G/aggs.VA
aggs['TB_S_frac'] = 100*aggs.tb_S/aggs.VA
aggs['TB_M_G_frac'] = 100*aggs.tb_M_G/aggs.VA
aggs['TB_F_G_frac'] = 100*aggs.tb_F_G/aggs.VA
aggs['TB_M_S_frac'] = 100*aggs.tb_M_S/aggs.VA
aggs['TB_F_S_frac'] = 100*aggs.tb_F_S/aggs.VA

# merge on inflation
prices = pd.read_csv('forpython.csv')
aggs = pd.merge(left=aggs,right=prices[['year','cpi']],how='left',on='year')
usgdp0 = aggs['VA'][np.logical_and(aggs.region=='USA',aggs.year==1995)].values[0]
aggs['RGDP'] = (100*aggs['VA']/usgdp0)/(aggs['cpi']/100.0)

# merge on aggregate RER data
prices['region']='USA'
aggs = pd.merge(left=aggs,right=prices[['year','region','rer']],how='left',on=['year','region'])

# grab US and ROW data
aggs_us = aggs[aggs.region=='USA'].reset_index(drop=True)

# write the us data to a pickle file for python plots, tables, etc.
aggs_us.to_pickle('output/key_us_data.pik')

# save some text files for the C program
cols=['RGDP','C_frac','I_frac','TB_frac','TB_M_frac','TB_F_frac','TB_G_frac','TB_S_frac','TB_M_G_frac','TB_F_G_frac','TB_M_S_frac','TB_F_S_frac','rer']
for c in cols:
    aggs_us[c].to_csv('output/'+c+'.txt',index=False)

if flag==0:
    aggs_rw = aggs[aggs.region=='ROW']
    for c in cols[:3]:
        aggs_rw[c].to_csv('output/'+c+'_ROW.txt',index=False)
else:
    aggs_ch = aggs[aggs.region=='CHN']
    for c in cols[:3]:
        aggs_ch[c].to_csv('output/'+c+'_CHN.txt',index=False)

    aggs_rw2 = aggs[aggs.region=='ROW']
    for c in cols[:3]:
        aggs_rw2[c].to_csv('output/'+c+'_ROW2.txt',index=False)

##################################################################################
# save aggregated data

m.to_pickle('output/wiod_m_'+suff+'.pik')
f.to_pickle('output/wiod_f_'+suff+'.pik')
vg.to_pickle('output/wiod_vg_'+suff+'.pik')
trd.to_pickle('output/wiod_trd_'+suff+'.pik')
