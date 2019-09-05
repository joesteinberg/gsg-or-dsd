##########################################################################################
# This script does the following:
# - Constructs the IO matrices from processed WIOD data files
# - Writes the matrices to csv files and latex tables
##########################################################################################

##########################################################################################
# Imports, constants, etc.
##########################################################################################
import pandas as pd
import numpy as np
import itertools
import locale
locale.setlocale(locale.LC_ALL,'en_US.utf8')

year=1995
ns=3
sectors=['Goods','Srvcs','Const']

flag = 1
nc=0
regions=[]
suff = ''

if flag==1:
    nc=2
    regions = {'USA':0,'ROW':1}
    countries=['USA','ROW']
    suff='2c'
else:
    nc=3
    regions = {'USA':0,'CHN':1,'ROW':2}
    countries=['USA','CHN','ROW']
    suff='3c'

inpath = 'output/'
outpath = 'output/'

m = pd.read_pickle('output/wiod_m_'+suff+'.pik')
f = pd.read_pickle('output/wiod_f_'+suff+'.pik')
vg = pd.read_pickle('output/wiod_vg_'+suff+'.pik')

def assign_region_num(rstr):
    return regions[rstr]

##########################################################################################
# load the data
##########################################################################################

vg = vg[vg.year==year]
f = f[f.year==year]
m = m[m.year==year]

vg['col_region_num'] = vg['col_region'].apply(assign_region_num)
f['col_region_num'] = f['col_region'].apply(assign_region_num)
f['row_region_num'] = f['row_region'].apply(assign_region_num)
m['col_region_num'] = m['col_region'].apply(assign_region_num)
m['row_region_num'] = m['row_region'].apply(assign_region_num)

# set consumption and intermediate use of construction to zero and make sure construction is completely nontraded
f.loc[f.row_sector==2,'C'] = 0.0
f.loc[np.logical_and(f.row_sector==2,f.row_region_num!=f.col_region_num),'I'] = 0.0
m.loc[m.row_sector==2,'M'] = 0.0

##########################################################################################
# construct the input-output matrix
##########################################################################################

vg = vg[vg.year==year].sort_values(by=['col_region_num','col_sector']).reset_index()
f = f[f.year==year].sort_values(by=['col_region_num','row_region_num','row_sector']).reset_index()
m = m[m.year==year].sort_values(by=['col_region_num','col_sector','row_region_num','row_sector']).reset_index()

vg['col_region_num'] = vg['col_region'].apply(assign_region_num)
f['col_region_num'] = f['col_region'].apply(assign_region_num)
f['row_region_num'] = f['row_region'].apply(assign_region_num)
m['col_region_num'] = m['col_region'].apply(assign_region_num)
m['row_region_num'] = m['row_region'].apply(assign_region_num)

rowsums = np.zeros( nc*ns + 1 )
colsums = np.zeros( nc*ns + nc*2 )

MM = m.pivot_table(values='M', index=['row_region_num','row_sector'], columns=['col_region_num','col_sector'])
VV = vg['VA'].values.reshape((1,nc*ns))
FF = f.pivot_table(values=['C','I'], index=['row_region_num','row_sector'], columns=['col_region_num'])
VV = np.hstack((VV,np.zeros((1,nc*2))))

iomat=np.vstack( ( np.hstack((MM,FF)) , VV ) )

for row in range(0,nc*ns + 1):
    rowsums[row] = np.sum(iomat[row,:])

for col in range(0,nc*ns + nc*2):
    colsums[col] = np.sum(iomat[:,col])

##########################################################################################
# rebalance
##########################################################################################

def coeffs(iomat):
    # Given world IO matrix (iomat), calculates IO coefficients and returs them in A

    A=np.zeros(iomat.shape)
    for col in range(0,A.shape[1]):
        A[:,col] = iomat[:,col]/np.sum(iomat[:,col])
    return A

def ras(iomat0,rowsums1,colsums1):
    # Given an initial IO matrix (iomat), and desired rowsums (rowsums1) and colsums (colsums1),
    # performs the RAS balancing procedure. Returns a new IO matrix (iomat) that is consistent
    # with desired row- and colsums.

    A0 = coeffs(iomat0)
    iomat = np.dot(A0,np.diag(colsums1))

    go=True
    iter=0
    maxit=10000
    tol=1.0e-15

    while go:
        iter=iter+1
        rowsums = np.sum(iomat,axis=1)
        r = np.divide(rowsums1,rowsums)
        iomat = np.dot(np.diag(r),iomat)
        colsums = np.sum(iomat,axis=0)
        s = np.divide(colsums1,colsums)
        iomat = np.dot(iomat,np.diag(s))
        colsums = np.sum(iomat,axis=0)
        rowsums = np.sum(iomat,axis=1)

        norm1 = max(np.divide(abs(rowsums-rowsums1),rowsums1))
        norm2 = max(np.divide(abs(colsums-colsums1),colsums1))
        if((norm1 <tol and norm2 <tol) or iter == maxit):
            go=False

    if iter==maxit:
        print 'RAS iteration did not converge!'
        print 'iter = ', iter, ' diff = ', max(norm1,norm2)
    else:
        print 'RAS converged after ',str(iter),' iterations'


    return iomat

# make sure it's balanced after imposing the restrictions on construction usage...
colsums[0:(nc*ns)] = rowsums[0:(nc*ns)] # make sure markets clear: gross output = total demand for each country/sector
rowsums[-1] = colsums[(nc*ns):].sum() # world value added must equal world final demand
iomat2 = ras(iomat,rowsums,colsums) # run RAS

##########################################################################################
# write output
##########################################################################################

def write_iomat_csv(iomat,fname):
    # Write world IO matrix (iomat) to csv file called filename

    usgdp = iomat[-1,0:3].sum()

    iomat2 = np.vstack((iomat,np.sum(iomat,axis=0).reshape((1,nc*ns+nc*2))))
    iomat2 = np.hstack((iomat2,np.sum(iomat2,axis=1).reshape((nc*ns+2,1))))
    iomat2 = 100*iomat2/usgdp

    np.savetxt(fname=outpath+fname,X=iomat2,fmt='%0.15f',delimiter=' ')

def write_iomat_latex(iomat,rowsums,colsums,caption,label,fname):
    # Given a world IO matrix (iomat), rowsums, colsums, creates a latex file
    # in location filename that contains a table with given caption and label.

    usgdp = iomat[-1,0:3].sum()
    iomat2 = 100*iomat[:,:]/usgdp
    rowsums2 = 100*rowsums/usgdp
    colsums2 = 100*colsums/usgdp

    M=iomat2[0:(nc*ns),0:(nc*ns)]
    V=iomat2[-1,0:(nc*ns)]
    Fc=iomat2[0:(nc*ns),(nc*ns):+((nc*ns)+nc)]
    Fx=iomat2[0:(nc*ns),((nc*ns)+nc):]

    with open(outpath + fname + '.tex','wb') as file:
        #file.write('\\begin{landscape}\n')
        file.write('\\begin{table}[p]\n')
        file.write('\\footnotesize\n')
        file.write('\\begin{center}\n')
        file.write('\\caption{'+caption+'}\n')
        file.write('\\label{tab:'+label+'}\n')

        file.write('\\begin{tabular}{cc')
        for i in range(0,nc*ns+nc*2):
            file.write('c')
        file.write('}\n')
        file.write('\\toprule\n')
            
        # use categories
        file.write('& & \\multicolumn{'+str(nc*ns)+'}{c}{Intermediate inputs}')
        file.write('&\\multicolumn{'+str(nc*2)+'}{c}{Final demand}\\\\\n')
        file.write('\\cmidrule(rl){3-'+str(2+nc*ns)+'}')
        file.write('\\cmidrule(rl){'+str(3+nc*ns)+'-'+str(3+nc*ns+2*nc-1)+'}')

        # country names
        file.write('&')
        for c in countries:
            file.write('& \\multicolumn{'+str(ns)+'}{c}{'+c+'}')
        for c in countries:
            file.write('& \\multicolumn{2}{c}{'+c+'}')
        file.write('\\\\\n')
        
        # underline country names
        x=3
        for i in range(nc):
            file.write('\\cmidrule(rl){'+str(x)+'-'+str(x+ns-1)+'}')
            x=x+ns
        for i in range(nc):
            file.write('\\cmidrule(rl){'+str(x)+'-'+str(x+1)+'}')
            x=x+2
        file.write('\n')

        # sector names
        file.write('&')
        for c in countries:
            for s in sectors:
                file.write('&' + s)
        for c in countries:
            file.write('& Cons & Inv')
        file.write('\\\\\n')
        file.write('\\midrule\n')

        for i in range(0,nc):
            file.write('\\multirow{'+str(ns)+'}{*}{\\begin{sideways}'+countries[i]+'\\end{sideways}}')
            for ii in range(0,ns):
                file.write('&'+sectors[ii])
                for j in range(0,nc):
                    for jj in range(0,ns):
                        if M[i*ns+ii,j*ns+jj]>1.0e-6:
                            tmpstr = locale.format('%0.2f',M[i*ns+ii,j*ns+jj],grouping=True)
                        else:
                            tmpstr='-'
                        file.write('&'+tmpstr)
                for j in range(0,nc):
                    if Fc[i*ns+ii,j]>1.0e-6:
                        tmpstr = locale.format('%0.2f',Fc[i*ns+ii,j],grouping=True)
                    else:
                        tmpstr = '-'
                    file.write('&'+tmpstr)
                    if Fx[i*ns+ii,j]>1.0e-6:
                        tmpstr = locale.format('%0.2f',Fx[i*ns+ii,j],grouping=True)
                    else:
                        tmpstr='-'
                    file.write('&'+tmpstr)
                file.write('\\\\\n')
            
        file.write('\\midrule\n')
        file.write('VA &')
        for i in range(0,nc):
            for ii in range(0,ns):
                tmpstr = locale.format('%0.2f',V[i*ns+ii],grouping=True)
                file.write('&'+tmpstr)

        file.write('\\\\\n')

        file.write('\\midrule\n')
        file.write('GO &')
        for i in range(0,nc):
            for ii in range(0,ns):
                tmpstr = locale.format('%0.2f',colsums2[i*ns+ii],grouping=True)
                file.write('&'+tmpstr)

        file.write('\\\\\n')
            
        file.write('\\bottomrule\n')
        file.write('\\end{tabular}\n')
        file.write('\\end{center}\n')
        file.write('\\normalsize\n')
        file.write('\\end{table}\n')
        #file.write('\\end{landscape}\n')

write_iomat_csv(iomat2,'iomat_'+suff+'.txt')
write_iomat_latex(iomat2,
                  rowsums,
                  colsums,
                  '1995 input-output table (U.S. GDP = 100)',
                  'iomat_'+suff,
                  'iomat_'+suff)

##########################################################################################
# construct the alternative no-IO matrix
##########################################################################################

m.loc[:,'M'] = 0.0

vg = vg[vg.year==year].sort_values(by=['col_region_num','col_sector']).reset_index()
f = f[f.year==year].sort_values(by=['col_region_num','row_region_num','row_sector']).reset_index()
m = m[m.year==year].sort_values(by=['col_region_num','col_sector','row_region_num','row_sector']).reset_index()

vg['col_region_num'] = vg['col_region'].apply(assign_region_num)
f['col_region_num'] = f['col_region'].apply(assign_region_num)
f['row_region_num'] = f['row_region'].apply(assign_region_num)
m['col_region_num'] = m['col_region'].apply(assign_region_num)
m['row_region_num'] = m['row_region'].apply(assign_region_num)

rowsums = np.zeros( nc*ns + 1 )
colsums = np.zeros( nc*ns + nc*2 )

MM = m.pivot_table(values='M', index=['row_region_num','row_sector'], columns=['col_region_num','col_sector'])
VV = vg['VA'].values.reshape((1,nc*ns))
FF = f.pivot_table(values=['C','I'], index=['row_region_num','row_sector'], columns=['col_region_num'])
VV = np.hstack((VV,np.zeros((1,nc*2))))

iomat=np.vstack( ( np.hstack((MM,FF)) , VV ) )

for row in range(0,nc*ns + 1):
    rowsums[row] = np.sum(iomat[row,:])

for col in range(0,nc*ns + nc*2):
    colsums[col] = np.sum(iomat[:,col])

# make sure it's balanced after imposing the restrictions on construction usage...
colsums[0:(nc*ns)] = rowsums[0:(nc*ns)] # make sure markets clear: gross output = total demand for each country/sector
rowsums[-1] = colsums[(nc*ns):].sum() # world value added must equal world final demand
iomat2 = ras(iomat,rowsums,colsums) # run RAS

write_iomat_csv(iomat2,'iomat_noio_'+suff+'.txt')
write_iomat_latex(iomat2,
                  rowsums,
                  colsums,
                  '1995 input-output table (No-IO alternative, U.S. GDP = 100)',
                  'iomat_noio_'+suff,
                  'iomat_noio_'+suff)
