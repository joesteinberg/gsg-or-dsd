##################################################################################
# setup

# imports
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

mpl.rc('font',**{'family':'serif','serif':['Palatino'],'size':20})
mpl.rc('font',size=20)
mpl.rc('text', usetex=True)
mpl.rc('lines',linewidth=1.5)
mpl.rc('savefig',bbox='tight')
mpl.rc('savefig',format='pdf')


colors=['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00']

###############################################################################

OM = [0.9,0.8,0.7,0.6,0.5]
labels=[r'$\omega=%0.1f$'%om for om in OM]
TAU = np.arange(1.0,1.5,0.01)

def R(tus,trw,om):
    top = ( (om*tus*(1.0+trw) + (1.0-om)*trw*(1.0+tus))**(-om) + 
            ((1.0-om)*tus*(1.0+trw) + om*trw*(1.0+tus))**(om-1) )
    bottom = ( (om*(1.0+trw) + (1.0-om)*(1.0+tus))**(-om) + ((1.0-om)*(1.0+trw) + om*(1.0+tus))**(om-1) )
    return top/bottom

def tby(tus,trw,om):
    return (1-om)(((1+tus)/(1+trw))**2.0*(1+tus)-(1+trw))*(om*(1+trw)+(1-om)*(1+tus))

cnt=0
fix,ax=plt.subplots()
for om in OM:
    ax.plot(TAU,[R(1.0,t,om) for t in TAU],colors[cnt])
    cnt=cnt+1

plt.legend(labels,loc='upper left',prop={'size':10})
plt.savefig('output/2period_test1.pdf')
plt.close()

cnt=0
fix,ax=plt.subplots()
for om in OM:
    ax.plot(TAU,[tby(1.0,t,om) for t in TAU],colors[cnt])
    cnt=cnt+1

plt.legend(labels,loc='upper left',prop={'size':10})
plt.savefig('output/2period_test2.pdf')
plt.close()


