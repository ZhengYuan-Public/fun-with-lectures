import numpy as np
import pylab as pl
import h5py
import sys

n = len(sys.argv)

fs = 15
lw = 2
F = pl.figure()
f1 = F.add_subplot(111)

legend = []
tmax = 0

for i in range(1,n):
    fname = sys.argv[i]
    legend.append(fname)
    h = h5py.File(fname,'r')
    pins = h['pincount'][:]
    fqs = h['frequency'][:]
    t = h['times'][:]
    h.close()
    density = pins / (fqs**2)
    f1.plot(t,density,color=(0,0,0),linewidth=lw)
    if(t[-1]>tmax): tmax=t[-1]

f1.plot([0,tmax],[np.pi,np.pi],'--')
f1.axis([0,tmax,0,100])
f1.set_aspect(tmax/100)
f1.set_yticks(np.hstack([f1.get_yticks(),np.pi]))
f1.set_xlabel('time',fontsize=fs)
f1.set_ylabel('pinwheel density',fontsize=fs)
f1.legend(legend,frameon=False)

pl.show()
