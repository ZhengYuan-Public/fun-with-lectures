import h5py
import numpy as np
import pylab as pl
import scipy as sc
from os import listdir
import sys

from PIL import Image

def getPic(fname):
    img = Image.open(fname)
    img.thumbnail((256, 256), Image.ANTIALIAS)
    p = np.array(img.getdata())/255
    try:
        p = p.reshape(img.size[0], img.size[1], 3)
        p = np.mean(p,axis=2)
    except:
        pass

    #img.show()
    p = np.flipud(p)
    return p.T.flatten()

P = np.array([],dtype=float)

if(len(sys.argv)<4):
    print('ERROR: supply directory, filetype, outputfile, e.g., python genTestPatterns.py shouval png ../natural.h5')
    exit()

dir = sys.argv[1] #'images/shouval'
typ = sys.argv[2] #'png'
outfile = sys.argv[3] # 'naturalImages.h5'

n=0
for f in listdir(dir):
    s = f.split('.')
    if(s[1]==typ):
        P = np.hstack([P, getPic(dir+'/'+f)])
        print('added: '+dir+'/'+f)
        n+=1


print(P.shape)
print(n)

Pscale = np.max(P)-np.min(P)
P = (P-np.min(P))/Pscale

h5f = h5py.File(outfile,'w')
h5f.create_dataset('P', data=P)
h5f.create_dataset('dims', data=np.array([256,256,n]))
h5f.close()
