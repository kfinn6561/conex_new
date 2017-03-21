import pylab as p
from root_numpy import root2array
import os

dirlist=os.listdir('.')
fnames=[x for x in dirlist if '.root' in x]

nfig=p.figure()
nax=nfig.add_subplot(111)
p.xlabel('X')
p.ylabel('N')

efig=p.figure()
eax=efig.add_subplot(111)
p.xlabel('X')
p.ylabel('dEdX')


for fname in fnames:
    lname=fname[:-5]
    data=root2array(fname,'Shower')
    print fname
   # print 'log E: %f' %data['lgE'][0]

    X=data['X'][0]
    N=data['N'][0]
    E=data['dEdX'][0]

    nax.plot(X,N,label=lname)
    eax.plot(X,E,label=lname)

nax.legend()
eax.legend()
p.show()
