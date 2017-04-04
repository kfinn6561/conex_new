import numpy as np

N0=0.03
Mss=[1,5,7,10,12,15]
Es=np.linspace(18,20,11)
#Edict={18:10,18.3:10,18.7:10,19:10,19.3:5,19.7:1,20:1}
#Es=Edict.keys()
nsamples=50
outfile='cmd_list.txt'
f=open(outfile,'w')

Mss=[10,12]
Es=[18.7,19.5]
nsamples=1000
'''
for Ms in Mss:
    for E in Es:
'''
for i in range(2):
    Ms=Mss[i]
    E=Es[i]
    cmd='bin/conex2r -E %g -e %g' %(E,E)
    if Ms==0:
        cmd +=' -Q -x no_cl '
    else:
        cmd += ' -C %d -N %g -c 0.17 -x cl_%d '%(Ms,N0,Ms)
    print cmd
    for i in range(nsamples):
        f.write(cmd+'\n')

f.close() 
print 'done'
