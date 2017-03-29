import numpy as np

N0=0.03
Ms=10
E=19.5
nsamples=10

Nsamples=60
groups=[1,2,3,4,5,6,10,20,30,60]

for group in groups:
    f=open('splittings_cmd_list_%d.txt' %group,'w')
    nsamples=Nsamples/group
    cmd='bin/conex2r -E %g -e %g -n %d' %(E,E,nsamples)
    if Ms==0:
        cmd +=' -Q -x no_cl '
    else:
        cmd += ' -C %d -N %g -x cl_%d '%(Ms,N0,Ms)
    for i in range(group):
        f.write(cmd+'\n')
    f.close() 
print 'done'
