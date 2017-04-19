import numpy as np

N0=0.03
Mss=[15,20,50,100]
Es=np.linspace(17.8,20,12)
#Edict={18:10,18.3:10,18.7:10,19:10,19.3:5,19.7:1,20:1}
#Es=Edict.keys()
nsamples=20
persim=10
outfile='cmd_list.txt'
f=open(outfile,'w')


for Ms in Mss:
    for E in Es:
        cmd='bin/conex2r -E %g -e %g -n %d' %(E,E,persim)
        if Ms==0:
            cmd +=' -Q -x no_cl '
        else:
            cmd += ' -C %d -N %g -c 0.17 -x cl_%d '%(Ms,N0,Ms)
        print cmd
        for i in range(nsamples):
            f.write(cmd+'\n')

f.close() 
print 'done'
