import numpy as np

N0=0.03
Mss=[0,1,5,10,20,50,100]
Es=np.linspace(18,20,11)
Edict={18:10,18.3:10,18.7:10,19:10,19.3:5,19.7:1,20:1}
Es=Edict.keys()
nsamples=10
outfile='cmd_list.txt'
f=open(outfile,'w')



for Ms in Mss:
    for E in Es:
        cmd='bin/conex2r -E %g -e %g -n %d' %(E,E,nsamples)
        if Ms==0:
            cmd +=' -Q -x no_cl '
        else:
            cmd += ' -C %d -N %g -x cl_%d '%(Ms,N0,Ms)
        print cmd
        for i in range(Edict[E]):
            f.write(cmd+'\n')

f.close() 
print 'done'