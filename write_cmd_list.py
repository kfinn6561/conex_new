import numpy as np

N0=0.03
Mss=[0,1,2,5,10,20,50,100]
Es=np.linspace(18,21,11)
nsamples=10
data_dir=''
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
        f.write(cmd+'\n')

f.close() 
print 'done'
