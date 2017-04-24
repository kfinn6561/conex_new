import numpy as np

N0=0.03
calc_no_cl=True
Mss=[10,50,100,150,200]
Mcs=[0.5,0.17,0.1,0.05]
Fs=[0.9,0.5,0.1]

Es=np.linspace(17.8,20,12)
nsamples=50
outfile='cmd_list.txt'
f=open(outfile,'w')

if calc_no_cl:
    


for Ms in Mss:
    for E in Es:
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
