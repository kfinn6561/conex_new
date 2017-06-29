import numpy as np

calc_no_cl=False

forward_threshold=0.1

Es=np.linspace(17.8,20,12)
nsamples=200
percmd= 1
outfile='cmd_list.txt'
f=open(outfile,'w')

if calc_no_cl:
    for E in Es:
        cmd='bin/conex2r -n %d -E %g -e %g -Q -x no_cl ' %(percmd,E,E)
        print cmd
        for i in range(nsamples):
            f.write(cmd+'\n')
    

for E in Es:
    cmd='bin/conex2r -E %g -e %g -n %d' %(E,E,percmd)
    cmd += ' -w %.3f -x forward_%d '%(forward_threshold,int(forward_threshold*100))
    print cmd
    for i in range(nsamples):
        f.write(cmd+'\n')

f.close() 
print 'done'
