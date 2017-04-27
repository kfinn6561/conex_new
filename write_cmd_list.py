import numpy as np

N0=0.03
calc_no_cl=True
Mss=[10,50,100,150,200]
Mcs=[0.5,0.17,0.1,0.05]
Fs=[0.9,0.5,0.1]

Mss=[]
Es=np.linspace(17.8,20,12)
nsamples=50
percmd= 10
outfile='cmd_list.txt'
f=open(outfile,'w')

if calc_no_cl:
    for E in Es:
        cmd='bin/conex2r -n %d -E %g -e %g -Q -x no_cl ' %(percmd,E,E)
        print cmd
        for i in range(nsamples):
            f.write(cmd+'\n')
    


for Ms in Mss:
    for Mc in Mcs:
        for F in Fs:
            for E in Es:
                cmd='bin/conex2r -E %g -e %g -n %d' %(E,E,percmd)
                cmd += ' -C %d -N %g -c %g -f %g -x cl_%d_%03d_%02d '%(Ms,N0,Mc,F,Ms,int(1000*Mc),int(100*F))
                print cmd
                for i in range(nsamples):
                    f.write(cmd+'\n')

f.close() 
print 'done'
