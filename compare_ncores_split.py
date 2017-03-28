import threading
import os
import subprocess
import time

def run_conex(cmd_list,logname):
    os.environ['ROOTSYS']=os.environ['HOME']+'/root-install'
    os.environ['PATH']+=os.pathsep+os.environ['ROOTSYS']+'/bin'
    os.environ['ROOT_OUT']='/export/ursa1/ktf243'#specific to ursa
    try:
        os.environ['LD_LIBRARY_PATH']+=os.pathsep+os.environ['ROOTSYS']+'/lib/root'
    except KeyError:
        os.environ['LD_LIBRARY_PATH']=os.environ['ROOTSYS']+'/lib/root'        
    logfile=open(logname,'a')
    for cmd in cmd_list:
        subprocess.call([cmd],stdout=logfile,shell=True)
    logfile.close()

    
class fakethread():
    def __init__(self):
        return
    def is_alive(self):
        return False
    def join(self):
        return True
        
Ns=range(4,15)
outfile=open('compare_ncores_split.txt','w')

for Ncores in Ns:
    print 'using %d cores' %Ncores
    start=time.time()
    f=open('cmd_list.txt','r')
    cmds=f.readlines()
    f.close()
    #cmds=cmds[:9]#debug
    
    threads=[fakethread() for i in range(Ncores)] #there may be a better way to do this
    no_jobs=len(cmds)
    jobs_per_core=int(no_jobs)/int(Ncores)
    jobs=[jobs_per_core for i in Ncores]
    for i in range(no_jobs-jobs_per_core*Ncores):
        jobs[i]+=1    
    job_no=0
    for i in range(Ncores):
        cmd_list=cmds[job_no:job_no+jobs[i]]
        job_no+=jobs[i]
        print 'issuing %d jobs to core %d' %(jobs[i],i+1)
        threads[i]=threading.Thread(target=run_conex,args=(cmd_list,'log_%d.txt' %(i+1),))
        threads[i].start()
    
    print 'All jobs assigned. Waiting for completion'
    for thread in threads:
        thread.join()
    finish=time.time()
    outfile.write('%d\t%f\n'%(Ncores,finish-start))
outfile.close()
print 'All jobs completed'
