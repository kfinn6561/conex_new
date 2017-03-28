import threading
import os
import subprocess
import time

def run_conex(cmd,logname):
    os.environ['ROOTSYS']=os.environ['HOME']+'/root-install'
    os.environ['PATH']+=os.pathsep+os.environ['ROOTSYS']+'/bin'
    os.environ['ROOT_OUT']='/export/ursa1/ktf243'#specific to ursa
    try:
        os.environ['LD_LIBRARY_PATH']+=os.pathsep+os.environ['ROOTSYS']+'/lib/root'
    except KeyError:
        os.environ['LD_LIBRARY_PATH']=os.environ['ROOTSYS']+'/lib/root'        
    logfile=open(logname,'a')
    subprocess.call([cmd],stdout=logfile,shell=True)
    logfile.close()

    
class fakethread():
    def __init__(self):
        return
    def is_alive(self):
        return False
    def join(self):
        return True
        
groups=[1,2,3,4,5,6,10,20,30,60]
outfile='compare_splittings_pool.txt'
Ncores=12
print 'using %d cores' %Ncores

for group in groups:
    print 'splitting into %d groups' %group
    start=time.time()
    f=open('splittings_cmd_list_%d.txt' %group,'r')
    cmds=f.readlines()
    f.close()
    #cmds=cmds[:9]#debug
    
    threads=[fakethread() for i in range(Ncores)] #there may be a better way to do this
    
    job_no=0
    while job_no<len(cmds):
        for i in range(Ncores):
            if job_no<len(cmds) and not threads[i].is_alive():
                cmd=cmds[job_no][:-2]
                print 'job %d of %d: running %s on core %d' %(job_no+1,len(cmds),cmd,i+1)
                threads[i]=threading.Thread(target=run_conex,args=(cmd,'log_%d.txt' %(i+1),))
                threads[i].start()
                job_no+=1
    
    print 'All jobs assigned. Waiting for completion'
    for thread in threads:
        thread.join()
    finish=time.time()
    f=open(outfile,'a')
    f.write('%d\t%f\n'%(group,finish-start))
    f.close()
print 'All jobs completed'
