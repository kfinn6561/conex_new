import threading
import os
import subprocess
import time

def run_conex(cmd,logname):
    os.environ['ROOTSYS']=os.environ['HOME']+'/root-install'
    os.environ['PATH']+=os.pathsep+os.environ['ROOTSYS']+'/bin'
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
        
Ncores=4

print 'using %d cores' %Ncores
start=time.time()
f=open('cmd_list.txt','r')
cmds=f.readlines()
f.close()
cmds=cmds[:9]#debug

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

print 'All jobs completed'
