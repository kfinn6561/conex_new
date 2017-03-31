import threading
import os
import subprocess
import time
import sys

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
 
    
def overprint(s):
    sys.stdout.write('\r')
    sys.stdout.flush()
    sys.stdout.write(s)
    
def hms(time): #converts a number in seconds into a string of the form HH:MM:SS
    hours=int(time/3600)
    time-=hours*3600
    minutes=int(time/60)
    seconds=int(time-minutes*60)
    out= '%d:%02d:%02d' %(hours,minutes,seconds)
    return out

class fakethread():
    def __init__(self):
        return
    def is_alive(self):
        return False
    def join(self):
        return True
    
def timing_update(start_times):
    now=time.time()
    out=''
    for i in range(len(start_times)):
        if start_times[i]:
            out+='core %d: %s. ' %(i+1,hms(now-start_times[i]))
        else:
            out+='core %d: stopped ' %(i+1)
    overprint(out)
        
Ncores=12

print 'using %d cores' %Ncores
f=open('cmd_list.txt','r')
cmds=f.readlines()
f.close()
#cmds=cmds[:9]#debug

threads=[fakethread() for i in range(Ncores)] #there may be a better way to do this
start_times=[False for i in range(Ncores)]


job_no=0
while job_no<len(cmds):
    for i in range(Ncores):
        if job_no<len(cmds) and not threads[i].is_alive():
            cmd=cmds[job_no][:-2]
            print '\njob %d of %d: running %s on core %d' %(job_no+1,len(cmds),cmd,i+1)
            threads[i]=threading.Thread(target=run_conex,args=(cmd,'log_%d.txt' %(i+1),))
            threads[i].start()
            start_times[i]=time.time()
            job_no+=1
    timing_update(start_times)
    time.sleep(5)

print '\nAll jobs assigned. Waiting for completion'

running_threads=range(Ncores)

while len(running_threads)>0:
    for i in range(Ncores):
        if i in running_threads and not threads[i].is_alive():
            print '\nCore %d has finished running: %d jobs remaining' %(i+1,len(running_threads)-1)
            running_threads.remove(i)
            start_times[i]=False
    timing_update(start_times)
    time.sleep(5)


for thread in threads:#this should be unecessary
    thread.join()

print '\nAll jobs completed'
