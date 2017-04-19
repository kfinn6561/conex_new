import multiprocessing
import os
import subprocess
import time
import sys

Ncores=12
max_runtime=5.*3600#5 hours should be long enough for any shower. If it takes longer than this it has probably crashed

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
            out+='%d: %s. ' %(i+1,hms(now-start_times[i]))
        else:
            out+='%d: stopped ' %(i+1)
    overprint(out)
    
def check_for_deaths(threads,start_times):
    now=time.time()
    out=[]
    for i in range(len(start_times)):
        if start_times[i] and (now-start_times[i]>max_runtime):
            print '\ncore %d has become stuck. Terminating' %(i+1)
            threads[i].terminate()
            out.append(i)
    return out
        
        
completed_fname='completed.dat'
completed_jobs=[]
try:
    if 'restart' in sys.argv:
        os.remove(completed_fname)
    else:
        f=open(completed_fname,'r')
        for line in f:
            completed_jobs.append(int(line))
except:
    pass



print 'using %d cores' %Ncores
f=open('cmd_list.txt','r')
cmds=f.readlines()
f.close()
#cmds=cmds[:9]#debug

threads=[fakethread() for i in range(Ncores)] #there may be a better way to do this
start_times=[False for i in range(Ncores)]
running_jobs=[False for i in range(Ncores)]

job_no=0
while job_no<len(cmds):
    for i in range(Ncores):
        while job_no in completed_jobs:
            job_no+=1
        if job_no<len(cmds) and not threads[i].is_alive():
            cmd=cmds[job_no][:-2]
            print '\njob %d of %d: running %s on core %d' %(job_no+1,len(cmds),cmd,i+1)
            #threads[i]=threading.Thread(target=run_conex,args=(cmd,'log_%d.txt' %(i+1),))
            threads[i]=multiprocessing.Process(target=run_conex,args=(cmd,'log_%d.txt' %(i+1),))
            threads[i].start()
            start_times[i]=time.time()
            if running_jobs[i]:
                f=open(completed_fname,'a')
                f.write('%d\n'%running_jobs[i])
                f.close()
            running_jobs[i]=job_no
            job_no+=1
    deaths=check_for_deaths(threads, start_times)
    for i in deaths:
        start_times[i]=False
        running_jobs[i]=False
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
            if running_jobs[i]:
                f=open(completed_fname,'a')
                f.write('%d\n'%running_jobs[i])
                f.close()
    deaths=check_for_deaths(threads, start_times)
    for i in deaths:
        start_times[i]=False
        running_jobs[i]=False
    timing_update(start_times)
    time.sleep(5)


for thread in threads:#this should be unecessary
    thread.join()
os.remove(completed_fname)
print '\nAll jobs completed'
