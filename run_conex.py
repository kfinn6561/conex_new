import os

def name_last_conex(fname):
    dirlist=os.listdir('.')
    rname=[x for x in dirlist if 'conex_eposlhc' in x][0]
    os.rename(rname,fname)
    return


os.system("bin/conex2r -e 17 -E 17 -C 1e20")
name_last_conex("no_classicalization.root")


for n in [0.001,0.01,0.1,1]:

    print '\n\n\nN=%f\n\n\n'%n
    os.system("bin/conex2r -e 17 -E 17 -C 10 -N %f" %n)
    name_last_conex('conex_%f.root' %n)

print 'done'
