#!/usr/bin/python

import sys
sys.dont_write_bytecode=True
import os
import shutil
import write_input_file as wif


print "Hello, World!. Lets start some awesome TIS."

def Read_Parameters():
    fparam = 'params.inp'
    global Params
    with open(fparam) as f:
        Params=[[str(x) for x in line.split()] for line in f]
    for row in Params:
        row.sort()

Read_Parameters()
print Params

def ensure_dir(f):
    if not os.path.exists(f):
        os.mkdir(f)
    return;

def warmstartgen():
   fjobid=open("jobid","r")
   a=fjobid.readline()
   a=a.split()
   jobid=a[len(a)-1]
   warmdir = 'bdtteps12-' + jobid
   os.chdir(warmdir)
   warmsucces=0
   #if os.path.exists("data.tar"):
   #    os.system("tar xvf data.tar dos_all.dat trajectory.out")
   if os.path.exists("dos_all.dat"):
       os.system("cp dos_all.dat ../dos_all.inp")
       warmsucces=1
   if os.path.exists("trajectory.out"):
       os.system("cp trajectory.out ../trajectory.inp")
       warmsucces=1
   os.chdir('..')
   if(warmsucces==1):
       shutil.copy("jobid","prevjobid")
   return warmsucces


for Eps in Params[0]:
    fdir = 'f'+str(Eps)
    ensure_dir(fdir)
    os.chdir(fdir)
    os.system("cp ../write_input_file.py ../submit.job ../lambda_u.inp ../lambda_b.inp ../bmd.run .")
    Ncycle1 = 1000
    Ncycle2 = 1000
    Type=0
    warmsucces=0
    if(Type==1):
        warmsucces=warmstartgen()
    if((warmsucces) or (Type==0)):
        wif.write_input(Eps,Type,Ncycle1,Ncycle2)
        #os.system("qsub submit.job")
    os.system("pwd")
    os.chdir('..')


