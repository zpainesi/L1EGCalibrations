#!/usr/bin/env python
from __future__ import print_function 
import os,sys,itertools


MAXEVENTS_RATE=-1
MAXEVENTS_EFF =2500000
printConfig=True
runTmplName='misc/runStep3and4.tpl.sh'
cfgTmplName='misc/Par_ApplyIsolation.dat.tpl.cfg'
optPerParSet=50
optionFile='/grid_mnt/t3storage3/athachay/l1egamma/isolation/CMSSW_12_0_2/src/L1EGCalibrations/Isolation/test/caliberationFiles/HistgramFile_step1step2_eg_rawet_12X_withQuality_optfile.root'

pwd=os.environ['PWD']
proxy_path=os.environ['X509_USER_PROXY']
HOME=os.environ['HOME']

if len(sys.argv) > 1:
    if int(sys.argv[1]) < 1:
        printConfig=False
tag=""
if len(sys.argv) > 2:
    tag=sys.argv[2]


def loadConfigTemplate(fname):
    f=open(fname,'r')
    txt=f.readlines()
    f.close()
    cfgString="".join(txt)
    return cfgString    


condorScriptString="\
executable = $(filename)\n\
output = $Fp(filename)run.$(Cluster).stdout\n\
error = $Fp(filename)run.$(Cluster).stderr\n\
log = $Fp(filename)run.$(Cluster).log\n\
+JobFlavour = \"espresso\"\n\
"

runScriptTxt=loadConfigTemplate(runTmplName)
cfgTmplate=loadConfigTemplate(cfgTmplName)

etMins=[5.0,7.5,10.0,12.5,15.0,17.5,20.0,22.5,25.0,27.5,30.0,32.5,35.0,37.5,40.0,42.5,45.0,50.0]
effMins=[0.1,0.2,0.4,0.5,0.7,0.8,0.9]
etMaxs=[5.0,7.5,10.0,12.5,15.0,17.5,20.0,22.5,25.0,27.5,30.0,32.5,35.0,37.5,40.0,42.5,45.0,50.0]

allOptsProd=itertools.product(*[etMins,effMins,etMaxs])
allOpts=[t for t in allOptsProd]
nOptsMax=len(allOpts)
print('Options in the cartesian product : ',nOptsMax)

optId=0
parIdx=0
njobs=0
if printConfig:
    optstr=''
    for i in range(nOptsMax):
        if allOpts[i][0] > allOpts[i][2]:
            continue
        optId+=1
        optstr+=str(optId)+'_'+str(allOpts[i][0]).replace('.','p')+'_'+str(allOpts[i][1]).replace('.','p')+'_'+str(allOpts[i][2]).replace('.','p')
        optstr+=':'+str(allOpts[i][0])+':'+str(allOpts[i][1])+':'+str(allOpts[i][2])    
        optstr+=','
    print(optstr[:-1],end='')
else :
    head='Condor/Jobs'+tag
    if not os.path.exists(head ):
        os.system('mkdir -p '+head)
    condorScriptName=head+'/job'+tag+'.sub'
    condorScript=open(condorScriptName,'w')
    condorScript.write(condorScriptString)

    optInParset=0
    optstr=''
    for i in range(nOptsMax):
        if allOpts[i][0] > allOpts[i][2]:
            continue
        optId+=1
        optInParset+=1
        optstr+=str(optId)+'_'+str(allOpts[i][0]).replace('.','p')+'_'+str(allOpts[i][1]).replace('.','p')+'_'+str(allOpts[i][2]).replace('.','p')
        optstr+=','
        if optInParset==optPerParSet:
            parIdx+=1
            dirName= pwd+'/'+head+'/Job_'+str(parIdx)
            if not os.path.exists(dirName):
                os.system('mkdir '+dirName)
            tmp=cfgTmplate.replace('@@OPTIONS',optstr[:-1])
            tmp=tmp.replace("@@OPT_FNAME",optionFile)
            tmp=tmp.replace("@@PWD",pwd)
            tmp=tmp.replace("@@TAG",tag+'p'+str(parIdx))
            tmp=tmp.replace("@@IDX",str(parIdx))
            tmp=tmp.replace("@@MAXEVENTS_EFF",str(MAXEVENTS_EFF))
            tmp=tmp.replace("@@MAXEVENTS_RATE",str(MAXEVENTS_RATE))
            cfgFileName=dirName+'/Par_ApplyIsolation'+tag+'p'+str(parIdx)+'.dat'
            cfgFile=open(cfgFileName,'w')
            cfgFile.write(tmp)
            cfgFile.close()
            
            runScriptName=dirName+'/'+tag+'run'+str(parIdx)+'.sh'
            if os.path.exists(runScriptName+'.sucess'):
                os.system('rm '+runScriptName+'.sucess')
            runScript=open(runScriptName,'w')
            tmp=runScriptTxt.replace("@@DIRNAME",dirName)
            tmp=tmp.replace("@@PROXY",proxy_path)
            tmp=tmp.replace("@@HOME",HOME)
            tmp=tmp.replace("@@PWD",pwd)
            tmp=tmp.replace("@@IDX",str(parIdx))
            tmp=tmp.replace("@@CFGFILENAME",cfgFileName)
            tmp=tmp.replace("@@RUNSCRIPT",runScriptName)
            runScript.write(tmp)
            runScript.close()
            os.system('chmod +x '+runScriptName)
            condorScript.write("queue filename matching ("+runScriptName+")\n")
            njobs+=1

            optInParset=0
            optstr=''
    if(optInParset > 0 ):
        parIdx+=1
        dirName= pwd+'/'+head+'/Job_'+str(parIdx)
        if not os.path.exists(dirName):
            os.system('mkdir '+dirName)
        tmp=cfgTmplate.replace('@@OPTIONS',optstr[:-1])
        tmp=tmp.replace("@@TAG",tag+'p'+str(parIdx))
        tmp=tmp.replace("@@PWD",pwd)
        tmp=tmp.replace("@@OPT_FNAME",optionFile)
        tmp=tmp.replace("@@IDX",str(parIdx))
        tmp=tmp.replace("@@MAXEVENTS_EFF",str(MAXEVENTS_EFF))
        tmp=tmp.replace("@@MAXEVENTS_RATE",str(MAXEVENTS_RATE))
        cfgFileName=dirName+'/Par_ApplyIsolation'+tag+'p'+str(parIdx)+'.dat'
        cfgFile=open(cfgFileName,'w')
        cfgFile.write(tmp)
        cfgFile.close()
        
        runScriptName=dirName+'/'+tag+'run'+str(parIdx)+'.sh'
        if os.path.exists(runScriptName+'.sucess'):
            os.system('rm '+runScriptName+'.sucess')
        runScript=open(runScriptName,'w')
        tmp=runScriptTxt.replace("@@DIRNAME",dirName)
        tmp=tmp.replace("@@PWD",pwd)
        tmp=tmp.replace("@@IDX",str(parIdx))
        tmp=tmp.replace("@@CFGFILENAME",cfgFileName)
        tmp=tmp.replace("@@RUNSCRIPT",runScriptName)
        runScript.write(tmp)
        runScript.close()
        os.system('chmod +x '+runScriptName)
        condorScript.write("queue filename matching ("+runScriptName+")\n")
        njobs+=1

        optInParset=0
        optstr=''
    
    print()
    print(' Number of options made : ',optId)
    print(" Number of jobs made : ", njobs)
    print(" Condor  submit file name : ", condorScriptName)
    condorScript.close()
