from __future__ import print_function 
import os,sys,itertools

printConfig=True
if len(sys.argv) > 1:
    if int(sys.argv[1]) < 1:
        printConfig=False

etMins=[5.0,10.0,15.0,20.0,25.0,30.0,35.0,40.0,45.0,50.0]
effMins=[0.1,0.5,0.9]
etMaxs=[5.0,10.0,15.0,20.0,25.0,30.0,35.0,40.0,45.0,50.0]

allOptsProd=itertools.product(*[etMins,effMins,etMaxs])
allOpts=[t for t in allOptsProd]
nOptsMax=len(allOpts)
print('Options in the cartesian product : ',nOptsMax)

optId=0
for i in range(nOptsMax):
    if allOpts[0] > allOpts[2]:
        continue
    optId+=1
    optstr=str(optId)+'_'+str(allOpts[i][0]).replace('.','p')+'_'+str(allOpts[i][1]).replace('.','p')+'_'+str(allOpts[i][2]).replace('.','p')
    if printConfig:
        optstr+=':'+str(allOpts[i][0])+':'+str(allOpts[i][1])+':'+str(allOpts[i][2])    
    optstr+=','
    print(optstr,end='')

print()
print()
print('Number of options made : ',optId)
print()
