#!/usr/bin/env python3
# coding: utf-8

import ROOT
import array,math
import time
import numpy as np
import sys

from Util import *


workarea='workarea/'
PRINT_EVT_INTERVAL = 50000

# ## Reading LUT
maxEntries=-10000
scaleFactorLUT=readLUT('swetasRun3MCRegressionFinalCorrectedLUT_V0.txt')
energyLUT = readLUT('data/egCompressELUT_4bit_v4.txt')
shapeLUT  = readLUT('data/egCompressShapesLUT_calibr_4bit_v4.txt')
ietaLUT   = readLUT('data/egCompressEtaLUT_4bit_v4.txt')
outputFileName=workarea+"ntupleBasedRemulatedFile_V0_NTuple_crab_4841files_fromSweta_scale1p0.root"

if len(sys.argv) > 1:
    scaleFactorLUT=readLUT(sys.argv[1])
if len(sys.argv) > 2:
    outputFileName=workarea + sys.argv[2]

# Setting up the TandPNtuple
TandPNtuple = ROOT.TFile("../store/l1egamma/NTuple_crab_4841files_fromSweta.root","READ")
TandPNtuple.Ntuplizer.cd()
TandPNtupleTree = TandPNtuple.Get("Ntuplizer/TagAndProbe")


# Setting Up Branches
source_branches=['eleProbeTriggerBits','eleTagTriggerBits','eleProbePt','eleProbeEta',
                 'eleProbePhi','eleProbeSclEt','eleProbeCharge','eleTagPt','eleTagEta',
                 'eleTagPhi','eleTagCharge','Mee','hltPt','hltEta','hltPhi','l1tPt',
                 'l1tEta','l1tPhi','l1tQual','l1tIso','l1tEmuPt','l1tEmuEta','l1tEmuPhi',
                 'l1tEmuQual','l1tEmuIso','l1tEmuNTT','l1tEmuTowerIEta','l1tEmuTowerIPhi',
                 'l1tEmuRawEt','l1tEmuIsoEt','isProbeLoose','isProbeTight','isProbeMedium',
                 'UnpackedDR','EmulatedDR','seedTowerEt','TowerHoE','iem','ihad','shape',
                 'isTagHLTmatched','isProbeHLTmatched','isOS','foundTag','Nvtx']

target_branches = source_branches
target_branches.remove("l1tEmuPt")

TandPNtupleTree.SetBranchStatus("*",0)
for br in source_branches:
    TandPNtupleTree.SetBranchStatus(br,1)


# In[13]:

def addTriggerBranch(outTree,templateName,startPt,endPt,step=1):
    triggerDict={}
    pt=startPt
    while pt<= endPt:
        triggerDict.update({pt:array.array('I',[0])})
        outTree.Branch(templateName.format(pt), triggerDict[pt],templateName.format(pt)+"/I");
        pt+=step
    return triggerDict


dataTypeDict={'Bool_t':'B','Float_t':'f','Int_t':'i','ULong64_t':'Q'}
dataTypeDict_root={'Bool_t':'O','Float_t':'F','Int_t':'I','ULong64_t':'l'}
for br in target_branches:
    dType=TandPNtupleTree.GetLeaf(br).GetTypeName()
    if dType not in dataTypeDict:
        print(dType," not foud in abbrev. dict")


outData={}
outDataType={}
for br in target_branches:
    dType=TandPNtupleTree.GetLeaf(br).GetTypeName()
    dTypeAbbv=dataTypeDict[dType]
    outData.update({br:array.array(dTypeAbbv,[0])})
    outDataType.update({br:dataTypeDict_root[dType]})

# ### Copying branches to new files

for br in target_branches:
    if br in source_branches:
        TandPNtupleTree.SetBranchAddress(br,outData[br])


# ### Output Tree
target_branches=source_branches
outputFile = ROOT.TFile.Open(outputFileName, "RECREATE")
outputFile.cd()
outputTree = ROOT.TTree("TagAndProbe","Ntuple based reemulation of L1")

for name, a in outData.items():
    outputTree.Branch(name, a, "{0}/{1}".format(name,outDataType[name] ))

newData={}
reEmulatedEt=array.array('f',[0])
outputTree.Branch('l1tEmuPt', reEmulatedEt,"l1tEmuPt/F");

triggerDict= addTriggerBranch(outputTree,"hasL1ReEmu_{}",1,60,1)

# In[ ]:


nentries = TandPNtupleTree.GetEntries()
if maxEntries>0 : 
    nentries= min(nentries,maxEntries)
print("Number of entries to process : ",nentries)

startT = time.time()
for e in range(nentries):
    if(e%PRINT_EVT_INTERVAL==0) :
        deltaT= time.time() -startT
        tLeft=(nentries-e)*deltaT/(e +1e-9)
        print("Processing event ",e," /",nentries,"  [ ",100.0*e/nentries  ," ] ","  Elapsed time : ","{0:0.2f}".format(deltaT) ,
              ' secs , Est. time left : ',"{0:0.2f}".format(tLeft))
    TandPNtupleTree.GetEntry(e)
    
    if TandPNtupleTree.l1tEmuRawEt <0 :
        continue

    cd_shape =  shapeLUT[ TandPNtupleTree.shape ] 
    cd_ieta  =  int(math.copysign(ietaLUT[abs(TandPNtupleTree.l1tEmuTowerIEta)], 
                                  TandPNtupleTree.l1tEmuTowerIEta))
    cd_E     =  energyLUT[ min(TandPNtupleTree.l1tEmuRawEt ,255 )]
    LUT_index     = getIndex(cd_shape,cd_E,abs(cd_ieta))
    calibLUTvalue = scaleFactorLUT[LUT_index]
    scaleFactor   = getScaleFactor(calibLUTvalue)
    reEmulatedEt[0] = min(scaleFactor*TandPNtupleTree.l1tEmuRawEt,255)
    reEmulatedEt[0] =int(reEmulatedEt[0])
    #if ( (int(reEmulatedEt[0]*10))%10 < 5 ) : reEmulatedEt[0] = int(reEmulatedEt[0]) + 0.5
    #else : reEmulatedEt[0] =int(reEmulatedEt[0])

    for pt in triggerDict:
        if reEmulatedEt[0] > pt :
            triggerDict[pt][0]=1
        else:
            triggerDict[pt][0]=0
    outputTree.Fill()

outputTree.Write()
outputFile.Purge()
outputFile.Close()
print("Saved file to  : ",outputFileName)
