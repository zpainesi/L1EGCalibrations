#!/usr/bin/env python
import ROOT
import array
import math
import operator
import time
import sys



compressedIetaFile  = "compressionLuts/tauCompressEtaLUT_2bit_v8.txt"
compressedEFile     = "compressionLuts/tauCompressELUT_5bit_v8.txt"
compressedNTTFile   = "compressionLuts/tauCompressnTTLUT_5bit_v8.txt"

inputFileName       = "/grid_mnt/t3storage3/athachay/l1egamma/data/run3MC/L1Ntuple.root"
outputFileName      = "EGIsolationStudyNtuple_tauBasedCopressionLutsForL1Ntuple.root"
treeName = "L1UpgradeTree"
maxEntries= 10+ 0*-1*int(1e6)

if len(sys.argv)>1:
    inputFileName=sys.argv[1]
if len(sys.argv)>2:
    maxEntries=int(sys.argv[2])
if len(sys.argv)>3:
    outputFileName=sys.argv[3]

def readLUT(lutFileName,superCompressionFactor=1):
    lut = {}
    with open(lutFileName, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if line[0]=='#': continue
            line=line.split("#")[0]
            tokens = line.split()
            if len(tokens)!=2: continue
            lut[int(tokens[0])] = int(int(tokens[1])/superCompressionFactor)
    return lut

## Read Ieta,E,shape compression mapping
compressedIeta  = readLUT(compressedIetaFile)
compressedE     = readLUT(compressedEFile)
compressedNTT   = readLUT(compressedNTTFile)

## Reading input trees 
## + filling compressed shape histo
inputFile = ROOT.TFile.Open(inputFileName)
inputTree = inputFile.Get(treeName)
inputTree.__class__ = ROOT.TTree


#inputTree.SetBranchStatus("*",0);
#
#inputTree.SetBranchStatus("RunNumber",1);
#inputTree.SetBranchStatus("EventNumber",1);
#inputTree.SetBranchStatus("eleProbeSclEt",1);
#inputTree.SetBranchStatus("eleProbeEta",1);
#inputTree.SetBranchStatus("eleProbePhi",1);
#inputTree.SetBranchStatus("l1tEmuTowerIEta",1);
#inputTree.SetBranchStatus("l1tEmuRawEt",1);
#inputTree.SetBranchStatus("l1tEmuPt",1);
#inputTree.SetBranchStatus("l1tEmuNTT",1);
#inputTree.SetBranchStatus("l1tEmuIsoEt",1);
#inputTree.SetBranchStatus("l1tEmuIso",1);


data = {
        "Run"    : array.array('i',[0]),
        "Event"  : array.array('l',[0]),
        "offlineEt"    : array.array('f',[0]),
        "offlineEta"   : array.array('f',[0]),
        "offlinePhi"   : array.array('f',[0]),
        "l1RawE"       : array.array('i',[0]),
        "l1Et"         : array.array('f',[0]),
        "ieta"         : array.array('i',[0]),
        "nTT"          : array.array('i',[0]),
        "isoEt"        : array.array('i',[0]),
        "isoCat"       : array.array('i',[0]),
        "compressedieta" : array.array('i',[0]),
        "compressedE"    : array.array('i',[0]),
        "compressedNTT"  : array.array('i',[0]),
       }

outputFile = ROOT.TFile.Open(outputFileName, "RECREATE")
outputFile.cd()
outputTree = ROOT.TTree("EGIsoCalibration","EGIsoCalibration")

for name, a in data.items():
    outputTree.Branch(name, a, "{0}/{1}".format(name, a.typecode.upper()))

nentries = inputTree.GetEntries()
if maxEntries>0 : 
    nentries= min(nentries,maxEntries)
print "Number of entries to process : ",nentries

startT = time.time()
for e in xrange(nentries):
    if(e%20000==0) :
        deltaT= time.time() -startT
        tLeft=(nentries-e)*deltaT/(e +1e-9)
        print "Processing event ",e," Elapsed time : ",deltaT ,' secs , Est. time left : ',tLeft
    inputTree.GetEntry(e)
    print e
    print inputTree.nEGs
#    data["Run"][0]         =  int(inputTree.RunNumber       )
#    data["Event"][0]       =  int(inputTree.EventNumber     )
#    data["offlineEt"][0]   =  int(inputTree.eleProbeSclEt   )
#    data["offlineEta"][0]  =  int(inputTree.eleProbeEta     )
#    data["offlinePhi"][0]  =  int(inputTree.eleProbePhi     )
#    data["l1RawE"][0]      =  int(inputTree.l1tEmuRawEt     )
#    data["l1Et"][0]        =  int(inputTree.l1tEmuPt        )
# #   print(data["l1Et"][0])
#    if inputTree.l1tEmuRawEt < 0 :
#        continue
#    data["ieta"][0]        =  int(inputTree.l1tEmuTowerIEta )
#    data["nTT"][0]         =  int(inputTree.l1tEmuNTT       )
#    data["isoEt"][0]       =  int(inputTree.l1tEmuIsoEt     )
#    data["isoCat"][0]      =  int(inputTree.l1tEmuIso       )
#    #data["compressedieta"][0] = int(math.copysign(compressedIeta[abs(data["ieta"][0])], data["ieta"][0]))
#    data["compressedieta"][0]  = compressedIeta[abs(data["ieta"][0])]
#    data["compressedE"][0]     = compressedE[min(data["l1RawE"][0],255)]
#    #data["compressedNTT"][0] = data["nTT"][0]
#    data["compressedNTT"][0] = min(compressedNTT[data["nTT"][0]],11)
#    outputTree.Fill()

outputTree.Write()
outputFile.Close()
inputFile.Close()
