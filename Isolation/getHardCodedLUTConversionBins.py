#!/usr/bin/env python3

compressedIetaFile  = "data/egCompressEtaLUT_4bit_v4.txt"
compressedEFile     = "data/egCompressELUT_4bit_v4.txt"
compressedNTTFile   = "data/egCompressNTTLUT_4bit_v4.txt"

def readLUT(lutFileName,superCompressionFactor=1):
    lut = {}
    with open(lutFileName, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if line[0]=='#': continue
            tokens = line.split()
            if len(tokens)!=2: continue
            lut[int(tokens[0])] = int(int(tokens[1])/superCompressionFactor)
    return lut

## Read Ieta,E,shape compression mapping
compressedIEta  = readLUT(compressedIetaFile,8)
compressedE     = readLUT(compressedEFile,1)
compressedNTT   = readLUT(compressedNTTFile,2)

prevVal =-1
nBin=0
strOut="const Int_t hardcodednTTBins[NbinsnTT] = "
for idx in compressedE:
    if prevVal != compressedE[idx]:
        if nBin==0:
            strOut+='{'
        else :
            strOut+=' , '
        strOut+=str(idx)
        prevVal=compressedE[idx]
        nBin+=1;
strOut+=" };"
print("const Int_t NbinsIEt = ",nBin," ; ")
print(strOut)

prevVal =-1
nBin=0
strOut="const Int_t hardcodedIetBins[NbinsIEt]  = "
for idx in compressedNTT:
    if prevVal != compressedNTT[idx]:
        if nBin==0:
            strOut+='{'
        else :
            strOut+=' , '
        strOut+=str(idx)
        prevVal=compressedNTT[idx]
        nBin+=1;
strOut+=" };"
print("const Int_t NbinsnTT = ",nBin," ; ")
print(strOut)

prevVal =-1
nBin=0
strOut="const Int_t hardcodedIetaBins[NbinsIEta] = " 
print(compressedIEta)
for idx in compressedIEta:
    if prevVal != compressedIEta[idx]:
        if nBin==0:
            strOut+='{'
        else :
            strOut+=' , '
        strOut+=str(idx)
        prevVal=compressedIEta[idx]
        nBin+=1;
strOut+=" };"
print("const Int_t NbinsIEta = ",nBin," ; ")
print(strOut)

