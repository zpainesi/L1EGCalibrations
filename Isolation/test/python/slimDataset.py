#!/usr/bin/env python
# coding: utf-8

# # Augumentation Techqnique For Making realistic data at higher nTTs
# ## Idea

# Make an augumented dataset with flat distribution in nTT . This will allow derivation of Isolation thersholds for higher nTT bins
# 
# - Steps
#   - Define `base dataset` [`BD`] as all the events in the re emulated sample with `l1tEmuRawEt > 0`
#   - Bin `BD` in nTT bins [ 20 to 200 in bins of width 20 ]  // kind of arbirary choise as of now
#   - For each bin 
#       - take the data that is already there in `BD` as `existingData` 
#       - Deriven number of events to be drawn to attain the global flat distribution in nTT : 
#       - Define the dataset `dataToSampleFrom` from  as events in `BD` after removing `existingData`
#       - Identify the PDF in NTT that need to be drawn so that the final distribution is flat
#       - Draw `nTT` for the augumented dataset from this distribution : `augmented_nTTVals`
#       - Draw `nEventsToSample` from `dataToSampleFrom`
#       - Derive  IsoEt of the `augmented_dataSet` as 
#           - `augmented_dataSet[ISO ET] =  augmented_dataSet[ISO ET] + slope * ( augmented_dataSet[NTT] -augmented_nTTVals )`
#       - Set `augmented_dataSet[ISO ET] = augmented_nTTVals`
#    - Good to go

# #### Branches required by step 1 and step 2





import matplotlib.pyplot as plt
import mplhep as hep
import uproot as urt
import json,os,argparse
import numpy as np
hep.style.use("CMS")
hep.cms.label("Work In Progress",data=True,com=13.6,year=2022)
def main():
 
    parser = argparse.ArgumentParser()
    parser.add_argument('-i',"--inputFile", help="Input File",default='workarea/TagAndProbe_ReEmul_EGEraG_Remul_caloParam_v0_6_recalib_v0p2.root')
    parser.add_argument("-o","--dest", help="destination To Use", default='workarea/results/dataAugmentation/' )
    parser.add_argument("--export", help="export the augmented dataset",default=False,action='store_true')

    args = parser.parse_args()

    branchList=[
            "l1tEmuPt",
            "l1tEmuNTT",
            "l1tEmuRawEt",
            "l1tEmuTowerIEta",
            "eleProbeSclEt",
            "l1tEmuIsoEt",
            "eleProbeEta",
            "eleProbePhi",
            "eleTagEta",
            "eleTagPhi",
            "isProbeLoose",
            "Nvtx",
            "hasL1Emu_26",
            "hasL1Emu_24",
            "hasL1Emu_22",
            "hasL1Emu_looseiso26",
            "hasL1Emu_looseiso24",
            "hasL1Emu_looseiso22",
            "hasL1Emu_tightiso26",
            "hasL1Emu_tightiso24",
            "hasL1Emu_tightiso22"
        ]
    
    prefixBase= args.dest
    os.system('mkdir -p '+prefixBase)
    
    # #### Open the base datafile here ]
    
    # - We use the EraG files
    inputFile=args.inputFile
    print("Opening Input file : ",inputFile)
    tAndPTree=urt.open(inputFile+':Ntuplizer/TagAndProbe')
    
    
    # #### Reading the branches and defining the safety mask to ignore bad events 
    
    
    
    baseDSet=tAndPTree.arrays(branchList, library="np")
    
    for PU_MIN,PU_MAX in zip( [0 ,40,50,60,70,80] , [ 40 , 50 , 60 , 70 , 80 , 1999]):
        
        safetyMask = baseDSet['l1tEmuRawEt'] > 0.0
        safetyMask = np.logical_and(safetyMask, baseDSet['Nvtx'] >= PU_MIN)
        safetyMask = np.logical_and(safetyMask, baseDSet['Nvtx'] <  PU_MAX)

        baseSet={ ky : baseDSet[ky][safetyMask]  for ky in baseDSet }
        
        
        # #### Checking basic histograms
        
        
        f,ax=plt.subplots(2,2,figsize=(16,14.5))
        ax=np.ndarray.flatten(ax)
        _=ax[0].hist(baseSet['l1tEmuRawEt'],bins=np.arange(-1.5,200.5,1))
        ax[0].set_xlabel('l1tEmuRawEt')  ; ax[0].set_ylabel('Events')
        _=ax[1].hist(baseSet['l1tEmuNTT']  ,bins=np.arange(-0.5,200.5,1))
        ax[1].set_xlabel('l1tEmuNTT')  ; ax[1].set_ylabel('Events')
        _=ax[2].hist(baseSet['l1tEmuTowerIEta'],bins=np.arange(-35.5,35.5,1))
        ax[2].set_xlabel('l1tEmuTowerIEta')  ; ax[2].set_ylabel('Events')
        _=ax[3].hist(baseSet['l1tEmuIsoEt']  ,bins=np.arange(-0.5,40.5,1))
        ax[3].set_xlabel('l1tEmuIsoEt')  ; ax[3].set_ylabel('Events')
        f.savefig(prefixBase+'/default_l1Histograms.png',bbox_inches='tight')
        
        
        if args.export :       
            export_filename=prefixBase+"/baseDSet"+f"_puIn{PU_MIN}To{PU_MAX}"+".root"
            export_file = urt.recreate(export_filename)
            print("\n Base Input file file exported as  : ",export_filename)
            allBranchesToWrite={ ky : baseSet[ky].dtype for ky in baseSet}
            export_file.mktree("Ntuplizer/TagAndProbe",allBranchesToWrite)
            export_file["Ntuplizer/TagAndProbe"].extend(baseSet)
            export_file.close()

            print("\n","- -"*20)
    
if __name__=='__main__':
    main( )

    

