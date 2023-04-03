#!/usr/bin/env python
# coding: utf-8



import uproot as urt
import json,os,argparse
import numpy as np
hep.cms.label("Work In Progress",data=True,com=13.6,year=2022)
def main():
 
    parser = argparse.ArgumentParser()
    parser.add_argument('-i',"--inputFile", help="Input File",default='workarea/TagAndProbe_ReEmul_EGEraG_Remul_caloParam_v0_6_recalib_v0p2.root')
    parser.add_argument("-o","--dest", help="destination To Use", default='workarea/results/dataAugmentation/' )
    parser.add_argument("-t","--tag", help="Tag for the ", default='' )
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
    
    #for PU_MIN,PU_MAX in zip( [0 ,40,50,60,70,80] , [ 40 , 50 , 60 , 70 , 80 , 1999]):
    if True:    
        safetyMask = baseDSet['l1tEmuRawEt'] > 0.0

        baseSet={ ky : baseDSet[ky][safetyMask]  for ky in baseDSet }
        
        
        # #### Checking basic histograms
        
        
        #f,ax=plt.subplots(2,2,figsize=(16,14.5))
        #ax=np.ndarray.flatten(ax)
        #_=ax[0].hist(baseSet['l1tEmuRawEt'],bins=np.arange(-1.5,200.5,1))
        #ax[0].set_xlabel('l1tEmuRawEt')  ; ax[0].set_ylabel('Events')
        #_=ax[1].hist(baseSet['l1tEmuNTT']  ,bins=np.arange(-0.5,200.5,1))
        #ax[1].set_xlabel('l1tEmuNTT')  ; ax[1].set_ylabel('Events')
        #_=ax[2].hist(baseSet['l1tEmuTowerIEta'],bins=np.arange(-35.5,35.5,1))
        #ax[2].set_xlabel('l1tEmuTowerIEta')  ; ax[2].set_ylabel('Events')
        #_=ax[3].hist(baseSet['l1tEmuIsoEt']  ,bins=np.arange(-0.5,40.5,1))
        #ax[3].set_xlabel('l1tEmuIsoEt')  ; ax[3].set_ylabel('Events')
        #f.savefig(prefixBase+'/default_l1Histograms.png',bbox_inches='tight')
        
        
        if args.export :       
            export_filename=prefixBase+"/slimmedDSet_"+args.tag+".root"
            export_file = urt.recreate(export_filename)
            print("\n Base Input file  exported as  : ",export_filename)
            allBranchesToWrite={ ky : baseSet[ky].dtype for ky in baseSet}
            export_file.mktree("Ntuplizer/TagAndProbe",allBranchesToWrite)
            export_file["Ntuplizer/TagAndProbe"].extend(baseSet)
            export_file.close()

            print("\n","- -"*20)
    
if __name__=='__main__':
    main( )

    

