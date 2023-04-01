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
    parser.add_argument("--noAug", help="dont augment dataset",default=False,action='store_true')
    parser.add_argument("--export", help="export the augmented dataset",default=False,action='store_true')
    parser.add_argument("--exportOnlyBase", help="export the base dataset",default=False,action='store_true')

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
    
    
    
    baseSet=tAndPTree.arrays(branchList, library="np")
    
    safetyMask= baseSet['l1tEmuRawEt'] > 0.0
    
    baseSet={ ky : baseSet[ky][safetyMask]  for ky in baseSet }
    
    
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
    
    
    
    
    #### Bin edges for the nTT binnig
    
    
    
    """
    binEdges=np.arange(40,200,20)
    counts,edg=np.histogram(baseSet['l1tEmuNTT'],bins=binEdges)
    nTargetPerBin=max(counts)
    print("Bin Edges  : ",binEdges)
    print("Counts per bin : ",counts)
    print("events to have per bin : ",nTargetPerBin)
    """
    binEdgesForNTTScan=np.arange(-0.5,199.5,1)
    countsNTTScan,_=np.histogram(baseSet['l1tEmuNTT'],bins=binEdgesForNTTScan)
    nTargetPerNTT=max(countsNTTScan)
    
    binSize=20
    binEdges=np.arange(-0.5,199.5,binSize)
    counts,edg=np.histogram(baseSet['l1tEmuNTT'],bins=binEdges)
    nTargetPerBin=nTargetPerNTT*binSize
    print("Bin Edges  : ",binEdges)
    print("Counts per bin : ",counts)
    print("events to have per bin : ",nTargetPerBin)
    
    nTTMin=binEdges[4]
    nTTMax=binEdges[5]
    print("Making dataset for nTT bin [  ",nTTMin,nTTMax," ]")
    
    for bidx in range(binEdges.shape[0]-1):
        if args.noAug:
            break
        nTTMin=binEdges[bidx]
        nTTMax=binEdges[bidx+1]
        saveBase=prefixBase+'/nTTBin_'+str(nTTMin)+"To"+str(nTTMax)
        os.system('mkdir -p '+saveBase)
        print()
        print(" Making dataset for nTT bin [  ",nTTMin,nTTMax," ]")
        
        
        # #### Making the `dataToSampleFrom` and `existingData`
        
        
        
        baseMaskInBin =np.logical_and( baseSet['l1tEmuNTT'] < nTTMax , baseSet['l1tEmuNTT'] >= nTTMin ) 
        existingData={ ky : baseSet[ky][baseMaskInBin] for ky in baseSet } 
        
        baseMaskNotInBin = np.logical_not(baseMaskInBin)
        dataToSampleFrom={ky : baseSet[ky][baseMaskNotInBin] for ky in baseSet}
        
        nEventsToSample = nTargetPerBin - np.sum(baseMaskInBin)
        nEventsInBaseToSampleFrom= np.sum(baseMaskNotInBin)
        print("\t Number of events in  existingData : ",np.sum(baseMaskInBin))
        print("\t Number of events to be sampled : ",nEventsToSample)
        print("\t Number of events in dataToSampleFrom from   : ",nEventsInBaseToSampleFrom)
        
        
        # #### Deriving the PDF of nTT distribution to sample from
        
        nEvts,nTTs,_=plt.hist(existingData['l1tEmuNTT'],bins=np.arange(nTTMin,nTTMax+1,1))
        nTTs=0.5*(nTTs[:-1]+nTTs[1:])  
        plt.title(f'NTT in [{nTTMin},{nTTMax})')
        plt.xlabel('nTT')
        plt.ylabel('Events')
        plt.savefig(saveBase+'/nTTinBin.png',bbox_inches='tight')
        plt.close()
        
        print(f"\t Deriving target nTT PDF in the bin [ {nTTMin} , {nTTMax} ]")
        
        nEventsToSamplePerBin=int(nTargetPerBin/nEvts.shape[0])
        eventCountToSample = np.ones(nEvts.shape[0],int)*int(nEventsToSamplePerBin)-nEvts # number of events to be sampled removing the existing events
        
        eventCountToSample[eventCountToSample < 0 ]=0  # don't sample events with nTT where we already have statistics
        eventCountToSample=(eventCountToSample)
        
        probaPerBin=eventCountToSample/np.sum(eventCountToSample)
        #print("Derived PDF : \n\t",end="")
        #for i,j in zip(nTTs,probaPerBin):
        #    print(f"({i},{j:.3f})  ,  ",end="")
        print("\n\t sum of the values in pdf : ",np.sum(probaPerBin))
        
        # #### Drawing `nTT` from the derived distribution
        
        print("\t Drawing augmented_nTTVals from derived pdf !")
        augmented_nTTVals =  np.random.choice(nTTs,size=nEventsToSample,p=probaPerBin)
        
        _=plt.hist(augmented_nTTVals,bins=np.arange(nTTMin,nTTMax,1))
        plt.title(f'NTT in [{nTTMin},{nTTMax})')
        plt.xlabel('nTT')
        plt.ylabel('Events')
        plt.savefig(saveBase+'/ntt_augumented.png',bbox_inches='tight')
        plt.close()
        
        
        # #### Drawing the events from the sample distribution
        
        
        print("\t Drawing events from dataToSampleFrom !")
        sample_idx=np.random.choice(nEventsInBaseToSampleFrom ,size=nEventsToSample,replace= False )
        augmented_dataSet={ky : dataToSampleFrom[ky][sample_idx] for ky in dataToSampleFrom}
        
        
        
        
        _=plt.hist(augmented_dataSet['l1tEmuNTT'],bins=np.arange(-0.5,200,1),density=False)
        plt.title(f'NTT in [{nTTMin},{nTTMax})')
        plt.xlabel('nTT')
        plt.savefig(saveBase+'/originalNTT_inArgumentedSet.png',bbox_inches='tight')
        plt.close()
        
        
        # #### Obtaining the IsoEt as per the augumented nTT Distibution
        
        # The variation of mean IsoEt vs  NTT stays same even with the noise levels increasing. We can use this to estimate the nTT distribution for the events with augumented nTT 
        # 
        # ![image.png](attachment:image.png)
        # 
        # 
        # **The mean of the slope in era G is `0.0946`. We can offset the original IsoEt in the event using this slope and the difference of nTT in the event and the augumented_ntt**
        
        
        print("\t Offsetting the IsoEt to the augumented ntt values")
        
        dNTT =  augmented_nTTVals - augmented_dataSet['l1tEmuNTT']
        augmented_isoEts = np.round( augmented_dataSet['l1tEmuIsoEt'] + dNTT* 0.0946)
        augmented_dataSet['l1tEmuNTT']   = np.array(augmented_nTTVals, np.int32)
        augmented_dataSet['l1tEmuIsoEt'] = np.array(augmented_isoEts, np.int32)
        
        _=plt.hist(dNTT,bins=np.arange(-200.5,200,1),density=False)
        plt.title(f'NTT in [{nTTMin},{nTTMax})')
        plt.xlabel('nTT_orig - nTT_augmented')
        plt.savefig(saveBase+'/dNTT_forArgumentedSet.png',bbox_inches='tight')
        plt.close()
        
        
        # #### Merge the `augmented_dataSet` and `existingData`
        
        
        print("\t Making the final augmented dataSet !")
        
        final_dataset={ky : np.concatenate([existingData[ky],augmented_dataSet[ky]]) for ky in augmented_dataSet}
        
        
        
        
        
        
        print("\t Making Comparison plots !")
        
        fd_mean=np.average(final_dataset['l1tEmuIsoEt'])
        fd_rms= np.std(final_dataset['l1tEmuIsoEt'])
        print(f"\t Mean [rms] of IsoEt Augmented dataset : {fd_mean:.3f} [{fd_rms:.3f}]")
        
        ed_mean=np.average(existingData['l1tEmuIsoEt'])
        ed_rms= np.std(existingData['l1tEmuIsoEt'])
        print(f"\t Mean [rms] of IsoEt Original dataset : {ed_mean:.3f},[{ed_rms:.3f}]")
        
        # _=plt.hist(augmented_isoEts,bins=np.arange(0.0,50.0,1),histtype='step'  , density=True , label ='AUG')
        _=plt.hist(final_dataset['l1tEmuIsoEt'],bins=np.arange(0.0,50.0,1),histtype='step'  , density=True , label ='AUG+Data')
        _=plt.hist(existingData['l1tEmuIsoEt'],bins=np.arange(0.0,50.0,1),histtype='step', density=True, label = 'DATA')
        plt.title(f'NTT in [{nTTMin},{nTTMax})')
        plt.xlabel('IsoEt')
        plt.annotate('AUG+DATA <IsoEt> : '+str(np.round(fd_mean,3))  ,fontsize=15,xy=(0.5,0.6),xycoords='axes fraction')
        plt.annotate('DATA <IsoEt> : '+str(np.round(ed_mean,3))      ,fontsize=15,xy=(0.5,0.55),xycoords='axes fraction')
        plt.annotate('AUG+DATA rms(IsoEt) : '+str(np.round(fd_rms,3)),fontsize=16,xy=(0.5,0.45),xycoords='axes fraction')
        plt.annotate('DATA rms(IsoEt) : '+str(np.round(ed_rms,3))    ,fontsize=16,xy=(0.5,0.40),xycoords='axes fraction')
        plt.legend()
        
        plt.savefig(saveBase+'/isoEt_comparison.png',bbox_inches='tight')
        plt.close()
        
        fd_mean=np.average(final_dataset['l1tEmuIsoEt']/final_dataset['l1tEmuNTT'])
        fd_rms= np.std(final_dataset['l1tEmuIsoEt']/final_dataset['l1tEmuNTT'])
        print(f"\t Mean [rms] of IsoEt/nTT Augmented dataset : {fd_mean:.3f},{fd_rms:.3f}")
        
        ed_mean=np.average(existingData['l1tEmuIsoEt'])
        ed_rms= np.std(existingData['l1tEmuIsoEt'])
        print(f"\t Mean [rms] of IsoEt/nTT Original dataset : {ed_mean:.3f},{ed_rms:.3f}")
        
        # _=plt.hist(augmented_isoEts,bins=np.arange(0.0,50.0,1),histtype='step'  , density=True , label ='AUG')
        _=plt.hist(final_dataset['l1tEmuIsoEt']/final_dataset['l1tEmuNTT'],bins=np.arange(0.0,1.0,0.01),histtype='step'  , density=False , label ='AUG+Data')
        _=plt.hist(existingData['l1tEmuIsoEt']/existingData['l1tEmuNTT'],bins=np.arange(0.0,1.0,0.01),histtype='step', density=False, label = 'DATA')
        plt.title(f'NTT in [{nTTMin},{nTTMax})')
        plt.xlabel('IsoEt')
        plt.annotate('AUG+DATA <IsoEt/Ntt> : '+str(np.round(fd_mean,3))  ,fontsize=15,xy=(0.5,0.6),xycoords='axes fraction')
        plt.annotate('DATA <IsoEt/NTT> : '+str(np.round(ed_mean,3))      ,fontsize=15,xy=(0.5,0.55),xycoords='axes fraction')
        plt.annotate('AUG+DATA rms(IsoEt/NTT) : '+str(np.round(fd_rms,3)),fontsize=16,xy=(0.5,0.45),xycoords='axes fraction')
        plt.annotate('DATA rms(IsoEt/NTT) : '+str(np.round(ed_rms,3))    ,fontsize=16,xy=(0.5,0.40),xycoords='axes fraction')
        plt.legend()
        
        plt.savefig(saveBase+'/isoEtOverNTT_comparison.png',bbox_inches='tight')
        plt.close()
        
        
        
        
        # _=plt.hist(augmented_dataSet['l1tEmuNTT'],bins=np.arange(0.0,200.0,1),histtype='step'  , density=True , label ='AUG')
        _=plt.hist(final_dataset['l1tEmuNTT'],bins=np.arange(0.0,200.0,1),histtype='step'  , density=False , label ='AUG+Data')
        _=plt.hist(existingData['l1tEmuNTT'],bins=np.arange(0.0,200.0,1),histtype='step', density=False, label = 'DATA')
        plt.title(f'NTT in [{nTTMin},{nTTMax})')
        plt.xlabel('nTT')
        plt.legend()
        
        plt.savefig(saveBase+'/nTT_comparison.png',bbox_inches='tight')
        plt.close()
        
        
        
        
        # _=plt.hist(augmented_dataSet['l1tEmuNTT'],bins=np.arange(0.0,200.0,1),histtype='step'  , density=True , label ='AUG')
        _=plt.hist(final_dataset['l1tEmuRawEt'],bins=np.arange(0.0,200.0,1),histtype='step'  , density=True , label ='AUG+Data')
        _=plt.hist(existingData['l1tEmuRawEt'],bins=np.arange(0.0,200.0,1),histtype='step', density=True, label = 'DATA')
        plt.title(f'NTT in [{nTTMin},{nTTMax})')
        plt.xlabel('l1tEmuRawEt')
        plt.legend()
        
        plt.savefig(saveBase+'/rawEt_comparison.png',bbox_inches='tight')
        plt.close()
        
        
        
        
        # _=plt.hist(augmented_dataSet['l1tEmuNTT'],bins=np.arange(0.0,200.0,1),histtype='step'  , density=True , label ='AUG')
        _=plt.hist(final_dataset['l1tEmuTowerIEta'],bins=np.arange(-40.0,40.0,1),histtype='step'  , density=True , label ='AUG+Data')
        _=plt.hist(existingData['l1tEmuTowerIEta'],bins=np.arange(-40.0,40.0,1),histtype='step', density=True, label = 'DATA')
        plt.title(f'NTT in [{nTTMin},{nTTMax})')
        plt.xlabel('l1tEmuTowerIEta')
        plt.legend()
        
        plt.savefig(saveBase+'/towerEta_comparison.png',bbox_inches='tight')
        plt.close()
        
        
        # ## Exporting the augmented dataset
        if args.export:       
            export_filename=saveBase+"/augDSet.root"
            export_file = urt.recreate(export_filename)
            print("\n\t Augmented file exported as  : ",export_filename)
            allBranchesToWrite={ ky : final_dataset[ky].dtype for ky in final_dataset}
            export_file.mktree("Ntuplizer/TagAndProbe",allBranchesToWrite)
            export_file["Ntuplizer/TagAndProbe"].extend(final_dataset)
            export_file.close()
    
        print("\n","- -"*20)
    if args.export or args.xportOnlyBase:       
        export_filename=prefixBase+"/baseDSet.root"
        export_file = urt.recreate(export_filename)
        print("\n Base Input file file exported as  : ",export_filename)
        allBranchesToWrite={ ky : baseSet[ky].dtype for ky in baseSet}
        export_file.mktree("Ntuplizer/TagAndProbe",allBranchesToWrite)
        export_file["Ntuplizer/TagAndProbe"].extend(baseSet)
        export_file.close()

        print("\n","- -"*20)
    
if __name__=='__main__':
    main( )

    

