import uproot
import ROOT
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
import argparse

def getEff(passCounts,failCounts):
    return passCounts/(failCounts +1e-9)



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-t',"--tag", help="Tag for file",default='')
    parser.add_argument('-i',"--input", help="input file")
    parser.add_argument('-d',"--dest", help="output destination prefix ",default='./')
    parser.add_argument("--doDefault", help="plot default Iolation effi. ",action='store_true')
    parser.add_argument("--opt", help="Option name ",default='ISO_LUT_Progression_v2_702_5p0_0p1_7p5')
    parser.add_argument("--L1Et_Threshold", help="Min Value of Offline Et ",default=12.0,type=float)
    parser.add_argument("--OfflineET_Min", help="Min Value of Offline Et ",default=0.0,type=float)
    parser.add_argument("--OfflineET_Max", help="Max Value of Offline Et ",default=1e9,type=float)
    parser.add_argument("--OfflineEta_Min", help="Min Value of Offline Eta ",default= 0.0,type=float)
    parser.add_argument("--OfflineEta_Max", help="Max Value of Offline Eta ",default=+2.5,type=float)
    args = parser.parse_args()
    
    L1Et_Threshold = args.L1Et_Threshold
    OfflineET_Min = args.OfflineET_Min
    OfflineET_Max = args.OfflineET_Max

    OfflineEta_Min = args.OfflineEta_Min
    OfflineEta_Max = args.OfflineEta_Max

    fileIn=uproot.open(args.input)
    data=fileIn['TagAndProbe']

    dataStore=data.arrays(['l1tEmuPt','l1tEmuNTT','l1tEmuRawEt',
                           'eleProbeSclEt','isProbeLoose','l1tEmuIso',
                           'eleProbeEta','eleProbePhi',
                           'eleTagEta','eleTagPhi',args.opt])
    
  #  print("Reading  : ",dataStore)

    dataStore['dR']= np.sqrt( 
                             ( dataStore['eleProbeEta'] - dataStore['eleTagEta'] )**2 + 
                             ( dataStore['eleProbePhi'] - dataStore['eleTagPhi'] )**2 
                            ) 
    
    preselection_mask = dataStore['eleProbeEta'] < 2.5
    preselection_mask = np.logical_and( preselection_mask , dataStore['dR'] > 0.6 )
    preselection_mask = np.logical_and( preselection_mask , dataStore['isProbeLoose'] > 0.5 ) 
    print("Preliminary selection : ",np.sum(preselection_mask))
    ## NTT SPECIFIC
    preselection_mask = np.logical_and( preselection_mask , dataStore['eleProbeSclEt'] > OfflineET_Min ) 
    print(f"selection OfflineET > {OfflineET_Min}  : ",np.sum(preselection_mask))
    preselection_mask = np.logical_and( preselection_mask , dataStore['eleProbeSclEt'] < OfflineET_Max ) 
    print(f"selection OfflineET < {OfflineET_Max}  : ",np.sum(preselection_mask))
    preselection_mask = np.logical_and( preselection_mask , np.abs(dataStore['eleProbeEta']) > OfflineEta_Min ) 
    print(f"selection | OfflineEta | > {OfflineEta_Min}  : ",np.sum(preselection_mask))
    preselection_mask = np.logical_and( preselection_mask , np.abs(dataStore['eleProbeEta']) < OfflineEta_Max ) 
    print(f"selection | OfflineEta | < {OfflineEta_Max}  : ",np.sum(preselection_mask))
    
    pass_mask= np.logical_and(preselection_mask , dataStore['l1tEmuPt']  > L1Et_Threshold     )
    pass_mask= np.logical_and(pass_mask         , dataStore['l1tEmuPt']  > 0.5    )
    iso_mask = dataStore[args.opt]    > 0.5    
    if args.doDefault:
        iso_mask =  dataStore['l1tEmuIso']   > 1.5 
    pass_mask= np.logical_and(pass_mask         , iso_mask  )
    
    
    allEvts = dataStore['l1tEmuNTT'][preselection_mask]
    passEvts= dataStore['l1tEmuNTT'][pass_mask]
    print("Events in after preselection      : ",np.sum(preselection_mask))
    print("Events in after trigger selection : ",np.sum(pass_mask))
    print("Efficiency :  ", np.sum(pass_mask)*1.0/np.sum(preselection_mask)  )
    
    
    Bins= np.arange(0,200,2)
    BinCenters=0.5*(Bins[:-1]+Bins[1:])
    count_all ,_=np.histogram(allEvts ,bins=Bins)
    count_pass,_=np.histogram(passEvts,bins=Bins)
    
    eff=getEff(count_pass,count_all)
    
    f=plt.figure(figsize=(8,4))
    plt.scatter(BinCenters,eff,c='b')
    hep.cms.label('internal',data=True,year='2022, EraG',com=13.6)
    plt.xlim([10,200])
    plt.ylim([-0.02,1.05])
    plt.ylabel('$\epsilon$')
    plt.xlabel('NTT',fontweight='bold',size=15)
    plt.grid(which='both')
    plt.text(50,0.2,"L1 E$_{T}$ >"+f" {L1Et_Threshold} "+"\n"+f"{OfflineET_Min} < "+"E$_{T}^{ Offline }$"+f"< {OfflineET_Max}",
            fontsize=15,fontweight='bold')

    fname=args.dest+'/'+args.tag+'effVsNTT.png'

    print("Output plot saved as  ",fname)
    plt.savefig(fname,bbox_inches='tight')


if __name__=='__main__':
    main()
