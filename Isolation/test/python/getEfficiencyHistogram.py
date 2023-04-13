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
    parser.add_argument("--opt", help="Option name ",default='ISO_LUT_Progression_v2_2_5p0_0p1_7p5')
    args = parser.parse_args()


    fileIn=uproot.open(args.input)
    data=fileIn['TagAndProbe']

    dataStore=data.arrays(['l1tEmuPt','l1tEmuNTT','l1tEmuRawEt',
                       'eleProbeSclEt','isProbeLoose',
                       'eleProbeEta','eleProbePhi',
                       'eleTagEta','eleTagPhi',args.opt])


    dataStore['dR']= np.sqrt( 
                             ( dataStore['eleProbeEta'] - dataStore['eleTagEta'] )**2 + 
                             ( dataStore['eleProbePhi'] - dataStore['eleTagPhi'] )**2 
                            ) 
    
    preselection_mask = dataStore['eleProbeEta'] < 2.5
    preselection_mask = np.logical_and( preselection_mask , dataStore['dR'] > 0.6 )
    preselection_mask = np.logical_and( preselection_mask , dataStore['isProbeLoose'] > 0.5 ) 
    
    
    pass_mask= np.logical_and(preselection_mask , dataStore['l1tEmuPt'] > 25 )
    pass_mask= np.logical_and(pass_mask         , dataStore[args.opt] > 0.5    )
    
    
    pT_all = dataStore['eleProbeSclEt'][preselection_mask]
    pT_pass= dataStore['eleProbeSclEt'][pass_mask]
    
    
    
    pTBins= np.arange(0,100,1)
    pTBinCenters=0.5*(pTBins[:-1]+pTBins[1:])
    count_all ,_=np.histogram(pT_all ,bins=pTBins)
    count_pass,_=np.histogram(pT_pass,bins=pTBins)
    
    eff=getEff(count_pass,count_all)
    
    f=plt.figure(figsize=(3.5,4))
    plt.scatter(pTBinCenters,eff,c='b')
    hep.cms.label('internal',data=True,year='2022, EraG',cm=13.6)
    plt.ylim([-0.01,1.02])
    plt.ylabel('$\epsilon$')
    plt.xlabel('E$_{T}$')
    plt.grid(which='both')
    fname=args.dest+'/'+args.tag+'effVsPt.png'
    print("Output plot saved as  ",fname)
    plt.savefig(fname,bbox_inches='tight')


if __name__=='__main__':
    main()
