import ROOT
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
import json,sys,os,argparse
hep.style.use("CMS")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-t',"--tag", help="Tag for file",default='')
    parser.add_argument('-i',"--input", help="input file")
    parser.add_argument('-d',"--dest", help="output destination prefix ",default='./')
    parser.add_argument("--eta", help="eta bin to make the plot for",default='0')
    parser.add_argument("--et" , help="et  bin to make the plot for",default='8')
    args = parser.parse_args()
    
    inputWPFile = ROOT.TFile(args.input)
    tag=args.tag
    saveBase=args.dest+'/'
    ET=int(args.et)
    ETA=int(args.eta)

    x=np.arange(0,32,1)
    
    for ieta in [ETA]:
        for et in [ET]:
            f,axes=plt.subplots(5,5,figsize=(25,25))
            axes=np.ndarray.flatten(axes)
            axid=0
            pid=0
            for eff in range(1,100+1):
                if axid%10==0:
                    print("Processing Eff ",eff)
                ax=axes[axid]
                proj=f'pz_{eff}_eta{ieta}_e{et}'
                fName='fit_'+proj
    
                hist= inputWPFile.Get('Step1Histos/'+proj)
                fit = inputWPFile.Get('Step1Histos/'+fName)
                a=fit.GetParameter(0)
                b=fit.GetParameter(1)
                chi2=fit.GetChisquare()
                ndf=fit.GetNDF()
                hep.histplot(hist,histtype='errorbar',ax=ax)
                y=a+b*x
                ax.annotate(f'{ eff = }% , { ieta = } ,{ et = }'      , (0.1,0.85), xycoords='axes fraction', fontsize=17)
                ax.annotate(f'Chi2 : {np.round(chi2,2)} , NDF : {ndf}', (0.1,0.7), xycoords='axes fraction',  fontsize=15 )
                ax.annotate(f'{np.round(a,2)} + {np.round(b,2)} *x '  , (0.1,0.6), xycoords='axes fraction',  fontsize=15)
                ax.plot(x,y,c='r')
                axid+=1
                if axid==25 :
                
                    foutname=saveBase+f'{tag}DiagonosticPlotAtEt{et}Eta{ieta}_part{pid}.png'
                    print("Saving figure as  ",foutname)
                    f.savefig(foutname , bbox_inches='tight')
                    plt.close(f)
                    f,axes=plt.subplots(5,5,figsize=(25,25))
                    axes=np.ndarray.flatten(axes)
                    axid=0
                    pid+=1

if __name__=='__main__':
    main()
