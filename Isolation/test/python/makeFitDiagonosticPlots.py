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
    parser.add_argument("--effMin", help="min efficiency bin to make the plot for",default=1,type=int)
    parser.add_argument("--effMax", help="max efficiency bin to make the plot for",default=100,type=int)
    parser.add_argument('--ietaMin',help="Min value for iEta",default=0,type=int)
    parser.add_argument('--ietMin' ,help="Min value for iEt ",default=0,type=int)
    parser.add_argument('--ietaMax',help="Max value for iEta",default=15,type=int)
    parser.add_argument('--ietMax' ,help="Max value for iEt ",default=15,type=int)
    args = parser.parse_args()
    
    inputWPFile = ROOT.TFile(args.input)
    tag=args.tag
    saveBase=args.dest+'/'
    if not os.path.exists(saveBase):
        os.system('mkdir -p '+saveBase)


    x=np.arange(0,32,1)
    problematicProjectsionsAndFits={'vals':[] , 'fit' : [] ,'projections' : []}

    f,axes=plt.subplots(4,4,figsize=(20,20))
    axes=np.ndarray.flatten(axes)
    axid=0
    count=0
    cMax= (args.effMax - args.effMin+1)*( args.ietaMax - args.ietaMin +1  )*(args.ietMax - args.ietMin +1)
    badFitId=0
    fid=0
    pid=0
    for eff in range(args.effMin,args.effMax+1):
        print("Processing Eff ",eff)
        for ieta in range(args.ietaMin,args.ietaMax+1):
            for et in range(args.ietMin,args.ietMax+1):
                count+=1
                if count%200==0:
                    print(f"\t working on {eff=} {ieta=} , {et=}  {count} / {cMax}")
                fid+=1    
                ax=axes[axid]
                proj=f'pz_{eff}_eta{ieta}_e{et}'
                fName='fit_'+proj
                hist= inputWPFile.Get('Step1Histos/'+proj)
                fit = inputWPFile.Get('Step1Histos/'+fName)
                try:
                    a=fit.GetParameter(0)
                    b=fit.GetParameter(1)
                except:
                    print(proj , hist)
                    print(fName , fit)
                    exit(1)

                if b>0.01:
                    continue
                badFitId+=1
                integral = hist.Integral()
                problematicProjectsionsAndFits['vals'].append( np.array([ eff , ieta, et ]) )
                print(f"\t\t{badFitId = } | Obtained a problematic fit at {eff=} {ieta=} {et=} {integral=} | {a = :.3} {b = :.2}"  )
                chi2=fit.GetChisquare()
                ndf=fit.GetNDF()
                hep.histplot(hist,histtype='errorbar',ax=ax)
                y=a+b*x
                ax.annotate(f'{ eff = }% , { ieta = } ,{ et = }'      , (0.1,0.85), xycoords='axes fraction', fontsize=17)
                ax.annotate(f'Chi2 : {np.round(chi2,2)} , NDF : {ndf}', (0.1,0.7), xycoords='axes fraction',  fontsize=15 )
                ax.annotate(f'{np.round(a,2)} + {np.round(b,2)} *x '  , (0.1,0.6), xycoords='axes fraction',  fontsize=15)
                ax.plot(x,y,c='r')
                axid+=1
                if axid==16:
                    foutname=saveBase+f'{tag}DiagonosticPlot_part{pid}.png'
                    print("Saving figure as  ",foutname)
                    f.savefig(foutname , bbox_inches='tight')
                    plt.close(f)
                    f,axes=plt.subplots(4,4,figsize=(20,20))
                    axes=np.ndarray.flatten(axes)
                    axid=0
                    pid+=1
#                if badFitId > 40:
#                    break    
#            if badFitId > 40:
#                break    
    f.savefig(foutname , bbox_inches='tight')
    axes=np.ndarray.flatten(axes)
    plt.close(f)
    
    f=plt.figure(figsize=(15,12))
    probalematicVals = np.array(problematicProjectsionsAndFits['vals'])
    axe=plt.subplot(2,1,1)
    axe.hist( probalematicVals[:,0] , bins=np.arange(-0.5,99.5 , 1.0 ))
    axe.set_xlabel("Efficiency [ % ]")
    axe.set_ylabel("n bad Fit")
    axe=plt.subplot(2,2,3)
    axe.hist( probalematicVals[:,1] , bins=np.arange(-0.5,args.ietaMax+0.5 , 1.0 ))
    axe.set_xlabel("iEta")
    axe.set_ylabel("n bad Fit")
    axe=plt.subplot(2,2,4)
    axe.hist( probalematicVals[:,2] , bins=np.arange(-0.5,args.ietMax +0.5 , 1.0 ))
    axe.set_xlabel("iEt")
    axe.set_ylabel("n bad Fit")
    foutname=saveBase+f'{tag}DiagonosticPlot_summaryDistribution.png'
    f.savefig(foutname , bbox_inches='tight')

if __name__=='__main__':
    main()
