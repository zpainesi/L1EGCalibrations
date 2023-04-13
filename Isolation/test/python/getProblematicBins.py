import ROOT 
import matplotlib.pyplot as plt
import numpy as np
import json,sys,os,argparse


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-i',"--input", help="input file",default=None)
    parser.add_argument("--exportPatch", help="Export the patched LUT",action='store_true')
    parser.add_argument("--printFit", help="Print the corresponding nTT dim in the LUT",action='store_true')
    args = parser.parse_args()
    
    allLUTS={'LUT' : args.input}
    for ky in allLUTS:
        fname=allLUTS[ky]
        print(f"processing {ky} {fname}")
        with open(fname) as f :
            txt=f.readlines()
    
        LUT=np.zeros((16,16,32),int)
        LUTIDX=np.zeros((16,16,32),int)
    
        for l in txt:
            if "#"== l[0]:
                continue
            if len(l) < 2 :
                continue
            lutIdx=int(l[:-1].split(' ')[0])
            val=int(l[:-1].split(' ')[1])
            items = l[:-1].split("=")[1].split(":")
    
            LUT[int(items[0])][int(items[1])][int(items[2])] = val
            LUTIDX[int(items[0])][int(items[1])][int(items[2])] = lutIdx
    
        problematicBINS={'nTT':[],'Et':[]}
        pibinNtt=0
        pibinEt =0
        LUT_Modified=np.copy(LUT)
        problemIDX=[]
        print("Processing : ",allLUTS[ky])
        wrongFits={}
        for ieta in range(16):
            for iEt in range(16):
                for inTT in range(31):
                    if LUT[ieta][iEt][inTT] > LUT[ieta][iEt][inTT+1] :
                        print(f'LUT[{ieta}][{iEt}][{inTT}] > LUT[{ieta}][{iEt}][{inTT+1}]  ie. ',end='')
                        print(f'{LUT[ieta][iEt][inTT]} > {LUT[ieta][iEt][inTT+1]}   [ check {LUTIDX[ieta][iEt][inTT]} ] ')
                        problematicBINS['nTT'].append(
                            np.array( [ ieta, iEt, inTT, LUT[ieta][iEt][inTT],LUT[ieta][iEt][inTT+1] ] )
                        )
                        pibinNtt+=1
                        problemIDX.append(LUTIDX[ieta][iEt][inTT])
                        if args.printFit:
                            fitvals=[]
                            for i in range(32):
                                fitvals.append( LUT[ieta][iEt][i]  )
                            wrongFits.update({(ieta,iEt):fitvals})
        

        if args.printFit:
            count=0
            for ky in wrongFits:
                print(count ,"\t| " ,ky," : ",wrongFits[ky])
                count+=1

        if args.exportPatch:
            f=open('patchExport.txt','w')
            for l in txt:
                if "#"== l[0]:
                    f.write(l)
                    continue
                if len(l) < 2 :
                    f.write(l)
                    continue
                lutIdx=int(l[:-1].split(' ')[0])
                val=int(l[:-1].split(' ')[1])
                cmt=' '.join(l[:-1].split(' ')[2:])
                if lutIdx in problemIDX:
                    f.write(f"{lutIdx} {val-1} # patched {val} --> {val-1} "+ cmt+'\n')
                else:
                    f.write(l)
        
            f.close()
            print("patch saved as patchExport.txt")
            

        print("\t problematic bins : ",pibinNtt)


if __name__=='__main__':
    main()
