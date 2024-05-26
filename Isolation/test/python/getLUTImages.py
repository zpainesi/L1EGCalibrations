import ROOT 
import matplotlib.pyplot as plt
import numpy as np
import json,sys,os,argparse
"""
Usage
  python3 -t looeIsoDev_prelims -i EG_Iso_LUT_Flat_WP_582_10p0_0p7_40p0.txt 

"""
parser = argparse.ArgumentParser()
parser.add_argument('-t',"--tag", help="Tag for file, this will become the name of the output png")
parser.add_argument('-i',"--input", help="input LUT file")
args = parser.parse_args()

allLUTS={}
if args.tag:
    allLUTS={args.tag : args.input}

for ky in allLUTS:
    fname=allLUTS[ky]
    print(f"processing {ky} [{fname}]")
    with open(fname) as f :
        txt=f.readlines()

    LUT=np.zeros((16,16,32),int)

    for l in txt:
        if "#"== l[0]:
            continue
        if len(l) < 2 :
            continue
        val=int(l[:-1].split(' ')[1])
        items = l[:-1].split("=")[1].split(":")
        #print(items)
        LUT[int(items[0])][int(items[1])][int(items[2])] = val

    f,axlist=plt.subplots(8,2,figsize=(10,18))
    axlist=np.ndarray.flatten(axlist)
    for et in range(16):
        ax=axlist[et]
        lut_proj=LUT[:,et,:]
    #     f,ax=plt.subplots(figsize=(32/4,16/4))
        maxIsoEt=40
        if et > 8:
            maxIsoEt=120
        lut_proj[lut_proj>maxIsoEt]=maxIsoEt
        c=ax.imshow(lut_proj,cmap='tab20c',origin='lower',vmin=0, vmax=maxIsoEt)
        t=ax.set_xticks(np.arange(0.0,32,4)-0.5,np.arange(0.0,32,4),minor=True)
        t=ax.set_yticks(np.arange(0.0,16,2)-0.5,np.arange(0.0,16,2),minor=True)
        t=ax.set_xticks(np.arange(0.0,32,4)-0.5,np.arange(0.0,32,4))
        t=ax.set_yticks(np.arange(0.0,16,4)-0.5,np.arange(0.0,16,4),color='w')
        ax2 = f.add_axes([ax.get_position().x1+0.005,ax.get_position().y0,0.02,ax.get_position().height])
        ax.grid(color='k',alpha=1,which='both')
        t=ax.text(25.2,12.3,"E$_{T}$ = "+str(et) ,fontweight='bold',c='r',alpha=1)
        t.set_bbox(dict(facecolor='w',alpha=1.0, edgecolor='r'))
        plt.colorbar(c,cax=ax2)
    f.suptitle('LUT : '+ky+"\n"+fname.split('/')[-1],y=0.91)
    fname=ky+'.png'
    print("Saving : ",fname)
    f.savefig(fname,bbox_inches='tight')
    #     if et > 2:
    #         break
