import ROOT 
import matplotlib.pyplot as plt
import numpy as np

allLUTS={
    
    'tau_Iso_2022'   : 'LUTs/Tau_Iso_LUT_2022_calibThr1p7_rate14kHz_V11gs_effMin0p9_G3.txt'
}

for ky in allLUTS:
    fname=allLUTS[ky]
    print(f"processing {ky} [{fname}]")
    with open(fname) as f :
        txt=f.readlines()

    LUT=np.zeros((4,32,32),int)

    for l in txt:
        if "#"== l[0]:
            continue
        if len(l) < 2 :
            continue
        val=int(l[:-1].split(' ')[1])
        items = l[:-1].split("=")[1].split(":")
        LUT[int(items[0])][int(items[1])][int(items[2])] = val

    f,axlist=plt.subplots(2,2,figsize=(16,16))
    axlist=np.ndarray.flatten(axlist)
    for eta in range(4):
        ax=axlist[eta]
        lut_proj=LUT[eta,:,:]
    #     f,ax=plt.subplots(figsize=(32/4,16/4))
        maxIsoEt=100
        #if et > 8:
        #    maxIsoEt=120
        lut_proj[lut_proj>maxIsoEt]=maxIsoEt

        c=ax.imshow(lut_proj,cmap='tab20c',origin='lower',vmin=0, vmax=maxIsoEt)
        t=ax.set_xticks(np.arange(0.0,32,4)-0.5,np.arange(0.0,32,4),minor=True)
        t=ax.set_yticks(np.arange(0.0,32,2)-0.5,np.arange(0.0,32,2),minor=True)
        t=ax.set_xticks(np.arange(0.0,32,4)-0.5,np.arange(0.0,32,4))
        t=ax.set_yticks(np.arange(0.0,32,4)-0.5,np.arange(0.0,32,4),color='w')
        ax2 = f.add_axes([ax.get_position().x1+0.005,ax.get_position().y0,0.02,ax.get_position().height])
        ax.grid(color='k',alpha=1,which='both')
        t=ax.text(25.2,12.3,"$eta$ = "+str(eta) ,fontweight='bold',c='r',alpha=1)
        t.set_bbox(dict(facecolor='w',alpha=1.0, edgecolor='r'))
        plt.colorbar(c,cax=ax2)
    f.suptitle('LUT : '+ky+"\n"+fname.split('/')[-1],y=0.91)
    f.savefig('LUTs/'+ky+'.jpeg',bbox_inches='tight')
    #     if et > 2:
    #         break
