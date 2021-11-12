
#      Calibration vs |ieta|,shape,E
#      Index is compressedShape+compressedE<<4+compressedIeta<<8.
#      The LUT output is (ET_off/ET_L1) between 0 and 2, encoded on 9 bits
#      

def getIndex(cd_sh,cd_E,cd_ieta):
    return cd_sh+ (cd_E<<4) + (cd_ieta<<8)  

def getShapeEnergyIeta(index):
    cd_sh  =( index>>(0) & 0xf   )
    cd_E   =( index>>(4) & 0xf   )
    cd_ieta=( index>>(4+4) & 0xf )
    return cd_sh,cd_E,cd_ieta

def getScaleFactor(scaleLUTvalue):
    return  2.0/2**9 * (scaleLUTvalue & 0x1ff) / 2.0


def readLUT(lutFileName):
    lut = {}
    with open(lutFileName, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if line[0]=='#': continue
            line=line.split('#')[0]
            tokens = line.split()
            if len(tokens)!=2: continue
            lut[int(tokens[0])] = int(tokens[1])
    return lut


