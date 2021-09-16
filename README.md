# L1EG_Calibrations
Clone from charis's L1 Repo
## How to derive the new calibrations recipe:

To produce the new LUTs for the L1 EG Calibrations we need to follow the steps below:
1. Emulation of L1 and making relavant TagAndProbe Ntuples
2. Setting up Regression Training and obtaining  new caliberation
3. Re-Emulating sample with new dataset

### Emulation of L1 to produse TagAndProbeNtuples
    See [EGTagAndProbe](https://github.com/ats2008/EGTagAndProbe) for making Ntuples with emulated data.
    
    Make sure the following variables are stored in the output ntuples : 
    * RunNumber
    * Event Number
    * Lumi
    * l1tEmuTowerIEta
    * l1tEmuRawEt
    * shape
    * eleProbeSclEt
    * TowerHoE
    * l1tEmuNTT
    * isProbeLoose


3. Make sure you have updated the files with the new TPs and/or Layer-1 corrections (L1Trigger/L1TCalorimeter/python/caloParams\_/.cfi.py &  L1Trigger/L1TCalorimeter/data/)


## Setting up the environment to run the regression training:
```
cmsrel CMSSW_7_6_0
cd CMSSW_7_6_0/src
cmsenv
cp -r /afs/cern.ch/work/p/pakontax/public/forMarina/EG_Calibrations/ .
rm -rf EG_Calibrations/HiggsAnalysis  # HiggsAnalysis must be under src
cp -r /afs/cern.ch/work/p/pakontax/public/forMarina/EG_Calibrations/HiggsAnalysis .
scram b -j 16
cd L1EG_Calibrations_modified/RegressionTraining
make -j 16                                                 # To generate regression.exe file
cd ..
```

**Backup location of the files:**  /afs/cern.ch/user/m/mkolosov/public/forCharis/EG_Calibrations

1. **Open EG_Calibrations/L1EG_Calibrations_modified/eletree_maker.C** and update the name of the input file to be used.
2. **Run EG_Calibrations/L1EG_Calibrations_modified/eletree_maker.C** to produce a new root file with the L1 quantities needed for the regression training. These variables are:
    * Run Number
    * Event Number
    * Lumi
    * iEta (l1tEmuTowerIEta from the TagAndProbe tree)
    * E (l1tEmuRawEt from the TagAndProbe tree)
    * Shape (shape from the TagAndProbe tree)
    * Target (target is equal to:    eleProbeSclEt / (0.5 * l1tEmuRawEt))
    * H/E  (TowerHoE from the TagAndProbe tree)
    * nTT  (l1tEmuNTT from the TagAndProbe tree)
    * Run2IdLevel  (isProbeLoose from the TagAndProbe tree)

3. Use the script **produceTreeWithCompressedIetaEShape.py** to compress the variables. This file uses the following .txt files containing the compressions of each variable:
    * compressedSortedShapes.txt
    * egCompressShapesLUT_calibr_4bit_v4.txt
    * egCompressELUT_4bit_v4.txt
    * egCompressEtaLUT_4bit_v4.txt
    
4. Derive the calibrations:


```
cd RegressionTraining     # Create a new config file, similar to GBRLikelihood_Trigger_Stage2_Trimming20_Layer1_Pantelis_compressedInputs.config
./regression.exe GBRLikelihood_Trigger_Stage2_Trimming20_Layer1_Pantelis_compressedInputs.config    # A root file with the results is saved
python makeLUT_Stage2_GBRLikelihood_Trimming_compressed_PU40bx25_Pantelis_test_newcal.py    # Produce the calibration LUT without the shapes:
```


To produce the calibration LUT with the shapes use the script **corrections/write_calibrLUT3.C** to add the shape according to nTT, Et and Eta. 


Compare the two LUTs, they should start from the same calibration value.


Add some information on the output LUT file, similar to:


>  
>  #  
>
>              # Calibration vs |ieta|,shape,E. Derived from Run 283478 data, with semi-parametric regression
>              # The LUT output is (ET_off/ET_L1) between 0 and 2, encoded on 9 bits
>              # Index is compressedShape+compressedE<<4+compressedIeta<<8.
>              # Compression version is v4 
>              #anything after # is ignored with the exception of the header
>              #the header is first valid line starting with #<header> versionStr(unused but may be in future) nrBitsAddress nrBitsData </header>
>              #<header> V8 12 10 </header>



Once you check that everything looks OK in the LUT with shapes, use it to produce a new emulated sample, and compare efficiencies and resolutions with the 
previous emulated sample. 
