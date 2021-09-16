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
girt clone git@github.com:ats2008/L1EGCalibrations.git
scram b -j 8
cd L1EGCalibrations/Calibration/RegressionTraining/
make -j 4  # To generate regression.exe file

```
## Preprocessing of data files for training
1. Make the reduced tree with relevant branches. Modify  `L1EGCalibrations/Calibration/RegressionTraining/makeRegressionTree.cc` to point to the TagAndProbe Ntuple toi train on.
   ```bash
    $ root makeRegressionTree.cc
   ```

   This step will produce a new root file with the L1 quantities needed for the regression training. These variables stored in it are:
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

2. Use the script `produceTreeWithCompressedIetaEShape.py` to compress the variables. This file uses the following .txt files containing the compressions of each variable:
    * data/compressedSortedShapes.txt
    * data/egCompressShapesLUT_calibr_4bit_v4.txt
    * data/egCompressELUT_4bit_v4.txt
    * data/egCompressEtaLUT_4bit_v4.txt
    
    This step will produce a new root file  which will be used for regression of caliberation LUT.

3. Setup the `regressionConfig.cfg` for regression training 

4. Derive the calibrations:
```bash
$ ./regression.exe regressionConfig.cfg
```
    This step should make a *result* in the designated output directory in `regressionConfig.cfg`

5. Use `produceCalibrationLUTwithoutShapes.py` to 
```bash
$ python produceCalibrationLUTwithoutShapes.py    # Produce the calibration LUT without the shapes:
```

6. To produce the calibration LUT with the shapes use the script `makeCorrectedLUT.cc` add the shape according to nTT, Et and Eta. 
```
$ root makeCorrectedLUT.cc
```

Compare the two LUTs, they should start from the same calibration value.

Add some information on the output LUT file, similar to:

```
    # Calibration vs |ieta|,shape,E. Derived from Run 283478 data, with semi-parametric regression
    # The LUT output is (ET_off/ET_L1) between 0 and 2, encoded on 9 bits
    # Index is compressedShape+compressedE<<4+compressedIeta<<8.
    # Compression version is v4 
    #anything after # is ignored with the exception of the header
    #the header is first valid line starting with #<header> versionStr(unused but may be in future) nrBitsAddress nrBitsData </header>
    #<header> V8 12 10 </header>
```

Once you check that everything looks OK in the LUT with shapes, use it to produce a new emulated sample, and compare efficiencies and resolutions with the 
previous emulated sample. 
