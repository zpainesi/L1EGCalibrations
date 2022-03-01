## Step 1 and 2 : Making the Working points and filling the Isolations
  * For compiling :
  ```
  make isolationAnalysis
  ```
  
## Step 3 and 4 : Obtaining Rate and Efficincies
  * For compiling :
  ```
  make applyIsolation
  ```

## Evaluating options 
 * For making the plots with area 
 ```
 make evalOptions
 ```
 * Files from step3/4 are can be analysed using the above script
 ```bash
 # For fixed rate area calulations
 ./evalOptions.exe \
    <step3step4 file> \
    FixedRateTurnons \
    <baseline Filename> \
    <baseline TGraph path in file > \
    <prefix to the saved pic>
 ```
 Eg.
 ```bash
 # For fixed rate area calulations
 ./evalOptions.exe \
    caliberationFiles/optionTurnOns/v2_withTheQuality/HistgramFile_step3step4_12X_Scanv2p10.root \
    FixedRateTurnons \
    /grid_mnt/t3storage3/athachay/l1egamma/isolation/CMSSW_12_0_2/src/L1EGCalibrations/Isolation/test/results/run3MCLegacyTurnons/run3MC_12_0_2_Recalib_turnon.root \
    Graph_from_hist_L1Et24_tightisoEfficiency \
    results/optScan/v2_withTheQuality/
 ```

