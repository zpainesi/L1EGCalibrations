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
 ```
 # For fixed rate area calulations
 ./evalOptions.exe data/HistgramFile_step3step4_12X.root FixedRateTurnons workarea/optionScanPlots/
 ```
 ```
 # For fixed et area calulations
 ./evalOptions.exe data/HistgramFile_step3step4_12X.root FixedEtTurnons workarea/optionScanPlots/
 ```

