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
    <fnames list : file with the list of root files to analyze>
    <chkQuality : 0/1 |  if 0, implies the area is calulated and written to result even if the turnon fails isGoodTurnon() or eT<eT_base+2 >  \
    <saveGraph : 0/1 weather to save the areaGraphs as png , note : the graph is always saved if its better than baseline > \
    <prefix To save : result.txt and the graphs will be saved here>  

#example
# will only write the good turnons to results and graphs
./evalOptions.exe fnamesGridAV0 1 0 results/GridBV0
```

