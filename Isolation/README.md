# Documentation for Isolation Studies for EGamma @ L1

## TODO
1. Change Var names to EG specific variables

## Work flow
### Step 0 : Making EG - Isolation Study Ntuples
  The EG -Isolation Ntuples are a reduced set of branches from the EG-TagAndProbeNtuples used in the context of L1 Et efficiency/resolution studies.
  The compression schemes used are stored in the `data/*LUT*` files. If it is requiered to have more statistics per bin , one may use custom the LUTs for making the comprssion. The current compression scheme uses the same LUTs used for the EG-Regression LUTs ( eta ans et ) and Tau's nTT LUT.
  ```bash=
  #./produceL1EGIsolationStudyNtuple.py <input file path> <max_events>  <output file name>
  ./produceL1EGIsolationStudyNtuple.py Input_TandPNtuple.root 50000 EGIsolationStudyNtuple.root
  ```
  
### Step 1 :  Computing the isolation cuts
  Use `Build_Isolation_WPs.cc` to make a root file which contains histograms labelled by (i,j,k) corresponding to the value of the isocut for each bin (ieta, iEt, nTT). 
  ```bash=
  egIsoStudy $ root -l -q 'Build_Isolation_WPs.cc("EGIsolationStudyNtuple.root","Iso_LUTs_Distributions.root",1000)'
  ```
### Step 2 : Relaxation of the isocut with pT and filling the _"Options"_
  We relax the isolation cut with pT to reach a 100% plateau efficiency while controlling the trigger rate.  We do the relaxation linearly, starting from a flat efficiency A up to a certain pT B, then relaxing linearly up to a certain pT C, when the 100% efficiency is reached (see picture). We test different isolation options corresponding to different combinations of the values A/B/C (this is what we call Option 1, 2…). 
  ```
  egIsoStudy $ root -l -q 'Fill_Isolation.C("Iso_LUTs_Distributions.root","Iso_LUTs_Options.root")' 
  ```
### Step 3 : Applying the isolation cuts
  Now we can apply the isolation (relaxed) cuts onto your signal MC samples, and get the turnons.We use `ApplyIsolationForTurnOns.C` for this purpose. We can  specify the ET threshold to compute the turnons in this macro. It takes as input your signal MC `EG - Isolation Study Ntuples`  as well as the `Iso_LUTs_Options.root` file.
```bash=
egIsoStudy $ root -l -q 'ApplyIsolationForTurnOns.C("EGIsolationStudyNtuple.root",".Iso_LUTs_Options.root",100000,"Iso_LUTs_Caliberated")
```
