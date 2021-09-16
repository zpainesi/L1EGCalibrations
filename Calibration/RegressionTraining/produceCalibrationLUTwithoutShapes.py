import ROOT
from GBR2LUT import GBR2LUT
from GBR2LUTEmulator_test_newcal import GBR2LUTEmulator_test_newcal

inputFile = "/grid_mnt/t3storage3/athachay/l1egamma/emulationstuff/CMSSW_7_6_0/src/EG_Calibrations/L1EGCalibrations/run3MCL1EtRegression_results.root"
version = "v17.04.04"
sortedShapesFile = "data/compressedSortedShapesForMakingLUT.txt"

header = """# Calibration vs |ieta|,shape,E. Derived from Run 283478 data, with semi-parametric regression
# The LUT output is (ET_off/ET_L1) between 0 and 2, encoded on 9 bits
# Index is compressedShape+compressedE<<4+compressedIeta<<8. 
# Compression version is v4
#anything after # is ignored with the exception of the header
#the header is first valid line starting with #<header> versionStr(unused but may be in future) nrBitsAddress nrBitsData </header>
#<header> V8 12 10 </header>
"""

ietapoints = []
shapepoints = []
Epoints = []
for i in range(0,16):
    ietapoints.append(i)

for i in range(0,16):
    shapepoints.append(i)

for i in range(0,16):
    Epoints.append(i)

gbr2luts = []

gbr2luts.append(GBR2LUTEmulator_test_newcal())
gbr2luts[-1].name = "Run3MCRegression_compressedsortedshape_results"
gbr2luts[-1].inputFileName = inputFile
gbr2luts[-1].outputFileName = "Run3MCRegression_corrections_Trimming10_compressedieta_compressedE_compressedshape_mode_PANTELIS_v2_{}.txt".format(version)
gbr2luts[-1].sortedShapes = sortedShapesFile
gbr2luts[-1].variablePoints.append(("abs(compressedieta)",ietapoints))
gbr2luts[-1].variablePoints.append(("compressedE",Epoints))
gbr2luts[-1].variablePoints.append(("compressedsortedshape",shapepoints))
gbr2luts[-1].retrieveForest()
gbr2luts[-1].createLUT(header)


