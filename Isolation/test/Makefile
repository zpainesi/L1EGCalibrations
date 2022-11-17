CXX 		=g++
LIBS		=`root-config --libs`

CXXFLAGS 	=-g3 -std=c++11  `root-config --cflags` 
INCPATH     = include/

TARGET		= isolationAnalysis
all:	$(TARGET) 

isolationAnalysis : src/IsolationAnalysis.cc 
	$(CXX) $(CXXFLAGS) -I$(INCPATH) -o $@.exe $<  $(LIBS)

evalOptions : src/evalOptions.cc 
	$(CXX) $(CXXFLAGS) -I$(INCPATH) -o $@.exe $<  $(LIBS)

applyIsolation : src/ApplyIsolation.cc 
	$(CXX) $(CXXFLAGS) -I$(INCPATH) -o $@.exe $<  $(LIBS)

L1EGRatePlotter : src/L1EGRatePlotter.cc 
	$(CXX) $(CXXFLAGS) -I$(INCPATH) -o $@.exe $<  $(LIBS)

L1EGRatePlotterData : src/L1EGRatePlotterData.cc 
	$(CXX) $(CXXFLAGS) -I$(INCPATH) -o $@.exe $<  $(LIBS)

L1EGRatePlotterPUReweig : src/L1EGRatePlotterPUReweig.cc
	$(CXX) $(CXXFLAGS) -I$(INCPATH) -o $@.exe $<  $(LIBS)

clean :
	@rm *.exe
