CXX 		=g++
LIBS		=`root-config --libs`

CXXFLAGS 	=-g3 -std=c++11  `root-config --cflags` 
INCPATH     = include/

TARGET		= isolationAnalysis
all:	$(TARGET) 

isolationAnalysisSlim : src/IsolationAnalysisSlim.cc 
	$(CXX) $(CXXFLAGS) -I$(INCPATH) -o $@.exe $<  $(LIBS)

isolationAnalysis : src/IsolationAnalysis.cc 
	$(CXX) $(CXXFLAGS) -I$(INCPATH) -o $@.exe $<  $(LIBS)

evalOptions : src/evalOptions.cc 
	$(CXX) $(CXXFLAGS) -I$(INCPATH) -o $@.exe $<  $(LIBS)

applyIsolation : src/ApplyIsolation.cc 
	$(CXX) $(CXXFLAGS) -I$(INCPATH) -o $@.exe $<  $(LIBS)

clean :
	@rm *.exe
