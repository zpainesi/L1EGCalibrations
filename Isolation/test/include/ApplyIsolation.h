//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Dec  1 07:58:18 2021 by ROOT version 6.24/06
// from TTree L1UpgradeTree/L1UpgradeTree
// found on file: L1Ntuple.root
//////////////////////////////////////////////////////////

#ifndef ApplyIsolation_h
#define ApplyIsolation_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <map>
#include <TFile.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH3F.h>
#include <vector>
#include <map>
#include <TGraphAsymmErrors.h>


// Header file for the classes stored in the TTree if any.

#define nBins_fine 100
class ApplyIsolation {
 public :
  TChain          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain
  
  TChain          *fChain1;   //!pointer to the analyzed TTree or TChain                                             
  Int_t           fCurrent1; //!current Tree number in a TChain                                                                                                                                            
  float frate; 
  
  // Declaration of leaf types
  
  //Rate
  UShort_t        nEGs;
  std::vector<float>   egEt;
  std::vector<short>   egBx;
  std::vector<short>   egIsoEt;
  std::vector<short>   egNTT;
  std::vector<short> egRawEt;  
  std::vector<short> egTowerIEta;
    
  std::map<TString,TH1F*> th1fStore;

  //Turn-On
   Float_t         l1tEmuPt;
   Int_t           l1tEmuTowerIEta;
   Int_t           l1tEmuNTT;
   Int_t           l1tEmuRawEt;    
   Float_t         eleProbeSclEt;
   Int_t           l1tEmuIsoEt;
  
   float eleProbeEta	;
   float eleProbePhi	;
   float eleTagEta   ;
   float eleTagPhi	;
   Int_t isProbeLoose;
   
   Int_t reportEvery;
   Long64_t maxEntriesForRate;
   Long64_t maxEntriesForEfficiency;


   ApplyIsolation(std::string& inputFileName);
   ~ApplyIsolation();
   
   void accessTree(std::string & input_filelist_rate, std::string & input_filelist_turn_on);
   void readParameters(const std::string jfile);
   void readTree();
   void loops();
   void bookHistograms();
   void saveHistograms();
   void tokenize(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters);
   
   
   Long64_t nEntries_;
   Long64_t nEntries1_;
   
   std::string ntupleFileNameRate_;
   std::string ntupleFileNameTurnOn_;
   std::string outputFileName_;
   std::string optionsFileName_;    
   std::vector<std::string> lutProgOptVec_;
   

   TFile* optionsFile_;
   TFile* outputFile_;

   bool bookedHistograms_;
   bool nTTRange_;

   void readLUTTable(std::string& file_name, unsigned int& nbin, std::map<short, short>& lut_map);

   UInt_t nBinsIEta;
   UInt_t nBinsIEt;
   UInt_t nBinsNTT;
   UInt_t ET_MAX = 255;
   UInt_t bunches =2544;
   
   std::map<short, short> lutMapEt;
   std::map<short, short> lutMapEta;
   std::map<short, short> lutMapNTT;

   std::map<TString,TH1F*> ptMap_;
   std::map<TString,TH1F*> rateMap_;
   std::map<TString,TH1F*> pt_pass_Map_;
   std::map<TString,TGraphAsymmErrors*> turnOn_Map_;
   TH1F* pT_all;
   std::vector<UInt_t> et_option;

   bool check_pt_rate_dir = false;
   bool check_rate_dir = false;
   bool check_pt_turn_on_dir = false;
   bool check_turn_on_dir = false;                                                                                             

   TDirectory* td;
   TDirectory* td1;
   TDirectory* td2;
   TDirectory* td3;
   Double_t binning[22] = {18,20,22,24,26,28,30,32,35,40,45,50,60,70,90,110,210,350,500,700,1000,2000};
   
   Double_t xEdges_fine[nBins_fine+1];
};

#endif












