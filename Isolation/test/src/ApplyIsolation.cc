#include "ApplyIsolation.h"
#include<stdlib.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <iostream>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TROOT.h>
#include <sstream>
#include <TBranchElement.h>
#include <fstream>
#include <map>
#include <sstream>
#include <iostream>
#include <string>
#include <TGraphAsymmErrors.h>

ApplyIsolation::ApplyIsolation(std::string& inputFileName){
  maxEntriesRate_=-1;
  maxEntriesEff_=-1;
  readParameters(inputFileName);
  if (ntupleFileNameRate_.size() == 0) {
    std::cout << " Inputfile for rate list missing !!!" << ntupleFileNameRate_ << std::endl;
    return;
  }

  if (ntupleFileNameTurnOn_.size() == 0) {
    std::cout << " Inputfile for turn on list missing !!!" << ntupleFileNameTurnOn_ << std::endl;
    return;
  }

  accessTree(ntupleFileNameRate_ , ntupleFileNameTurnOn_);
  
  assert(fChain);

  bookedHistograms_ = false;

}

ApplyIsolation::~ApplyIsolation() {
}

void ApplyIsolation::accessTree(std::string & input_filelist_rate , std::string & input_filelist_turn_on) {
  std::ifstream myFileRate;
  myFileRate.open(input_filelist_rate.c_str(), std::ios::in);
  if (!myFileRate) {
    std::cout << "Input File: " << input_filelist_rate << " could not be opened!" << std::endl;
    fChain = 0;
  } else {
    fChain = new TChain("L1UpgradeTree");
    static constexpr int BUF_SIZE = 256;
    char buf[BUF_SIZE];
    while (myFileRate.getline(buf, BUF_SIZE, '\n')) {  // Pops off the newline character                                     
      std::string line(buf);
      fChain->AddFile(line.c_str(),-1);
      std::cout << "Adding file " << line << " Entries " << fChain->GetEntries() <<  std::endl;
    }
    //    fChain->Print();                                                                                                             
  }

  std::ifstream myFileTurnOn;
  myFileTurnOn.open(input_filelist_turn_on.c_str(), std::ios::in);

  if (!myFileTurnOn) {
    std::cout << "Input File: " << input_filelist_turn_on << " could not be opened!" << std::endl;
    fChain1 = 0;
  } 
  else {
    fChain1 = new TChain("Ntuplizer/TagAndProbe");
    static constexpr int BUF_SIZE = 256;
    char buf[BUF_SIZE];
    while (myFileTurnOn.getline(buf, BUF_SIZE, '\n')) {  // Pops off the newline character                                
      std::string line(buf);
      fChain1->AddFile(line.c_str(),-1);
      std::cout << "Adding file " << line << " Entries " << fChain1->GetEntries() <<  std::endl;
    }
    //    fChain1->Print();
  }
 }

void ApplyIsolation::loops() {
  short in_compressediEta;
  short in_compressediEt;
  short in_compressedNTT;
  Long64_t nbytes = 0, nb = 0;

  optionsFile_ = new TFile(optionsFileName_.c_str(), "READ");
 
    ////////////////////////////////////////////////////////////////////////////////////////////////
	
    std::map<TString, TH3F*>  ResultProgressionMap ;
    
    for(UInt_t i = 1 ; i < N_OPTIONS ; ++i ) {
      bool Filled_Progression= kFALSE;
      Int_t IsoCut_Progression;
      TString ResultProgressionName_= "LUT_Progression_";
      std::ostringstream convert;
      convert.clear();  
      convert << i;
      ResultProgressionName_ += convert.str();
      ResultProgressionMap[ResultProgressionName_]  = (TH3F*)gDirectory->Get(ResultProgressionName_.Data());
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////


 // optionsFile_->cd("Step2Histos");
  
  //Event loop for rate                                                                                                                                                                                  
  if (fChain == 0) return;
  nEntries_ = fChain->GetEntries();
  std::cout << " Entries " << nEntries_ << std::endl;
  int den=0;

  Long64_t maxEntriesRate = nEntries_;
  if(maxEntriesRate_ > 0)  maxEntriesRate = maxEntriesRate_ < nEntries1_ ? maxEntriesRate_ : nEntries1_;
  
  std::cout<<" Processing "<<maxEntriesRate<<" / "<<nEntries_<<" entries \n";
  for (Long64_t jentry=0; jentry < maxEntriesRate; jentry++) {
    Long64_t ientry = fChain->LoadTree(jentry);
    if(jentry%1000 == 0) std::cout << "Rate Loop Event  # " << jentry << std::endl;
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    
    if (nEGs < 1) continue;
    
    for(UInt_t i = 1 ; i < N_OPTIONS ; ++i ) {
      bool Filled_Progression= kFALSE;
      Int_t IsoCut_Progression;
      TString ResultProgressionName_= "LUT_Progression_";
      TString pt_Progression_= "pt_Progression_";
      std::ostringstream convert;
      convert.clear();  
      convert << i;
      ResultProgressionName_ += convert.str();
      pt_Progression_ += convert.str();
      
      TH3F* ResultProgressionName = ResultProgressionMap[ResultProgressionName_];
      
      for (UShort_t iEG=0; iEG < nEGs; ++iEG) {
	
	if (egBx[iEG]!=0)   continue;
	
	float EG_Et  = egEt[iEG];
	short EG_NTT = egNTT[iEG];
	short EG_IEt = egIEt[iEG];
	short EG_IEta = egIEta[iEG];
	short EG_Iso_Et = egIsoEt[iEG];
	
	std::map<short, short>::iterator EtaPos = lutMapEta.find(abs(EG_IEta));
	if (EtaPos != lutMapEta.end()) in_compressediEta = EtaPos->second;
	else in_compressediEta = nBinsIEta-1;
	
	std::map<short, short>::iterator EtPos = lutMapEt.find(int(EG_Et));
	if (EtPos != lutMapEt.end()) in_compressediEt = EtPos->second;
	else {
	  in_compressediEt = nBinsIEt-1;
	  //	  EG_Et = ET_MAX;
	  std::cout<<EG_IEt<<","<<EG_Et<<std::endl;
	  EG_Et = EG_Et > ET_MAX ? ET_MAX : EG_Et;
	}
	
	std::map<short, short>::iterator NTTPos = lutMapNTT.find(EG_NTT);
	if (NTTPos != lutMapNTT.end()) in_compressedNTT = NTTPos->second;
	else in_compressedNTT = nBinsNTT-1;
        
	IsoCut_Progression = ResultProgressionName->GetBinContent(in_compressediEta+1,in_compressediEt+1,in_compressedNTT+1);
	
	if(!Filled_Progression && EG_Iso_Et<=IsoCut_Progression) {
	  ptMap_[pt_Progression_]->Fill(EG_Et);
	  Filled_Progression = kTRUE;
	}
      }
    }
    den++;
  }
  std::cout<<den<<std::endl;
  double scale = (11.2456 * bunches)/den;

  //Filling Rate Histos  
  for(UInt_t i=0;i< ET_MAX;i++) {
    for(UInt_t j = 1 ; j < N_OPTIONS ; j++ ) {
      TString CurrentNameHisto = "pt_Progression_";
      TString CurrentNameHisto1= "rate_Progression_";
      std::ostringstream convert;
      convert.clear();  
      convert << j;
      CurrentNameHisto += convert.str();
      CurrentNameHisto1 += convert.str();
      rateMap_[CurrentNameHisto1]->SetBinContent(i+1, ptMap_[CurrentNameHisto]->Integral(i+1,ET_MAX)*scale);
    }
  }
  
  //Finding Et for fixed rate (frate)
  int xbin;
  for(UInt_t j = 1 ; j < N_OPTIONS ; j++ ) {
    TString CurrentNameHisto1= "rate_Progression_";
    std::ostringstream convert;
    convert.clear();
    convert << j;
    CurrentNameHisto1 += convert.str();
    float vl,bl; //last value, last bin
    float frate=16.7;
    for(UInt_t i=0;i< ET_MAX;i++) {
      float rate = rateMap_[CurrentNameHisto1]->GetBinContent(i+1);
      if(rate > frate) {
        bl =i;
	vl=rate;
      }
      
      if(rate <=frate){
	xbin = (rate - frate) < abs(rate-vl) ? i : bl;
	break;    
      }
    }     
    std:: cout<<"option: "<<j<<" , "<<"et: "<< xbin<<std::endl;
    et_option.push_back(xbin);
  }

  //  bool check_turn_on_dir = false;
  //booking histograms for fixed Et for all options
  for(UInt_t i = 1 ; i < N_OPTIONS ; ++i ) {
    TString PtPassName_fr_ = "pt_pass_option_";
    TString turnOn_Option_fr_ = "turnOn_Option_";

    std::ostringstream convert;
    convert << i;

    std::ostringstream convert1;
    convert1 << et_option[i-1];

    PtPassName_fr_ += convert.str();
    PtPassName_fr_ +="_Et_";
    PtPassName_fr_=PtPassName_fr_ + convert1.str() + "_fr";
    td2->cd();
    // optionsFile_ = new TFile(optionsFileName_.c_str(), "OPEN");
    //optionsFile_->cd("Step2Histos");    
    TH1F*PtPass_fr= new TH1F(PtPassName_fr_, PtPassName_fr_, 21,binning);
    pt_pass_Map_.insert(std::make_pair(PtPassName_fr_, PtPass_fr));


    turnOn_Option_fr_ += convert.str();
    turnOn_Option_fr_ +="_Et_";
    turnOn_Option_fr_ = turnOn_Option_fr_ + convert1.str() + "_fr";
    
    if(!check_turn_on_dir) {                                                                                                     
      td3 = outputFile_->mkdir("turn_on_progression");                                                                           
      check_turn_on_dir = true;                                                                                                  
    }                                                                                                                            
    td3->cd();
    TGraphAsymmErrors* turnOn_Option_fr;
    turnOn_Map_.insert(std::make_pair(turnOn_Option_fr_,turnOn_Option_fr));
  }
  
  optionsFile_ = new TFile(optionsFileName_.c_str(), "READ");
  optionsFile_->cd("Step2Histos");     
  
  //Event loop for turnon
  if (fChain1 == 0) return;
  nEntries1_ = fChain1->GetEntries();
  std::cout << " Entries " << nEntries1_ << std::endl;
  
  Long64_t maxEntriesEff = nEntries1_;
  if(maxEntriesEff_ > 0)  maxEntriesEff = maxEntriesEff_ < nEntries1_ ? maxEntriesEff_ : nEntries1_;
  
  std::cout<<" Processing "<<maxEntriesEff<<" / "<<nEntries1_<<" entries \n";
  Long64_t nbytes1 = 0, nb1 = 0;
  for (Long64_t jentry=0; jentry< maxEntriesEff; jentry++) {
    Long64_t ientry = fChain1->LoadTree(jentry);
    if(jentry%1000 == 0) std::cout << "Turnon Loop Event  # " << jentry << std::endl;
    if (ientry < 0) break;
    nb1 = fChain1->GetEntry(jentry);   nbytes1 += nb1;
    if(l1tEmuNTT<60 && nTTRange_) continue;  //define range
    if(l1tEmuRawEt < 0.) continue;
    pT_all->Fill(eleProbeSclEt);
    
    std::map<short, short>::iterator EtaPos = lutMapEta.find(abs(l1tEmuTowerIEta));
    if (EtaPos != lutMapEta.end()) in_compressediEta = EtaPos->second;
    else in_compressediEta = nBinsIEta-1;
  
    if(l1tEmuPt < 0.) continue;
    //    if(l1tEmuNTT < 0) continue; 
    std::map<short, short>::iterator EtPos = lutMapEt.find(int(l1tEmuPt));
    if (EtPos != lutMapEt.end()) in_compressediEt = EtPos->second;
    else in_compressediEt = nBinsIEt-1;

    std::map<short, short>::iterator NTTPos = lutMapNTT.find(l1tEmuNTT);
    if (NTTPos != lutMapNTT.end()) in_compressedNTT = NTTPos->second;
    else in_compressedNTT = nBinsNTT-1;
  

    Int_t IsoCut_Progression;
    for(UInt_t e = 20 ; e <= 50 ; e += 5) {

      for(UInt_t i = 1 ; i < N_OPTIONS ; ++i ) {
	TString ResultProgressionName_= "LUT_Progression_";
        TString PtPassName_= "pt_pass_option_";
	
	std::ostringstream convert;
	std::ostringstream convert1;
	//convert.clear();  
	convert << i;
	convert1 << e;
	ResultProgressionName_ += convert.str();
	PtPassName_ += convert.str();
	PtPassName_ +="_Et_";
	PtPassName_ += convert1.str();
	
	TH3F* ResultProgressionName = (TH3F*)gDirectory->Get(ResultProgressionName_.Data());
	IsoCut_Progression = ResultProgressionName->GetBinContent(in_compressediEta+1,in_compressediEt+1,in_compressedNTT+1);
	//	if(i==1 && e==30 && l1tEmuIsoEt >= IsoCut_Progression && eleProbeSclEt>60.)
	// std::cout<<"j: "<<jentry<<" , iso: "<<l1tEmuIsoEt<<" , IsoCut_Progression: "<<IsoCut_Progression<<" , eleProbeSclEt: "<<eleProbeSclEt<<std::endl;
	if(l1tEmuPt >= e && l1tEmuIsoEt <= IsoCut_Progression){
	  pt_pass_Map_[PtPassName_]->Fill(eleProbeSclEt);
	}
      }
    }
    
    for(UInt_t i = 1 ; i < N_OPTIONS ; ++i ) {
      TString PtPassName_fr= "pt_pass_option_";
      TString ResultProgressionName_= "LUT_Progression_";      
      std::ostringstream convert;
      convert << i;
      std::ostringstream convert1;
      convert1 << et_option[i-1];
      
      PtPassName_fr += convert.str();
      PtPassName_fr = PtPassName_fr + "_Et_" + convert1.str() +"_fr";
      ResultProgressionName_ += convert.str();
        //td2->cd();
      //optionsFile_ = new TFile(optionsFileName_.c_str(), "OPEN");
      // optionsFile_->cd("Step2Histos");
      TH3F* ResultProgressionName = (TH3F*)gDirectory->Get(ResultProgressionName_.Data());
      IsoCut_Progression = ResultProgressionName->GetBinContent(in_compressediEta+1,in_compressediEt+1,in_compressedNTT+1);
      if(l1tEmuPt >= et_option[i-1] && l1tEmuIsoEt <= IsoCut_Progression){
	pt_pass_Map_[PtPassName_fr]->Fill(eleProbeSclEt);
      }
    }
  }

  //Turnons for various Ets
  for(UInt_t e = 20 ; e <= 50 ; e += 5) {
    for(UInt_t i = 1 ; i < N_OPTIONS ; i++ ) {
      TString PtPassName_= "pt_pass_option_";
      TString turnOn_Option_="turnOn_Option_";
      std::ostringstream convert;
      convert << i;
      std::ostringstream convert1;
      convert1 << e;
      PtPassName_ += convert.str();
      PtPassName_ +="_Et_";
      PtPassName_ += convert1.str();
      
      turnOn_Option_ += convert.str();
      turnOn_Option_ +="_Et_";
      turnOn_Option_ += convert1.str();
            outputFile_->cd();
      if(!check_turn_on_dir) {
	td3 = outputFile_->mkdir("turn_on_progression");
        check_turn_on_dir = true;
      }
      td3->cd();
      turnOn_Map_[turnOn_Option_]=new TGraphAsymmErrors( pt_pass_Map_[PtPassName_],pT_all,"cp"); 
      turnOn_Map_[turnOn_Option_]->Write();
    }
  }

  //Turnons for fixed Rate
  for(UInt_t i = 1 ; i < N_OPTIONS ; i++ ) {
    TString PtPassName_fr= "pt_pass_option_";
    TString turnOn_Option_fr="turnOn_Option_";
    std::ostringstream convert;
    convert << i;
    std::ostringstream convert1;
    convert1 << et_option[i-1];
    
    PtPassName_fr += convert.str();
    PtPassName_fr = PtPassName_fr + "_Et_" + convert1.str() +"_fr";
    
    td3->cd();
    turnOn_Option_fr += convert.str();
    turnOn_Option_fr = turnOn_Option_fr + "_Et_" + convert1.str() +"_fr";;
    
    turnOn_Map_[turnOn_Option_fr]=new TGraphAsymmErrors( pt_pass_Map_[PtPassName_fr],pT_all,"cp");
    turnOn_Map_[turnOn_Option_fr]->Write();
  }
}
  

void ApplyIsolation::bookHistograms() {
  outputFile_ = new TFile(outputFileName_.c_str(), "RECREATE");
  outputFile_->cd();

  for(UInt_t i = 1 ; i < N_OPTIONS ; ++i ) {
    TString  CurrentNameHisto = "pt_Progression_";
    TString CurrentNameHisto1= "rate_Progression_";
    std::ostringstream convert;
    convert << i;
    CurrentNameHisto += convert.str();
    CurrentNameHisto1 += convert.str();

    if(!check_pt_rate_dir) { 
       td = outputFile_->mkdir("pt_progression_for_rate");
      check_pt_rate_dir = true;
    }
    td->cd();
    TH1F* pt_Progression=  new TH1F(CurrentNameHisto, CurrentNameHisto ,ET_MAX ,  0.0 , ET_MAX);
    ptMap_.insert(std::make_pair(CurrentNameHisto,pt_Progression));

    if(!check_rate_dir) {
       td1 = outputFile_->mkdir("rate_progression");
      check_rate_dir = true;
    }
    td1->cd();   
    TH1F* rate_Progression = new TH1F(CurrentNameHisto1, CurrentNameHisto1, ET_MAX ,  0.0 , ET_MAX );
    rateMap_.insert(std::make_pair(CurrentNameHisto1,rate_Progression));

    for(UInt_t e = 20 ; e <= 50 ; e += 5) {
      TString PtPassName_= "pt_pass_option_";
      TString turnOn_Option_="turnOn_Option_";
      std::ostringstream convert1;
      convert1 << e;    
      PtPassName_ += convert.str();
      PtPassName_ +="_Et_";
      PtPassName_ += convert1.str();  

      if(!check_pt_turn_on_dir) {
	td2 = outputFile_->mkdir("pt_progression_for_turnon");
	check_pt_turn_on_dir = true;
      }
      td2->cd();
      TH1F*PtPass= new TH1F(PtPassName_, PtPassName_, 21,binning);
      pt_pass_Map_.insert(std::make_pair(PtPassName_, PtPass));
      
      turnOn_Option_ += convert.str();
      turnOn_Option_ +="_Et_";
      turnOn_Option_ += convert1.str();
      
      TGraphAsymmErrors* turnOn_Option;
      turnOn_Map_.insert(std::make_pair(turnOn_Option_,turnOn_Option));
    }
  }
  //td2->cd();
  pT_all = new TH1F("pT_all","pT_all",21,binning);
  
  bookedHistograms_ = true;
}

void ApplyIsolation::readLUTTable(std::string& file_name, unsigned int& nbin, std::map<short, short>& lut_map) {
  std::cout << "Opening LUT file " << file_name << std::endl;
   std::ifstream lutFile(file_name.c_str());

  if (!lutFile) {
    std::cerr << "Input File: " << file_name << " could not be opened!" << std::endl;
    return;
  }
  std::string line;
  nbin = 0;
  if(lutFile.is_open()) {
    nbin = 0;
    while(std::getline(lutFile,line)) {
      if (line.substr(0,1) == "#" || line.substr(0,2) == "//") continue;
      if (line.length() == 0) continue;
      std::vector<std::string> tokens;
      tokenize(line,tokens," ");
      unsigned int key   = std::atoi(tokens.at(0).c_str());
      unsigned int value = std::atoi(tokens.at(1).c_str());
      lut_map.insert({ key, value });
      if (nbin < value) nbin = value;
    }
    nbin += 1;
  }
  std::cout << " nbin " << nbin << std::endl;
}

void ApplyIsolation::saveHistograms() {
  if (outputFile_ && bookedHistograms_) {
    outputFile_->Write();
    outputFile_->Close();
  }
}

void ApplyIsolation::readParameters(const std::string jfile) {
  std::cout << jfile << std::endl;
  std::ifstream jobcardFile(jfile.c_str());
   if (!jobcardFile) {
    std::cerr << "Input File: " << jfile << " could not be opened!" << std::endl;
    return;
  }
  std::string line;
  if(jobcardFile.is_open()) {
    while(std::getline(jobcardFile,line)) {
      // enable '#' and '//' style comments                                                                                                                                                   
      if (line.substr(0,1) == "#" || line.substr(0,2) == "//") continue;
      std::vector<std::string> tokens;
      tokenize(line,tokens,"=");
      std::cout << tokens[0] << ":" << tokens[1] << std::endl;
      std::string key   = tokens.at(0);
      std::string value = tokens.at(1);

      if(key=="NtupleFileNameRate")        ntupleFileNameRate_= value;
      else if(key=="NtupleFileNameTurnOn")        ntupleFileNameTurnOn_= value;
      else if (key=="OptionsFileName") optionsFileName_ = value.c_str();
      else if (key=="OutputFileName")  outputFileName_ = value.c_str();
      else if (key=="EtLUTFileName")    readLUTTable(value,nBinsIEt, lutMapEt);
      else if (key=="EtaLUTFileName")   readLUTTable(value,nBinsIEta, lutMapEta);
      else if (key=="ReportEvery")	    {
                reportEvery= atoi(value.c_str());
      }
      else if (key=="MaxEntriesEfficiency")	    {
                maxEntriesEff_= atoi(value.c_str());
      }
      else if (key=="MaxEntriesRate")	    {
                maxEntriesRate_= atoi(value.c_str());
      }
      else if (key=="NTTLUTFileName")   readLUTTable(value,nBinsNTT, lutMapNTT);
      else if (key=="nTTRange") nTTRange_ = std::stoi(value.c_str());
      else
	std::cout << " unknown option " << " key " << key << std::endl;
    }
  }
  jobcardFile.close();
  if(nTTRange_)
    std::cout<<nTTRange_<<"got yaaa"<<std::endl;
}


void ApplyIsolation::readTree() {
  fChain->SetMakeClass(1);
  fChain->SetBranchAddress("nEGs", &nEGs);
  fChain->SetBranchAddress("egEt", &egEt);
  fChain->SetBranchAddress("egIEt", &egIEt);
  fChain->SetBranchAddress("egIEta", &egIEta);
  fChain->SetBranchAddress("egNTT", &egNTT);
  fChain->SetBranchAddress("egBx", &egBx);
  fChain->SetBranchAddress("egIsoEt", &egIsoEt);


  /*
 fChain->SetBranchAddress("nEGs", &nEGs, &b_L1Upgrade_nEGs);
  fChain->SetBranchAddress("egEt", &egEt, &b_L1Upgrade_egEt);
  fChain->SetBranchAddress("egIEt", &egIEt, &b_L1Upgrade_egIEt);
  fChain->SetBranchAddress("egIEta", &egIEta, &b_L1Upgrade_egIEta);
  fChain->SetBranchAddress("egNTT", &egNTT, &b_L1Upgrade_egNTT);
  fChain->SetBranchAddress("egBx", &egBx, &b_L1Upgrade_egBx);
  fChain->SetBranchAddress("egIsoEt", &egIsoEt, &b_L1Upgrade_egIsoEt);

  */

  fChain1->SetBranchAddress("l1tEmuPt", &l1tEmuPt);
  fChain1->SetBranchAddress("l1tEmuEta",&l1tEmuEta);
  fChain1->SetBranchAddress("l1tEmuNTT",&l1tEmuNTT);
  fChain1->SetBranchAddress("l1tEmuRawEt",&l1tEmuRawEt);
  fChain1->SetBranchAddress("l1tEmuTowerIEta",&l1tEmuTowerIEta);
  fChain1->SetBranchAddress("eleProbeSclEt",&eleProbeSclEt);
  fChain1->SetBranchAddress("l1tEmuIsoEt",&l1tEmuIsoEt);

}


void ApplyIsolation::tokenize(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters) {
  // Skip delimiters at beginning.                                                                                          
  std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  
  // Find first "non-delimiter".                                                                                            
  std::string::size_type pos = str.find_first_of(delimiters, lastPos);
  
  while (std::string::npos != pos || std::string::npos != lastPos)  {
    // Found a token, add it to the vector.                                                                      
    tokens.push_back(str.substr(lastPos, pos - lastPos));
    
    // Skip delimiters.  Note the "not_of"                                                                       
    lastPos = str.find_first_not_of(delimiters, pos);
    
    // Find next "non-delimiter"                                                                                    
    pos = str.find_first_of(delimiters, lastPos);
  }
}

int main(int argc,char *argv[]) {
  if (argc == 1) {
    std::cout << " No option provided!!!" << std::endl;
    return 1;
  }
  std::string data_file = argv[1];
  ApplyIsolation treeReader(data_file);
  treeReader.readTree();
  treeReader.bookHistograms();
  std::cout << " Calling Loop" << std::endl;
  treeReader.loops();
  treeReader.saveHistograms();
  return 0;
}


