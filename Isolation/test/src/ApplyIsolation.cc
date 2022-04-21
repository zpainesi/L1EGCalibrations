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
#include <chrono>

ApplyIsolation::ApplyIsolation(std::string& inputFileName){
  maxEntriesForRate=-1;
  maxEntriesForEfficiency=-1;
  reportEvery=5000;
  frate=16.7;
  frateDoubleEG=3.9;
  hasPUReweitingFile=false;
  puReweightingFileName="";
  puReweightingHist=nullptr;
  for(Int_t i=0 ; i<= nBins_fine;i++)
  {
    xEdges_fine[i]=0.0+i*1.0;
  } 

  readParameters(inputFileName);
  if (ntupleFileNameRate_.size() == 0) {
    std::cout << " Inputfile for rate list missing !!!" << ntupleFileNameRate_ << std::endl;
    return;
  }

  if (ntupleFileNameTurnOn_.size() == 0) {
    std::cout << " Inputfile for turn on list missing !!!" << ntupleFileNameTurnOn_ << std::endl;
    return;
  }
  std::cout<<"hasPUReweitingFile : "<<hasPUReweitingFile<<"\n";

  if(hasPUReweitingFile)
  {
    
    auto puHistFile=new TFile(puReweightingFileName,"read");
    if(not puHistFile)
    {
        std::cout<<"PU Rewething file not found : "<<puHistFile;
        exit(2);
    }

    puReweightingHist=(TH1D*) puHistFile->Get("ratio_hist");

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
    fChain = new TChain("l1UpgradeEmuTree/L1UpgradeTree");

    if(hasPUReweitingFile) fChain_evtTree = new TChain("l1EventTree/L1EventTree");

    static constexpr int BUF_SIZE = 256;
    char buf[BUF_SIZE];
    while (myFileRate.getline(buf, BUF_SIZE, '\n')) {  // Pops off the newline character                                     
      std::string line(buf);
      fChain->AddFile(line.c_str(),-1);
      std::cout << "Adding file " << line << " Entries " << fChain->GetEntries() <<"  ";
      if(hasPUReweitingFile) {
        fChain_evtTree->AddFile(line.c_str(),-1);
        std::cout<<"( Evt Tree : "<<fChain_evtTree->GetEntries()<<" ";
      }
      
      std::cout<<std::endl;
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
  optionsFile_->cd("Step2Histos");
  std::map<TString, TH3F*> optionsMap;
    for (auto it :lutProgOptVec_) {
      bool Filled_Progression= kFALSE;
      Int_t IsoCut_Progression;
      TString ResultProgressionName_= "LUT_Progression_";
      TString pt_Progression_= "pt_Progression";
      ResultProgressionName_ +=it;
      pt_Progression_ +=it;
	  TH3F* ResultProgressionName = (TH3F*)optionsFile_->Get("Step2Histos/"+ResultProgressionName_);
      if( not ResultProgressionName) {
        std::cout<<"BAD : "<<ResultProgressionName_<<" !! gonna break later\n exiting";
        exit(1);
      }
      else 
      { 
        std::cout<<" Read option : "<<ResultProgressionName->GetName()<<"\n";
        optionsMap[ResultProgressionName_]=ResultProgressionName;
      }
    }

  //Event loop for rate                                                                                                                                                                                  
  if (fChain == 0) return;
  nEntries_ = fChain->GetEntriesFast();
  if(maxEntriesForRate >0 ) nEntries_= maxEntriesForRate < nEntries_ ? maxEntriesForRate : nEntries_;
  std::cout << " Entries " << nEntries_ << std::endl;
  int den=0;
  Double_t puWeight=1.0;
  auto t_start = std::chrono::high_resolution_clock::now();
  auto t_end = std::chrono::high_resolution_clock::now();
  Double_t nEvents(0);
    
  TH1F eT1ForDoubleEG("et1","et1", ET_MAX , -0.4 , ET_MAX - 0.5 );
  TH1F eTLooseForDoubleEG("etLoose","etLoose",ET_MAX,-0.4, ET_MAX - 0.5);

  for (Long64_t jentry=0; jentry < nEntries_; jentry++) {

    Long64_t ientry = fChain->LoadTree(jentry);
    puWeight=1.0;
    if(hasPUReweitingFile)
    {
        if(nVtx <= puReweightingHist->GetNbinsX()) {
        puWeight=puReweightingHist->GetBinContent(nVtx);
        }
        else
        {
            puWeight=0.0;
        }

        //std::cout<<"nVtx , puWeight : "<<nVtx<<" , "<<puWeight<<"\n";
    }
    nEvents+=puWeight;
 
    if(jentry%reportEvery == 0 )
       {
             t_end = std::chrono::high_resolution_clock::now();
             std::cout<<"Processing Entry in event loop : "<<jentry<<" / "<<nEntries_<<"  [ "<<100.0*jentry/nEntries_<<"  % ]  "
                      << " Elapsed time : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()/1000.0
                      <<"  Estimated time left : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()*( nEntries_ - jentry)/(1e-9 + jentry)* 0.001
                      <<std::endl;
       }


    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if(hasPUReweitingFile){
        nb = fChain_evtTree->GetEntry(jentry);   nbytes += nb;
    }
    if (nEGs < 1) continue;
       // RUN 2 LUT Rates
    float maxEt=-1e9,maxEtLoose=-1e9,maxEtTight=-1e9,maxDoubleEGIsoEt=-1e3;
    eT1ForDoubleEG.Reset();
    eTLooseForDoubleEG.Reset();
    
    for (UShort_t iEG=0; iEG < nEGs; ++iEG) {
	
    if (egBx[iEG]!=0)   continue;
	
	float EG_Et  = egEt[iEG];
	short EG_Iso  = egIso[iEG];
	EG_Et = EG_Et > ET_MAX ? ET_MAX : EG_Et;
    

    if(EG_Et < 0) continue;

    if( EG_Et > maxEt) maxEt=EG_Et;
    if( EG_Et > maxEtLoose and ( EG_Iso==2 or EG_Iso==3) ) maxEtLoose=EG_Et;
    if( EG_Et > maxEtTight and ( EG_Iso==3 or EG_Iso==1) ) maxEtTight=EG_Et;

    eT1ForDoubleEG.Fill( EG_Et );
    if( EG_Iso==2 or EG_Iso==3 )
        eTLooseForDoubleEG.Fill( EG_Et );
    }

    if(maxEt > 0){
        ptMap_["pt_Run2"]->Fill(maxEt,puWeight);
    }
    if(maxEtLoose >= 0.0){
        ptMap_["pt_Run2_loose"]->Fill(maxEtLoose,puWeight);
    }
    if(maxEtTight >=0.0) {
        ptMap_["pt_Run2_tight"]->Fill(maxEtTight,puWeight);
    }
    
    maxDoubleEGIsoEt=-1e3;

    for(UInt_t i=10;i< (ET_MAX ) ;i++) 
    {
           if( eT1ForDoubleEG.Integral(1 + i -10, ET_MAX ) >= 2 ) 
           {
              if(eTLooseForDoubleEG.Integral(1 + i , ET_MAX ) >=1)
              {
                    maxDoubleEGIsoEt=i;
              }
              else
              {
                 break;
              }
           }
           else
           {
                break;
           }
    }

    if( maxDoubleEGIsoEt >= 0.0)
    {
        ptMap_["pt_Run2_DoubleEGLoose_diff10"]->Fill(maxDoubleEGIsoEt,puWeight);
    }

    float maxOptEt;

     // std::cout<<"\nJentry = "<<jentry<<"\n";
    for (auto it :lutProgOptVec_) {
        bool Filled_ProgressionD=false;
      bool Filled_Progression= kFALSE;
      Int_t IsoCut_Progression;
      TString ResultProgressionName_= "LUT_Progression_";
      TString pt_Progression_= "pt_Progression";
      TString pt_DoubleEGProgression_= "pt_DoubleEGProgression";
      ResultProgressionName_ +=it;
      pt_Progression_ +=it;
      pt_DoubleEGProgression_ +=it;
      
      //TH3F* ResultProgressionName = (TH3F*)gDirectory->Get(ResultProgressionName_.Data());
	  //TH3F* ResultProgressionName = (TH3F*)optionsFile_->Get("Step2Histos/"+ResultProgressionName_);
	  TH3F* ResultProgressionName = optionsMap[ResultProgressionName_];
      
      maxOptEt=-1e9;
      eT1ForDoubleEG.Reset();
      eTLooseForDoubleEG.Reset();
      //std::cout<<" opt : "<<pt_DoubleEGProgression_<<" \n";
      for (UShort_t iEG=0; iEG < nEGs; ++iEG) {
	
	if (egBx[iEG]!=0)   continue;
	
	float EG_Et  = egEt[iEG];
	short EG_NTT = egNTT[iEG];
	short EG_TowerIEta = egTowerIEta[iEG];
	short EG_Iso_Et = egIsoEt[iEG];
	short EG_Raw_Et = egRawEt[iEG];
      //std::cout<<"\t"<<iEG<<" : EG Et : "<<EG_Et<<"\n";


	std::map<short, short>::iterator EtaPos = lutMapEta.find(int(abs(EG_TowerIEta)));
	if (EtaPos != lutMapEta.end()) in_compressediEta = EtaPos->second;
	else in_compressediEta = nBinsIEta-1;
	
	std::map<short, short>::iterator EtPos = lutMapEt.find(int(EG_Raw_Et));
	if (EtPos != lutMapEt.end()) in_compressediEt = EtPos->second;
	else {
	  in_compressediEt = nBinsIEt-1;
	  EG_Et = EG_Et > ET_MAX ? ET_MAX : EG_Et;
	}
	
	std::map<short, short>::iterator NTTPos = lutMapNTT.find(int(EG_NTT));
	if (NTTPos != lutMapNTT.end()) in_compressedNTT = NTTPos->second;
	else in_compressedNTT = nBinsNTT-1;
        
	IsoCut_Progression = ResultProgressionName->GetBinContent(in_compressediEta+1,in_compressediEt+1,in_compressedNTT+1);
	
	if(EG_Iso_Et<=IsoCut_Progression and EG_Et > maxOptEt ) {
       maxOptEt=EG_Et;
      }


    eT1ForDoubleEG.Fill( EG_Et );
    if( EG_Iso_Et<=IsoCut_Progression  )
        eTLooseForDoubleEG.Fill( EG_Et );
    }
    

     if(maxOptEt < 1e8)
     {
	        ptMap_[pt_Progression_]->Fill(maxOptEt,puWeight);
     }
    
    maxDoubleEGIsoEt=-1e3;
    for(UInt_t i=10;i< (ET_MAX ) ;i++) 
      {
             if( eT1ForDoubleEG.Integral(i+1 - 10, ET_MAX ) >= 2 ) 
             {
                if(eTLooseForDoubleEG.Integral(i + 1 , ET_MAX ) >=1)
                {
                      maxDoubleEGIsoEt = i ;
                }
                else
                {
                   break;
                }
             }
             else
             {
                  break;
             }
      }
    
     if( maxDoubleEGIsoEt >= 0.0)
      {
        //std::cout<<"\t\t\t  Fiering EGDouble : "<<maxDoubleEGIsoEt<<"\n";
        ptMap_[pt_DoubleEGProgression_]->Fill(maxDoubleEGIsoEt,puWeight);
      }
    }

    den++;
    }
  
  std::cout<<"Den , nEnries, nEvents : "<<den<<","<<nEvents<<std::endl;
  double scale = (11.2456 * bunches)/nEvents;

  //Filling Rate Histos  
  for(UInt_t i=0;i< ET_MAX;i++) {

      rateMap_["rate_Run2"]      ->SetBinContent(i+1, ptMap_["pt_Run2"]->Integral(i+1,ET_MAX)*scale);
      rateMap_["rate_Run2_loose"]->SetBinContent(i+1, ptMap_["pt_Run2_loose"]->Integral(i+1,ET_MAX)*scale);
      rateMap_["rate_Run2_tight"]->SetBinContent(i+1, ptMap_["pt_Run2_tight"]->Integral(i+1,ET_MAX)*scale);
      rateMap_["rate_Run2_DoubleEGLoose_diff10"]->SetBinContent(i+1, ptMap_["pt_Run2_DoubleEGLoose_diff10"]->Integral(i+1,ET_MAX)*scale);

      //std::cout<<" Et = "<<i<<" : Integral : "<<ptMap_["pt_Run2"]->Integral(i+1,ET_MAX)<<" rate : "<<ptMap_["pt_Run2"]->Integral(i+1,ET_MAX)*scale<<"\n";
      //std::cout<<"    loose "<<i<<" : Integral : "<<ptMap_["pt_Run2_loose"]->Integral(i+1,ET_MAX)<<" rate : "<<ptMap_["pt_Run2_loose"]->Integral(i+1,ET_MAX)*scale<<"\n";
      //std::cout<<"    tight "<<i<<" : Integral : "<<ptMap_["pt_Run2_tight"]->Integral(i+1,ET_MAX)<<" rate : "<<ptMap_["pt_Run2_tight"]->Integral(i+1,ET_MAX)*scale<<"\n";
    
    for (auto it :lutProgOptVec_) {
      TString CurrentNameHisto = "pt_Progression";
      TString CurrentNameHisto2 = "pt_DoubleEGProgression";
      TString CurrentNameHisto1= "rate_Progression";
      TString CurrentNameHisto3= "rate_DoubleEGProgression";
      CurrentNameHisto += it;
      CurrentNameHisto2 += it;
      CurrentNameHisto1 += it;
      CurrentNameHisto3 += it;
      rateMap_[CurrentNameHisto1]->SetBinContent(i+1, ptMap_[CurrentNameHisto]->Integral(i+1,ET_MAX)*scale);
      rateMap_[CurrentNameHisto3]->SetBinContent(i+1, ptMap_[CurrentNameHisto2]->Integral(i+1,ET_MAX)*scale);
    }

  }

  
  //Finding Et for fixed rate (frate)
  int xbin;
  for (auto it :lutProgOptVec_) {

    TString CurrentNameHisto1= "rate_Progression";
    CurrentNameHisto1 += it;
    float vl,bl,rate; //last value, last bin
    for(UInt_t i=0;i< ET_MAX;i++) {
      rate = rateMap_[CurrentNameHisto1]->GetBinContent(i+1);
      if(rate > frate) {
        bl =i;
	vl=rate;
      }
      
      if(rate <=frate){
	   xbin = (rate - frate) < abs(rate-vl) ? i : bl;
	break;    
      }
    }

    std:: cout<<"option: "<<it<<" , "<<"et: "<<xbin<<" +/- " <<vl/sqrt(vl/scale)<<std::endl;
    et_option.push_back(xbin);
 
    CurrentNameHisto1= "rate_DoubleEGProgression";
    CurrentNameHisto1 += it;
    for(UInt_t i=0;i< ET_MAX;i++) {
      rate = rateMap_[CurrentNameHisto1]->GetBinContent(i+1);
      if(rate > frateDoubleEG) {
        bl =i;
	vl=rate;
      }
      
      if(rate <=frateDoubleEG){
	xbin = (rate - frateDoubleEG) < abs(rate-vl) ? i : bl;
	break;    
      }
    }     
    std:: cout<<"option: "<<it<<" , Double EG "<<" et: "<<xbin<<"_"<<xbin-10<<" +/- " <<vl/sqrt(vl/scale)<<std::endl;
    etDoubleEG_option.push_back(xbin);
  }


  //booking histograms for Single EG fixed Et for all options
  int n=-1;
  for (auto it :lutProgOptVec_) {    
    n++;
    TString PtPassName_fr_ = "pt_pass_option";
    TString turnOn_Option_fr_ = "turnOn_Option";
    
    PtPassName_fr_ += it;
    PtPassName_fr_ +="_Et_";
    PtPassName_fr_=PtPassName_fr_ + std::to_string(et_option[n])  + "_fr";
    td2->cd();
    //TH1F*PtPass_fr= new TH1F(PtPassName_fr_, PtPassName_fr_, 21,binning);
    TH1F*PtPass_fr= new TH1F(PtPassName_fr_, PtPassName_fr_, nBins_fine,xEdges_fine);
    pt_pass_Map_.insert(std::make_pair(PtPassName_fr_, PtPass_fr));


    turnOn_Option_fr_ += it;
    turnOn_Option_fr_ +="_Et_";
    turnOn_Option_fr_ = turnOn_Option_fr_ + std::to_string(et_option[n]) + "_fr";
    
    if(!check_turn_on_dir) {                                                                                                     
      td3 = outputFile_->mkdir("turn_on_progression");                                                                           
      check_turn_on_dir = true;                                                                                                  
    }                                                                                                                            
    td3->cd();
    TGraphAsymmErrors* turnOn_Option_fr;
    turnOn_Map_.insert(std::make_pair(turnOn_Option_fr_,turnOn_Option_fr));
  }


  //booking histograms for Double EG fixed Rate for all options
  n=-1;
  for (auto it :lutProgOptVec_) {    
    n++;
    TString PtPassName_fr_ = "pt_DoubleEGIsoPass_option";
    TString turnOn_Option_fr_ = "turnOnDoubleEGIso_Option";
    
    PtPassName_fr_ += it;
    PtPassName_fr_ +="_Et_";
    PtPassName_fr_=PtPassName_fr_ + std::to_string(etDoubleEG_option[n])  + "_fr";
    td2->cd();
    //TH1F*PtPass_fr= new TH1F(PtPassName_fr_, PtPassName_fr_, 21,binning);
    TH1F*PtPass_fr= new TH1F(PtPassName_fr_, PtPassName_fr_, nBins_fine,xEdges_fine);
    pt_pass_Map_.insert(std::make_pair(PtPassName_fr_, PtPass_fr));


    turnOn_Option_fr_ += it;
    turnOn_Option_fr_ +="_Et_";
    turnOn_Option_fr_ = turnOn_Option_fr_ + std::to_string(etDoubleEG_option[n]) + "_fr";
    
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
  Long64_t nbytes1 = 0, nb1 = 0;
  double sum=0;
  std::cout << "Total available  Entries For Efficiency " << nEntries1_ << std::endl;
  if(maxEntriesForEfficiency >0 ) nEntries1_= maxEntriesForEfficiency < nEntries1_ ? maxEntriesForEfficiency : nEntries1_;
  std::cout<<"Processing a total of "<<nEntries1_<<" Entries \n";

  t_start = std::chrono::high_resolution_clock::now();
  for (Long64_t jentry=0; jentry< nEntries1_; jentry++) {
    Long64_t ientry = fChain1->LoadTree(jentry);
    
    if(jentry%reportEvery == 0 )
       {
             t_end = std::chrono::high_resolution_clock::now();
             std::cout<<"Processing Entry in event loop : "<<jentry<<" / "<<nEntries1_<<"  [ "<<100.0*jentry/nEntries1_<<"  % ]  "
                      << " Elapsed time : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()/1000.0
                      <<"  Estimated time left : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()*( nEntries1_ - jentry)/(1e-9 + jentry)* 0.001
                      <<std::endl;
       }
 

    if (ientry < 0) break;
    nb1 = fChain1->GetEntry(jentry);   nbytes1 += nb1;
 //   if(l1tEmuNTT<60 && nTTRange_) continue;  //define range
    if(l1tEmuRawEt < 0.) continue;
    
    if( isProbeLoose==0 ) continue;
    if( fabs(eleProbeEta) >= 2.5) continue;
    if( sqrt(pow(eleProbeEta-eleTagEta,2)+pow(eleProbePhi-eleTagPhi,2)) < 0.6 ) continue;

    pT_all->Fill(eleProbeSclEt);
    sum++;
    
    std::map<short, short>::iterator EtaPos = lutMapEta.find(abs(l1tEmuTowerIEta));
    if (EtaPos != lutMapEta.end()) in_compressediEta = EtaPos->second;
    else in_compressediEta = nBinsIEta-1;
    
    if(l1tEmuPt < 0.) continue;
    //    if(l1tEmuNTT < 0) continue; 
    std::map<short, short>::iterator EtPos = lutMapEt.find(l1tEmuRawEt);
    if (EtPos != lutMapEt.end()) in_compressediEt = EtPos->second;
    else in_compressediEt = nBinsIEt-1;

    std::map<short, short>::iterator NTTPos = lutMapNTT.find(l1tEmuNTT);
    if (NTTPos != lutMapNTT.end()) in_compressedNTT = NTTPos->second;
    else in_compressedNTT = nBinsNTT-1;
    Int_t IsoCut_Progression;
    for(UInt_t e = 15 ; e <= 50 ; e += 5) {
      for (auto it :lutProgOptVec_) {
	  
      TString ResultProgressionName_= "LUT_Progression_";
      TString PtPassName_= "pt_pass_option";
	  ResultProgressionName_ +=it;
	  PtPassName_ += it;
	  PtPassName_ +="_Et_";
	  PtPassName_ += std::to_string(e);
	  TH3F* ResultProgressionName=optionsMap[ResultProgressionName_];
      IsoCut_Progression = ResultProgressionName->GetBinContent(in_compressediEta+1,in_compressediEt+1,in_compressedNTT+1);
	  
	  if(l1tEmuPt >= e && l1tEmuIsoEt <= IsoCut_Progression){
	    pt_pass_Map_[PtPassName_]->Fill(eleProbeSclEt);
	  }

      }
    }
    int n=-1;
    // Fixed rate TurnOns
    for (auto it :lutProgOptVec_) {
      n++;
      TString PtPassName_fr= "pt_pass_option";
      TString ResultProgressionName_= "LUT_Progression_";      
      
      PtPassName_fr += it;
      PtPassName_fr = PtPassName_fr + "_Et_" + std::to_string(et_option[n]) +"_fr";
      ResultProgressionName_ +=it;
      
      //TH3F* ResultProgressionName = (TH3F*)gDirectory->Get(ResultProgressionName_.Data());
	  TH3F* ResultProgressionName=optionsMap[ResultProgressionName_];
      IsoCut_Progression = ResultProgressionName->GetBinContent(in_compressediEta+1,in_compressediEt+1,in_compressedNTT+1);
      
      if(l1tEmuPt >= et_option[n] && l1tEmuIsoEt <= IsoCut_Progression){
	      pt_pass_Map_[PtPassName_fr]->Fill(eleProbeSclEt);
      }

      PtPassName_fr= "pt_DoubleEGIsoPass_option";
      PtPassName_fr += it;
      PtPassName_fr = PtPassName_fr + "_Et_" + std::to_string(etDoubleEG_option[n] ) +"_fr";
      
      if( ( l1tEmuPt >= etDoubleEG_option[n] ) && (l1tEmuIsoEt <= IsoCut_Progression)){
	      pt_pass_Map_[PtPassName_fr]->Fill(eleProbeSclEt);
      }

    }

    // using Run2 LUT                                                                                                              
    if(hasL1Emu_32==1) pt_pass_Map_Run2_["pt_pass_Et_32"]->Fill(eleProbeSclEt);
    if(hasL1Emu_31==1) pt_pass_Map_Run2_["pt_pass_Et_31"]->Fill(eleProbeSclEt);
    if(hasL1Emu_30==1) pt_pass_Map_Run2_["pt_pass_Et_30"]->Fill(eleProbeSclEt);
    if(hasL1Emu_29==1) pt_pass_Map_Run2_["pt_pass_Et_29"]->Fill(eleProbeSclEt);
    if(hasL1Emu_28==1) pt_pass_Map_Run2_["pt_pass_Et_28"]->Fill(eleProbeSclEt);
    if(hasL1Emu_27==1) pt_pass_Map_Run2_["pt_pass_Et_27"]->Fill(eleProbeSclEt);
    if(hasL1Emu_26==1) pt_pass_Map_Run2_["pt_pass_Et_26"]->Fill(eleProbeSclEt);
    if(hasL1Emu_25==1) pt_pass_Map_Run2_["pt_pass_Et_25"]->Fill(eleProbeSclEt);
    if(hasL1Emu_24==1) pt_pass_Map_Run2_["pt_pass_Et_24"]->Fill(eleProbeSclEt);
    if(hasL1Emu_23==1) pt_pass_Map_Run2_["pt_pass_Et_23"]->Fill(eleProbeSclEt);
    if(hasL1Emu_22==1) pt_pass_Map_Run2_["pt_pass_Et_22"]->Fill(eleProbeSclEt);
    if(hasL1Emu_21==1) pt_pass_Map_Run2_["pt_pass_Et_21"]->Fill(eleProbeSclEt);
    if(hasL1Emu_20==1) pt_pass_Map_Run2_["pt_pass_Et_20"]->Fill(eleProbeSclEt);
    
    if(hasL1Emu_looseiso32==1) pt_pass_Map_Run2_["pt_pass_Et_32_loose"]->Fill(eleProbeSclEt);
    if(hasL1Emu_looseiso31==1) pt_pass_Map_Run2_["pt_pass_Et_31_loose"]->Fill(eleProbeSclEt);
    if(hasL1Emu_looseiso30==1) pt_pass_Map_Run2_["pt_pass_Et_30_loose"]->Fill(eleProbeSclEt);
    if(hasL1Emu_looseiso29==1) pt_pass_Map_Run2_["pt_pass_Et_29_loose"]->Fill(eleProbeSclEt);
    if(hasL1Emu_looseiso28==1) pt_pass_Map_Run2_["pt_pass_Et_28_loose"]->Fill(eleProbeSclEt);
    if(hasL1Emu_looseiso27==1) pt_pass_Map_Run2_["pt_pass_Et_27_loose"]->Fill(eleProbeSclEt);
    if(hasL1Emu_looseiso26==1) pt_pass_Map_Run2_["pt_pass_Et_26_loose"]->Fill(eleProbeSclEt);
    if(hasL1Emu_looseiso25==1) pt_pass_Map_Run2_["pt_pass_Et_25_loose"]->Fill(eleProbeSclEt);
    if(hasL1Emu_looseiso24==1) pt_pass_Map_Run2_["pt_pass_Et_24_loose"]->Fill(eleProbeSclEt);
    if(hasL1Emu_looseiso23==1) pt_pass_Map_Run2_["pt_pass_Et_23_loose"]->Fill(eleProbeSclEt);
    if(hasL1Emu_looseiso22==1) pt_pass_Map_Run2_["pt_pass_Et_22_loose"]->Fill(eleProbeSclEt);
    if(hasL1Emu_looseiso21==1) pt_pass_Map_Run2_["pt_pass_Et_21_loose"]->Fill(eleProbeSclEt);
    if(hasL1Emu_looseiso20==1) pt_pass_Map_Run2_["pt_pass_Et_20_loose"]->Fill(eleProbeSclEt);
    
    if(hasL1Emu_tightiso32==1) pt_pass_Map_Run2_["pt_pass_Et_32_tight"]->Fill(eleProbeSclEt);
    if(hasL1Emu_tightiso31==1) pt_pass_Map_Run2_["pt_pass_Et_31_tight"]->Fill(eleProbeSclEt);
    if(hasL1Emu_tightiso30==1) pt_pass_Map_Run2_["pt_pass_Et_30_tight"]->Fill(eleProbeSclEt);
    if(hasL1Emu_tightiso29==1) pt_pass_Map_Run2_["pt_pass_Et_29_tight"]->Fill(eleProbeSclEt);
    if(hasL1Emu_tightiso28==1) pt_pass_Map_Run2_["pt_pass_Et_28_tight"]->Fill(eleProbeSclEt);
    if(hasL1Emu_tightiso27==1) pt_pass_Map_Run2_["pt_pass_Et_27_tight"]->Fill(eleProbeSclEt);
    if(hasL1Emu_tightiso26==1) pt_pass_Map_Run2_["pt_pass_Et_26_tight"]->Fill(eleProbeSclEt);
    if(hasL1Emu_tightiso25==1) pt_pass_Map_Run2_["pt_pass_Et_25_tight"]->Fill(eleProbeSclEt);
    if(hasL1Emu_tightiso24==1) pt_pass_Map_Run2_["pt_pass_Et_24_tight"]->Fill(eleProbeSclEt);
    if(hasL1Emu_tightiso23==1) pt_pass_Map_Run2_["pt_pass_Et_23_tight"]->Fill(eleProbeSclEt);
    if(hasL1Emu_tightiso22==1) pt_pass_Map_Run2_["pt_pass_Et_22_tight"]->Fill(eleProbeSclEt);
    if(hasL1Emu_tightiso21==1) pt_pass_Map_Run2_["pt_pass_Et_21_tight"]->Fill(eleProbeSclEt);
    if(hasL1Emu_tightiso20==1) pt_pass_Map_Run2_["pt_pass_Et_20_tight"]->Fill(eleProbeSclEt);
  }
  t_end = std::chrono::high_resolution_clock::now();
  std::cout << " Turn Ons Done: Elapsed time : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()/1000.0<<std::endl;
  //Turnons for various Ets
  for(UInt_t e = 15 ; e <= 50 ; e += 5) {
    for (auto it :lutProgOptVec_) {
      outputFile_->cd();
      
      if(!check_turn_on_dir) {
	    td3 = outputFile_->mkdir("turn_on_progression");
        check_turn_on_dir = true;
      }
      td3->cd();

      TString PtPassName_= "pt_pass_option";
      TString turnOn_Option_="turnOn_Option";
      
      PtPassName_ += it;
      PtPassName_ +="_Et_";
      PtPassName_ += std::to_string(e);
      turnOn_Option_ += it;
      turnOn_Option_ +="_Et_";
      turnOn_Option_ += std::to_string(e);
      turnOn_Map_[turnOn_Option_]=new TGraphAsymmErrors( pt_pass_Map_[PtPassName_],pT_all,"cp"); 
      turnOn_Map_[turnOn_Option_]->Write();
      
      double acceptance = pt_pass_Map_[PtPassName_]->GetEntries()/nEntries1_ ;
      th1fStore["FixedEtTurnons"]->Fill(turnOn_Map_[turnOn_Option_]->GetName(),e);
      th1fStore["FixedEtTurnons_Acceptance"]->Fill(turnOn_Map_[turnOn_Option_]->GetName(),acceptance);
      
      PtPassName_= "pt_pass_DoubleEGOption";
      turnOn_Option_="turnOn_DoubleEGOption";
      PtPassName_ += it;
      PtPassName_ +="_Et_";
      PtPassName_ += std::to_string(e)+"_"+std::to_string(e-10);
      turnOn_Option_ += it;
      turnOn_Option_ +="_Et_";
      turnOn_Option_ += std::to_string(e)+"_"+std::to_string(e-10);
      turnOn_Map_[turnOn_Option_]=new TGraphAsymmErrors( pt_pass_Map_[PtPassName_],pT_all,"cp"); 
      turnOn_Map_[turnOn_Option_]->Write();
      acceptance = pt_pass_Map_[PtPassName_]->GetEntries()/nEntries1_ ;

      th1fStore["FixedEtDoubleEGTurnons"]->Fill(turnOn_Map_[turnOn_Option_]->GetName(),e);
      th1fStore["FixedEtDoubleEGTurnons_Acceptance"]->Fill(turnOn_Map_[turnOn_Option_]->GetName(),acceptance);
    }
  }

  //Turnons for fixed Rate 
  n=-1;
  for (auto it :lutProgOptVec_) {
    n++;
    TString PtPassName_fr= "pt_pass_option";
    TString turnOn_Option_fr="turnOn_Option";
    PtPassName_fr += it;
    PtPassName_fr = PtPassName_fr + "_Et_" + std::to_string(et_option[n]) +"_fr";
    td3->cd();
    turnOn_Option_fr += it;
    turnOn_Option_fr = turnOn_Option_fr + "_Et_" + std::to_string(et_option[n]) +"_fr";;
    turnOn_Map_[turnOn_Option_fr]=new TGraphAsymmErrors( pt_pass_Map_[PtPassName_fr],pT_all,"cp");
    turnOn_Map_[turnOn_Option_fr]->Write();

    double acceptance = pt_pass_Map_[PtPassName_fr]->GetEntries() ;
    th1fStore["FixedRateTurnons"]->Fill(turnOn_Map_[turnOn_Option_fr]->GetName(),et_option[n]);
    th1fStore["FixedRateTurnons_Acceptance"]->Fill(turnOn_Map_[turnOn_Option_fr]->GetName(),acceptance);

    PtPassName_fr= "pt_DoubleEGIsoPass_option";
    turnOn_Option_fr="turnOnDoubleEGIso_Option";
    PtPassName_fr += it;
    PtPassName_fr = PtPassName_fr + "_Et_" + std::to_string(etDoubleEG_option[n]) +"_fr";
    td3->cd();
    turnOn_Option_fr += it;
    turnOn_Option_fr = turnOn_Option_fr + "_Et_" + std::to_string(etDoubleEG_option[n]) +"_fr";;
    turnOn_Map_[turnOn_Option_fr]=new TGraphAsymmErrors( pt_pass_Map_[PtPassName_fr],pT_all,"cp");
    turnOn_Map_[turnOn_Option_fr]->Write();
    acceptance = pt_pass_Map_[PtPassName_fr]->GetEntries() ;
    th1fStore["FixedRateDoubleEGTurnons"]->Fill(turnOn_Map_[turnOn_Option_fr]->GetName(),etDoubleEG_option[n]);
    th1fStore["FixedRateDoubleEGTurnons_Acceptance"]->Fill(turnOn_Map_[turnOn_Option_fr]->GetName(),acceptance);

  }

  // using Run2 LUT                                                                                                                                      
  td4->cd();
  //std::vector<int> v = {24, 26};
  std::vector<int> v = {20,21,22,23,24,25,26,27,28,29,30,31,32};
  for(int e : v) {
    TString PtPassName_= "pt_pass";
    TString turnOnName_="turnOn";
    PtPassName_ +="_Et_";
    PtPassName_ += std::to_string(e);
    TString PtPassLooseIsoName_=PtPassName_+"_loose";
    TString PtPassTightIsoName_=PtPassName_+"_tight";

    turnOnName_ +="_Et_";
    turnOnName_ += std::to_string(e);

    TString turnOnLooseIsoName_=turnOnName_+"_loose";
    TString turnOnTightIsoName_=turnOnName_+"_tight";
    turnOn_Map_Run2_[turnOnName_]=new TGraphAsymmErrors( pt_pass_Map_Run2_[PtPassName_],pT_all,"cp");
    turnOn_Map_Run2_[turnOnLooseIsoName_]=new TGraphAsymmErrors( pt_pass_Map_Run2_[PtPassLooseIsoName_],pT_all,"cp");
    turnOn_Map_Run2_[turnOnTightIsoName_]=new TGraphAsymmErrors( pt_pass_Map_Run2_[PtPassTightIsoName_],pT_all,"cp");

    turnOn_Map_Run2_[turnOnName_]->Write();
    turnOn_Map_Run2_[turnOnLooseIsoName_]->Write();
    turnOn_Map_Run2_[turnOnTightIsoName_]->Write();

    double acceptance = pt_pass_Map_Run2_[PtPassName_]->GetEntries() ;
    double acceptance_loose = pt_pass_Map_Run2_[PtPassLooseIsoName_]->GetEntries() ;
    double acceptance_tight = pt_pass_Map_Run2_[PtPassTightIsoName_]->GetEntries() ;

    th1fStore["Run2Turnons_Acceptance"]->Fill(turnOn_Map_Run2_[turnOnName_]->GetName(),acceptance);
    th1fStore["Run2Turnons_LooseIso_Acceptance"]->Fill(turnOn_Map_Run2_[turnOnLooseIsoName_]->GetName(),acceptance_loose);
    th1fStore["Run2Turnons_TightIso_Acceptance"]->Fill(turnOn_Map_Run2_[turnOnTightIsoName_]->GetName(),acceptance_tight);
  }
  
}

void ApplyIsolation::bookHistograms() {
  outputFile_ = new TFile(outputFileName_.c_str(), "RECREATE");
  outputFile_->cd();
  
  TString CurrentNameHisto = "pt_Run2";
  TString CurrentNameHisto1= "rate_Run2";
  TString CurrentNameHisto2= "";
  TString CurrentNameHisto3= "";
  TH1F* pt_Progression;
  TH1F* rate_Progression;
  if(!check_pt_rate_dir) { 
      td = outputFile_->mkdir("pt_progression_for_rate");
      check_pt_rate_dir = true;
    }
  td->cd();
  CurrentNameHisto  = "pt_Run2";
  pt_Progression=  new TH1F(CurrentNameHisto, CurrentNameHisto ,ET_MAX ,  0.0-0.4 , ET_MAX-0.4);
  ptMap_.insert(std::make_pair(CurrentNameHisto,pt_Progression));
  CurrentNameHisto  = "pt_Run2_loose";
  pt_Progression=  new TH1F(CurrentNameHisto, CurrentNameHisto ,ET_MAX ,  0.0-0.4 , ET_MAX-0.4);
  ptMap_.insert(std::make_pair(CurrentNameHisto,pt_Progression));
  CurrentNameHisto  = "pt_Run2_tight";
  pt_Progression=  new TH1F(CurrentNameHisto, CurrentNameHisto ,ET_MAX ,  0.0-0.4 , ET_MAX-0.4);
  ptMap_.insert(std::make_pair(CurrentNameHisto,pt_Progression));
  CurrentNameHisto  = "pt_Run2_DoubleEGLoose_diff10";
  pt_Progression=  new TH1F(CurrentNameHisto, CurrentNameHisto ,ET_MAX ,  0.0-0.4 , ET_MAX-0.4);
  ptMap_.insert(std::make_pair(CurrentNameHisto,pt_Progression));
  

  if(!check_rate_dir) {
      td1 = outputFile_->mkdir("rate_progression");
      check_rate_dir = true;
    }
   td1->cd();   
   CurrentNameHisto1 = "rate_Run2";
   rate_Progression = new TH1F(CurrentNameHisto1, CurrentNameHisto1, ET_MAX ,  0.0-0.4 , ET_MAX-0.4 );
   rateMap_.insert(std::make_pair(CurrentNameHisto1,rate_Progression));
   CurrentNameHisto1 = "rate_Run2_loose";
   rate_Progression = new TH1F(CurrentNameHisto1, CurrentNameHisto1, ET_MAX ,  0.0-0.4 , ET_MAX-0.4 );
   rateMap_.insert(std::make_pair(CurrentNameHisto1,rate_Progression));
   CurrentNameHisto1 = "rate_Run2_tight";
   rate_Progression = new TH1F(CurrentNameHisto1, CurrentNameHisto1, ET_MAX ,  0.0-0.4 , ET_MAX-0.4 );
   rateMap_.insert(std::make_pair(CurrentNameHisto1,rate_Progression));
   CurrentNameHisto1 = "rate_Run2_DoubleEGLoose_diff10";
   rate_Progression = new TH1F(CurrentNameHisto1, CurrentNameHisto1, ET_MAX ,  0.0-0.4 , ET_MAX-0.4 );
   rateMap_.insert(std::make_pair(CurrentNameHisto1,rate_Progression));

  for (auto it : lutProgOptVec_) {    
    CurrentNameHisto = "pt_Progression";
    CurrentNameHisto1= "rate_Progression";
    CurrentNameHisto2= "pt_DoubleEGProgression";
    CurrentNameHisto3= "rate_DoubleEGProgression";
    CurrentNameHisto += it;
    CurrentNameHisto1 += it;
    CurrentNameHisto2 += it;
    CurrentNameHisto3 += it;
    
    if(!check_pt_rate_dir) { 
      td = outputFile_->mkdir("pt_progression_for_rate");
      check_pt_rate_dir = true;
    }

    td->cd();
    pt_Progression=  new TH1F(CurrentNameHisto, CurrentNameHisto ,ET_MAX ,  0.0-0.4 , ET_MAX-0.4);
    ptMap_.insert(std::make_pair(CurrentNameHisto,pt_Progression));
    pt_Progression=  new TH1F(CurrentNameHisto2, CurrentNameHisto2 ,ET_MAX ,  0.0-0.4 , ET_MAX-0.4);
    ptMap_.insert(std::make_pair(CurrentNameHisto2,pt_Progression));
    if(!check_rate_dir) {
      td1 = outputFile_->mkdir("rate_progression");
      check_rate_dir = true;
    }
    td1->cd();   

    rate_Progression = new TH1F(CurrentNameHisto1, CurrentNameHisto1, ET_MAX ,  0.0-0.4 , ET_MAX-0.4 );
    rateMap_.insert(std::make_pair(CurrentNameHisto1,rate_Progression));
    
    rate_Progression = new TH1F(CurrentNameHisto3, CurrentNameHisto3, ET_MAX ,  0.0-0.4 , ET_MAX-0.4 );
    rateMap_.insert(std::make_pair(CurrentNameHisto3,rate_Progression));
    
    for(UInt_t e = 15 ; e <= 50 ; e += 5) {

      
      if(!check_pt_turn_on_dir) {
	        td2 = outputFile_->mkdir("pt_progression_for_turnon");
	        check_pt_turn_on_dir = true;
      }

      td2->cd();
      
      // Single EG Turnons

      TString PtPassName_    = "pt_pass_option";
      TString turnOn_Option_ = "turnOn_Option";
      PtPassName_ += it;
      PtPassName_ +="_Et_";
      PtPassName_ += std::to_string(e);  
      turnOn_Option_ += it;
      turnOn_Option_ +="_Et_";
      turnOn_Option_ += std::to_string(e);

      TH1F*PtPass= new TH1F(PtPassName_, PtPassName_, nBins_fine,xEdges_fine);
      pt_pass_Map_.insert(std::make_pair(PtPassName_, PtPass));
      TGraphAsymmErrors* turnOn_Option;
      turnOn_Map_.insert(std::make_pair(turnOn_Option_,turnOn_Option));
      
      PtPassName_= "pt_pass_DoubleEGOption";
      turnOn_Option_="turnOn_DoubleEGOption";
      
      PtPassName_ += it;
      PtPassName_ +="_Et_";
      PtPassName_ += std::to_string(e)+"_"+std::to_string(e-10);  
      turnOn_Option_ += it;
      turnOn_Option_ +="_Et_";
      turnOn_Option_ += std::to_string(e)+"_"+std::to_string(e-10);;
      
      PtPass= new TH1F(PtPassName_, PtPassName_, nBins_fine,xEdges_fine);
      pt_pass_Map_.insert(std::make_pair(PtPassName_, PtPass));
      TGraphAsymmErrors* turnOn_Option2;
      turnOn_Map_.insert(std::make_pair(turnOn_Option_,turnOn_Option2));
      
    }
  }
  td2->cd();
  //pT_all = new TH1F("pT_all","pT_all",21,binning);
  pT_all = new TH1F("pT_all","pT_all",nBins_fine,xEdges_fine);
  
  th1fStore["FixedRateTurnons"] = new TH1F("FixedRateTurnons","",3,0.0,3.0);
  th1fStore["FixedRateTurnons"]->SetCanExtend(TH1::kAllAxes);
  th1fStore["FixedRateTurnons_Acceptance"] = new TH1F("FixedRateTurnons_Acceptance","",3,0.0,3.0);
  th1fStore["FixedRateTurnons_Acceptance"]->SetCanExtend(TH1::kAllAxes);  
  
  th1fStore["FixedRateDoubleEGTurnons"] = new TH1F("FixedRateDoubleEGTurnons","",3,0.0,3.0);
  th1fStore["FixedRateDoubleEGTurnons"]->SetCanExtend(TH1::kAllAxes);
  th1fStore["FixedRateDoubleEGTurnons_Acceptance"] = new TH1F("FixedRateDoubleEGTurnons_Acceptance","",3,0.0,3.0);
  th1fStore["FixedRateDoubleEGTurnons_Acceptance"]->SetCanExtend(TH1::kAllAxes);  
  
  th1fStore["FixedEtTurnons"]   = new TH1F("FixedEtTurnons","",3,0.0,3.0);
  th1fStore["FixedEtTurnons"]->SetCanExtend(TH1::kAllAxes);
  th1fStore["FixedEtTurnons_Acceptance"]   = new TH1F("FixedEtTurnons_Acceptance","",3,0.0,3.0);
  th1fStore["FixedEtTurnons_Acceptance"]->SetCanExtend(TH1::kAllAxes);

  th1fStore["FixedEtDoubleEGTurnons"]   = new TH1F("FixedEtDoubleEGTurnons","",3,0.0,3.0);
  th1fStore["FixedEtDoubleEGTurnons"]->SetCanExtend(TH1::kAllAxes);
  th1fStore["FixedEtDoubleEGTurnons_Acceptance"]   = new TH1F("FixedEtDoubleEGTurnons_Acceptance","",3,0.0,3.0);
  th1fStore["FixedEtDoubleEGTurnons_Acceptance"]->SetCanExtend(TH1::kAllAxes);

  //using Run2 LUT
  std::vector<int> v = {20,21,22,23,24,25,26,27,28,29,30,31,32};
  for(int e : v) {
    TString PtPassName_= "pt_pass";
    TString turnOnName_="turnOn";
    PtPassName_ +="_Et_";
    PtPassName_ += std::to_string(e);
    TString PtPassLooseIsoName_=PtPassName_+"_loose";
    TString PtPassTightIsoName_=PtPassName_+"_tight";
    
    if(!check_turn_on_dir_Run2) {
      td4 = outputFile_->mkdir("turnon_progression_Run2");
      check_turn_on_dir_Run2 = true;
    }
    td4->cd();
    //TH1F*PtPass= new TH1F(PtPassName_, PtPassName_, 38,binning);
    //TH1F*PtPassLoose= new TH1F(PtPassLooseIsoName_, PtPassLooseIsoName_, 38,binning);
    //TH1F*PtPassTight= new TH1F(PtPassTightIsoName_, PtPassTightIsoName_, 38,binning);
    TH1F*PtPass= new TH1F(PtPassName_, PtPassName_, nBins_fine,xEdges_fine);
    TH1F*PtPassLoose= new TH1F(PtPassLooseIsoName_, PtPassLooseIsoName_, nBins_fine,xEdges_fine);
    TH1F*PtPassTight= new TH1F(PtPassTightIsoName_, PtPassTightIsoName_, nBins_fine,xEdges_fine);
    
    pt_pass_Map_Run2_.insert(std::make_pair(PtPassName_, PtPass));
    pt_pass_Map_Run2_.insert(std::make_pair(PtPassLooseIsoName_, PtPassLoose));
    pt_pass_Map_Run2_.insert(std::make_pair(PtPassTightIsoName_, PtPassTight));
    
    turnOnName_="turnOn";
    turnOnName_ +="_Et_";
    turnOnName_ += std::to_string(e);

    TString turnOnLooseIsoName_=turnOnName_+"_loose";
    TString turnOnTightIsoName_=turnOnName_+"_tight";

    TGraphAsymmErrors* turnOn;
    TGraphAsymmErrors* turnOnLoose;
    TGraphAsymmErrors* turnOnTight;

    turnOn_Map_Run2_.insert(std::make_pair(turnOnName_,turnOn));
    turnOn_Map_Run2_.insert(std::make_pair(turnOnLooseIsoName_,turnOnLoose));
    turnOn_Map_Run2_.insert(std::make_pair(turnOnTightIsoName_,turnOnTight));
  }


  th1fStore["Run2Turnons_Acceptance"]   = new TH1F("Run2Turnons_Acceptance","",3,0.0,3.0);
  th1fStore["Run2Turnons_Acceptance"]->SetCanExtend(TH1::kAllAxes);
  th1fStore["Run2Turnons_LooseIso_Acceptance"]   = new TH1F("Run2Turnons_LooseIso_Acceptance","",3,0.0,3.0);
  th1fStore["Run2Turnons_LooseIso_Acceptance"]->SetCanExtend(TH1::kAllAxes);
  th1fStore["Run2Turnons_TightIso_Acceptance"]   = new TH1F("Run2Turnons_TightIso_Acceptance","",3,0.0,3.0);
  th1fStore["Run2Turnons_TightIso_Acceptance"]->SetCanExtend(TH1::kAllAxes);
  

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
    outputFile_->cd();
    for (std::map<TString,TH1F *>::iterator it=th1fStore.begin() ; it!=th1fStore.end(); ++it)
      {
        auto &ahist = *(it->second); 
        ahist.Write();
      }
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
      else if (key=="NTTLUTFileName")   readLUTTable(value,nBinsNTT, lutMapNTT);
      else if (key=="nTTRange") nTTRange_ = std::stoi(value.c_str());
      else if (key=="MaxEntriesForEfficency") maxEntriesForEfficiency = atoi(value.c_str());
      else if (key=="MaxEntriesForRate")   maxEntriesForRate = atoi(value.c_str());
      else if (key=="FixedRate") frate = std::stof(value.c_str());
      else if (key=="DoubleEGFixedRate") frateDoubleEG = std::stof(value.c_str());
      else if (key=="ReportEvery")    reportEvery = atoi(value.c_str());
      else if (key=="PUReweightingHistFile") {
         puReweightingFileName =  value.c_str();
         hasPUReweitingFile=true;
      }
      else if (key=="LUTProgressionOptions")
        {
	  std::string tmp_string = value;
	  std::vector<std::string> tmp_vec;
          tokenize(tmp_string,tmp_vec,",");
          for (auto it : tmp_vec) {
            lutProgOptVec_.push_back(it);
          }
        }
      else if (key=="ReportEvery")	    {
                reportEvery= atoi(value.c_str());
      }
      else if (key=="MaxEntriesForEfficency")	    {
                maxEntriesForEfficiency= atoi(value.c_str());
      }
      else if (key=="MaxEntriesForRate")	    {
                maxEntriesForRate= atoi(value.c_str());
      }
      else 
	std::cout << " unknown option " << " key " << key << std::endl;
    }
  }
  jobcardFile.close();
}


void ApplyIsolation::readTree() {
  fChain->SetMakeClass(1);
  //for rate
  fChain->SetBranchAddress("nEGs", &nEGs);
  fChain->SetBranchAddress("egEt", &egEt);
  fChain->SetBranchAddress("egTowerIEta", &egTowerIEta); 
  fChain->SetBranchAddress("egNTT", &egNTT);
  fChain->SetBranchAddress("egBx", &egBx);
  fChain->SetBranchAddress("egIsoEt", &egIsoEt);
  fChain->SetBranchAddress("egRawEt", &egRawEt);
  fChain->SetBranchAddress("egIso", &egIso);

  if(hasPUReweitingFile) fChain_evtTree->SetBranchAddress("nPV_True", &nVtx);

  //for turnon
  fChain1->SetBranchAddress("l1tEmuPt", &l1tEmuPt);
  fChain1->SetBranchAddress("l1tEmuNTT",&l1tEmuNTT);
  fChain1->SetBranchAddress("l1tEmuRawEt",&l1tEmuRawEt);
  fChain1->SetBranchAddress("l1tEmuTowerIEta",&l1tEmuTowerIEta);
  fChain1->SetBranchAddress("eleProbeSclEt",&eleProbeSclEt);
  fChain1->SetBranchAddress("l1tEmuIsoEt",&l1tEmuIsoEt);

  fChain1->SetBranchAddress("eleProbeEta"  ,&eleProbeEta			);
  fChain1->SetBranchAddress("eleProbePhi"  ,&eleProbePhi			);
  fChain1->SetBranchAddress("eleTagEta"    ,&eleTagEta		    	);
  fChain1->SetBranchAddress("eleTagPhi"    ,&eleTagPhi			    );
  fChain1->SetBranchAddress("isProbeLoose" ,&isProbeLoose			);
  //fChain1->SetBranchAddress("hasL1Emu_26",&hasL1Emu_26);
  //fChain1->SetBranchAddress("hasL1Emu_24",&hasL1Emu_24);
  //fChain1->SetBranchAddress("hasL1Emu_looseiso26",&hasL1Emu_looseiso26);
  //fChain1->SetBranchAddress("hasL1Emu_looseiso24",&hasL1Emu_looseiso24);
  //fChain1->SetBranchAddress("hasL1Emu_tightiso26",&hasL1Emu_tightiso26);
  //fChain1->SetBranchAddress("hasL1Emu_tightiso24",&hasL1Emu_tightiso24);

  fChain1->SetBranchAddress("hasL1Emu_20",&hasL1Emu_20);
  fChain1->SetBranchAddress("hasL1Emu_21",&hasL1Emu_21);
  fChain1->SetBranchAddress("hasL1Emu_22",&hasL1Emu_22);
  fChain1->SetBranchAddress("hasL1Emu_23",&hasL1Emu_23);
  fChain1->SetBranchAddress("hasL1Emu_24",&hasL1Emu_24);
  fChain1->SetBranchAddress("hasL1Emu_25",&hasL1Emu_25);
  fChain1->SetBranchAddress("hasL1Emu_26",&hasL1Emu_26);
  fChain1->SetBranchAddress("hasL1Emu_27",&hasL1Emu_27);
  fChain1->SetBranchAddress("hasL1Emu_28",&hasL1Emu_28);
  fChain1->SetBranchAddress("hasL1Emu_29",&hasL1Emu_29);
  fChain1->SetBranchAddress("hasL1Emu_30",&hasL1Emu_30);
  fChain1->SetBranchAddress("hasL1Emu_31",&hasL1Emu_31);
  fChain1->SetBranchAddress("hasL1Emu_32",&hasL1Emu_32);

  fChain1->SetBranchAddress("hasL1Emu_looseiso20",&hasL1Emu_looseiso20);
  fChain1->SetBranchAddress("hasL1Emu_looseiso21",&hasL1Emu_looseiso21);
  fChain1->SetBranchAddress("hasL1Emu_looseiso22",&hasL1Emu_looseiso22);
  fChain1->SetBranchAddress("hasL1Emu_looseiso23",&hasL1Emu_looseiso23);
  fChain1->SetBranchAddress("hasL1Emu_looseiso24",&hasL1Emu_looseiso24);
  fChain1->SetBranchAddress("hasL1Emu_looseiso25",&hasL1Emu_looseiso25);
  fChain1->SetBranchAddress("hasL1Emu_looseiso26",&hasL1Emu_looseiso26);
  fChain1->SetBranchAddress("hasL1Emu_looseiso27",&hasL1Emu_looseiso27);
  fChain1->SetBranchAddress("hasL1Emu_looseiso28",&hasL1Emu_looseiso28);
  fChain1->SetBranchAddress("hasL1Emu_looseiso29",&hasL1Emu_looseiso29);
  fChain1->SetBranchAddress("hasL1Emu_looseiso30",&hasL1Emu_looseiso30);
  fChain1->SetBranchAddress("hasL1Emu_looseiso31",&hasL1Emu_looseiso31);
  fChain1->SetBranchAddress("hasL1Emu_looseiso32",&hasL1Emu_looseiso32);

  fChain1->SetBranchAddress("hasL1Emu_tightiso20",&hasL1Emu_tightiso20);
  fChain1->SetBranchAddress("hasL1Emu_tightiso21",&hasL1Emu_tightiso21);
  fChain1->SetBranchAddress("hasL1Emu_tightiso22",&hasL1Emu_tightiso22);
  fChain1->SetBranchAddress("hasL1Emu_tightiso23",&hasL1Emu_tightiso23);
  fChain1->SetBranchAddress("hasL1Emu_tightiso24",&hasL1Emu_tightiso24);
  fChain1->SetBranchAddress("hasL1Emu_tightiso25",&hasL1Emu_tightiso25);
  fChain1->SetBranchAddress("hasL1Emu_tightiso26",&hasL1Emu_tightiso26);
  fChain1->SetBranchAddress("hasL1Emu_tightiso27",&hasL1Emu_tightiso27);
  fChain1->SetBranchAddress("hasL1Emu_tightiso28",&hasL1Emu_tightiso28);
  fChain1->SetBranchAddress("hasL1Emu_tightiso29",&hasL1Emu_tightiso29);
  fChain1->SetBranchAddress("hasL1Emu_tightiso30",&hasL1Emu_tightiso30);
  fChain1->SetBranchAddress("hasL1Emu_tightiso31",&hasL1Emu_tightiso31);
  fChain1->SetBranchAddress("hasL1Emu_tightiso32",&hasL1Emu_tightiso32);

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

 
