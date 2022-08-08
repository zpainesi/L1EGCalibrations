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
  } 
  else {

    //fChain = new TChain("L1UpgradeEmuTree");
    //fChain_1= new TChain("L1EventTree");
    fChain = new TChain("l1UpgradeEmuTree/L1UpgradeTree");
    fChain_1= new TChain("l1EventTree/L1EventTree");
    static constexpr int BUF_SIZE = 256;
    char buf[BUF_SIZE];
    while (myFileRate.getline(buf, BUF_SIZE, '\n')) {  // Pops off the newline character                                     

      std::string line(buf);
      fChain->AddFile(line.c_str(),-1);
      fChain_1->AddFile(line.c_str(),-1);
      std::cout << "Adding file " << line << " Entries " << fChain->GetEntries() <<  std::endl;
      //fChain->Print();    
      //fChain_1->Print();

    } 
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

void ApplyIsolation::createEfficiencyHistograms(TH1* h_P, TH1* h_F, TH1* h_E)
{
  Double_t num, den, ratio, eL, eU;
  double a  = 0.3173;
  double aeff = (1-a)/2;

  for (int k=1; k<=h_P->GetNbinsX(); ++k){
    num = h_P->GetBinContent(k);
    den = h_P->GetBinContent(k) + h_F->GetBinContent(k);
    eU = 0.0;
    eL = 0.0;
    h_E->SetBinErrorOption(TH1::kPoisson);
    if(den!=0 && num != 0){
      ratio = num/den;
      if (num > 100 || den*(1-ratio) > 100) {
	eU = sqrt(ratio*(1-ratio)/den);
	eL = sqrt(ratio*(1-ratio)/den);
      } else {
	eU = (1-BetaInverse(aeff,den-num,num+1))-ratio;
	eL = ratio-(1-BetaInverse(1-aeff,den-num+1,num));
      }

    } else {
      ratio = 0.0;
      eL = 0.0;
      eU = 0;
    }
    h_E->SetBinContent(k,ratio);
    h_E->SetBinError(k,eL);
    h_E->SetStats(0);
  }
}


void ApplyIsolation::loops() {
  short in_compressediEta;
  short in_compressediEt;
  short in_compressedNTT;
  Long64_t nbytes = 0, nbytes_1=0, nb = 0, nb_1=0;
  //  weightFile_ = new TFile(weightFileName_.c_str(), "READ");
  //TH1D*w= (TH1D*)weightFile_->Get("ratio_hist");
  optionsFile_ = new TFile(optionsFileName_.c_str(), "READ");
  optionsFile_->cd("Step2Histos");

  //Event loop for rate  ( SIngle EG + Double EG)                                                                                                  
  if (fChain == 0) return;
  nEntries_ = fChain->GetEntriesFast();
  std::cout << "Total available  Entries for Rate: " << nEntries_ << std::endl;
  if(maxEntriesForRate_ >0 ) nEntries_= maxEntriesForRate_ < nEntries_ ? maxEntriesForRate_ : nEntries_;
  std::cout<<"Processing a total of: "<<nEntries_<<" Entries \n";

  auto t_start = std::chrono::high_resolution_clock::now();
  auto t_end = std::chrono::high_resolution_clock::now();

  for (Long64_t jentry=0; jentry < nEntries_; jentry++) {
    Long64_t ientry = fChain->LoadTree(jentry);
    Long64_t i1entry = fChain_1->LoadTree(jentry);
    if (ientry < 0) break;
    if (i1entry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    nb_1 = fChain_1->GetEntry(jentry);   nbytes_1 += nb_1;
    
    if(jentry%reportEvery_ == 0 )
      {
	t_end = std::chrono::high_resolution_clock::now();
	std::cout<<"Processing Entry in event loop (Rate) : "<<jentry<<" / "<<nEntries_<<"  [ "<<100.0*jentry/nEntries_<<"  % ]  "
		 << " Elapsed time : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()/1000.0
		 <<"  Estimated time left : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()*( nEntries_ - jentry)/(1e-9 + jentry)* 0.001
		 <<std::endl;
      }
    //    weight=w->GetBinContent(nPV_True+1);
    //    std::cout<<"weight ="<<weight<<","<<nPV_True<<std::endl;
    //    sum_weight=sum_weight+weight; 
    if( nPV_True > 52) continue;
    if( nPV_True < 44 ) continue;

    if (nEGs < 1) continue;
	 
    for (auto it :lutProgOptVec_) {
      bool Filled_Progression= kFALSE;
      bool Filled_ProgressionD= kFALSE;

      Int_t IsoCut_Progression;
      TString ResultProgressionName_= "LUT_Progression_";
      TString pt_Progression_= "pt_Progression";
      TString pt_ProgressionD_= "pt_Progression_double";
      
      ResultProgressionName_ += it;
      pt_Progression_ +=it;
      pt_ProgressionD_ +=it;

      TH3F* ResultProgressionName = (TH3F*)gDirectory->Get(ResultProgressionName_.Data());
     /*std::cout<<" ResultProgressionName : "<<ResultProgressionName_<<"\n";
     for(UInt_t kk=0;kk< ResultProgressionName->GetNbinsX();kk++)
     for(UInt_t ll=0;ll< ResultProgressionName->GetNbinsY();ll++)
     for(UInt_t mm=0;mm< ResultProgressionName->GetNbinsZ();mm++)
     {   
	    IsoCut_Progression = ResultProgressionName->GetBinContent(kk,ll,mm);
        std::cout<<"\t\tIsoCut_Progression : "<<kk<<" : "<<ll<<" : "<<mm<<" | "<<IsoCut_Progression<<"\n";
     }*/
      
      for (UShort_t iEG=0; iEG < nEGs; ++iEG) {
	if (egBx[iEG]!=0)   continue;
	float EG_Et  = egEt[iEG];
	short EG_NTT = egNTT[iEG];
	short EG_TowerIEta = egTowerIEta[iEG];
	short EG_Iso_Et = egIsoEt[iEG];
	short EG_Raw_Et = egRawEt[iEG];

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
	if(!Filled_Progression && EG_Iso_Et<=IsoCut_Progression) {
	  ptMap_[pt_Progression_]->Fill(EG_Et);
	  Filled_Progression = kTRUE;
	}
	
	if(nEGs >=2 && !Filled_ProgressionD){ 
	  if(EG_Et > 10. && EG_Iso_Et<=IsoCut_Progression) {
	    if(iEG==0) {
	      Et_Double = std::min(egEt[iEG+1], (EG_Et-10));
	      ptMap_[pt_ProgressionD_]->Fill(Et_Double);
	      Filled_ProgressionD = kTRUE;
	    }	      
	    else  Et_Double = std::min(egEt[iEG+1], (EG_Et-10));
	    ptMap_[pt_ProgressionD_]->Fill(Et_Double);
	    Filled_ProgressionD = kTRUE;
	  }
	}
      }
    }
  }

  double scale = (11.2456 * bunches)/nEntries_;
  
  //Filling Rate Histos ( Single EG + Double EG)
  for(UInt_t i=0;i< ET_MAX;i++) {
    for (auto it :lutProgOptVec_) {
      TString CurrentNameHisto = "pt_Progression";
      TString CurrentNameHisto1= "rate_Progression";
      TString CurrentNameHistoD = "pt_Progression_double";
      TString CurrentNameHisto1D = "rate_Progression_double";
      
      CurrentNameHisto += it;
      CurrentNameHisto1 += it;
      CurrentNameHistoD += it;
      CurrentNameHisto1D += it;
      
      rateMap_[CurrentNameHisto1]->SetBinContent(i+1, ptMap_[CurrentNameHisto]->Integral(i+1,ET_MAX)*scale); //Single EG
      rateMap_[CurrentNameHisto1D]->SetBinContent(i+1, ptMap_[CurrentNameHistoD]->Integral(i+1,ET_MAX)*scale); // Double EG

    if(i==10 or i==22 or i==28) 
    {
     std::cout<<CurrentNameHisto1<<" Et "<<i<<ptMap_[CurrentNameHisto]->Integral(i+1,ET_MAX)<<"*"<<scale<<" = "<<ptMap_[CurrentNameHisto]->Integral(i+1,ET_MAX)*scale<<"\n";
    }
  }
  }
  
  optionsFile_ = new TFile(optionsFileName_.c_str(), "READ");
  optionsFile_->cd("Step2Histos");     
  
  //Event loop for turnon
  if (fChain1 == 0) return;
  nEntries1_ = fChain1->GetEntriesFast();
  Long64_t nbytes1 = 0, nb1 = 0;
  double sum=0;
  std::cout << "Total available  Entries For Efficiency " << nEntries1_ << std::endl;
  if(maxEntriesForEfficiency_ >0 ) nEntries1_= maxEntriesForEfficiency_ < nEntries1_ ? maxEntriesForEfficiency_ : nEntries1_;
  std::cout<<"Processing a total of "<<nEntries1_<<" Entries \n";
  t_start = std::chrono::high_resolution_clock::now();
  
  for (Long64_t jentry=0; jentry< nEntries1_; jentry++) {
    Long64_t ientry = fChain1->LoadTree(jentry);
    if (ientry < 0) break;
    nb1 = fChain1->GetEntry(jentry);   nbytes1 += nb1;
    
    if(jentry%reportEvery_ == 0 )
      {
	    t_end = std::chrono::high_resolution_clock::now();
	    std::cout<<"Processing Entry in event loop (TurnOn) : "<<jentry<<" / "<<nEntries1_<<"  [ "<<100.0*jentry/nEntries1_<<"  % ]  "
		     << " Elapsed time : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()/1000.0
		     <<"  Estimated time left : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()*( nEntries1_ - jentry)/(1e-9 + jentry)* 0.001
		     <<std::endl;
      }
    if(!( isProbeLoose==1 && fabs(eleProbeEta) < 2.5  && sqrt(pow(eleProbeEta-eleTagEta,2)+pow(eleProbePhi-eleTagPhi,2))>0.6)) continue;
    if(l1tEmuRawEt < 0.) continue;
    pT_all->Fill(eleProbeSclEt);
    sum++;
    
    std::map<short, short>::iterator EtaPos = lutMapEta.find(abs(l1tEmuTowerIEta));
    if (EtaPos != lutMapEta.end()) in_compressediEta = EtaPos->second;
    else in_compressediEta = nBinsIEta-1;
    
    std::map<short, short>::iterator EtPos = lutMapEt.find(l1tEmuRawEt);
    if (EtPos != lutMapEt.end()) in_compressediEt = EtPos->second;
    else in_compressediEt = nBinsIEt-1;

    std::map<short, short>::iterator NTTPos = lutMapNTT.find(l1tEmuNTT);
    if (NTTPos != lutMapNTT.end()) in_compressedNTT = NTTPos->second;
    else in_compressedNTT = nBinsNTT-1;
  

    Int_t IsoCut_Progression;
    for(UInt_t e = 5 ; e <= 55; e += 1) {
      for (auto it :lutProgOptVec_) {
	TString ResultProgressionName_= "LUT_Progression_";
        TString PtPassName_= "pt_pass_option";
	
	ResultProgressionName_ += it;
	PtPassName_ += it;
	PtPassName_ +="_Et_";
	PtPassName_ += std::to_string(e);
	
	TH3F* ResultProgressionName = (TH3F*)gDirectory->Get(ResultProgressionName_.Data());
	IsoCut_Progression = ResultProgressionName->GetBinContent(in_compressediEta+1,in_compressediEt+1,in_compressedNTT+1);
	
	if(l1tEmuPt >= e && l1tEmuIsoEt <= IsoCut_Progression)	  pt_pass_Map_[PtPassName_]->Fill(eleProbeSclEt);
	
	//For Nvtx with Iso
	if(eleProbeSclEt>30) {
	  if(l1tEmuPt >= e && l1tEmuIsoEt <= IsoCut_Progression) { 
	    TString NvtxPassIso_ ="Nvtx_Pass_Iso_";
	    NvtxPassIso_+= it;
	    NvtxPassIso_+="_Et_";
	    NvtxPassIso_+=std::to_string(e);
	    
	    TString EtaPassIso_ ="Eta_Pass_Iso_";
            EtaPassIso_+= it;
            EtaPassIso_+="_Et_";
            EtaPassIso_+=std::to_string(e);
	    
	    Nvtx_pass_Map_[NvtxPassIso_]->Fill(Nvtx,1);
	    eta_pass_Map_[EtaPassIso_]->Fill(eleProbeEta,1);
	  }
	  else {
	    TString NvtxFailIso_ ="Nvtx_Fail_Iso_";
            NvtxFailIso_+= it;
            NvtxFailIso_+="_Et_";
            NvtxFailIso_+=std::to_string(e);
	    
            TString EtaFailIso_ ="Eta_Fail_Iso_";
            EtaFailIso_+= it;
            EtaFailIso_+="_Et_";
            EtaFailIso_+=std::to_string(e);


	    Nvtx_fail_Map_[NvtxFailIso_]->Fill(Nvtx,1);
	    eta_fail_Map_[EtaFailIso_]->Fill(eleProbeEta,1);
	  }
	}
      }    //Option loop closing
      
      // if(e==24 && l1tEmuPt >= e) PtPass_inclusive->Fill(eleProbeSclEt); 
      
      //For Vertex  without Iso
	if(eleProbeSclEt>30) {                                                                                     
	  if(l1tEmuPt >= e) {
	    TString NvtxPass_ ="Nvtx_Pass_";
	    NvtxPass_+="_Et_";
	    NvtxPass_+=std::to_string(e);
	    
	    TString EtaPass_ ="Eta_Pass_";
	    EtaPass_+="_Et_";
	    EtaPass_+=std::to_string(e);
	    
	    Nvtx_pass_Map_[NvtxPass_]->Fill(Nvtx,1);
	    eta_pass_Map_[EtaPass_]->Fill(eleProbeEta,1);
	  }
	  else {
	    TString NvtxFail_ ="Nvtx_Fail_";
	    NvtxFail_+="_Et_";
	    NvtxFail_+=std::to_string(e);
	    
	    TString EtaFail_ ="Eta_Fail_";
	    EtaFail_+="_Et_";
	    EtaFail_+=std::to_string(e);
	    
	    Nvtx_fail_Map_[NvtxFail_]->Fill(Nvtx,1);
	    eta_fail_Map_[EtaFail_]->Fill(eleProbeEta,1);
	  }
	}
    }
    // using Run2 LUT [***Make it generalised]             
    if(hasL1Emu_26==1) pt_pass_Map_Run2_["pt_pass_Et_26"]->Fill(eleProbeSclEt);
    if(hasL1Emu_24==1)   pt_pass_Map_Run2_["pt_pass_Et_24"]->Fill(eleProbeSclEt);
    if(hasL1Emu_22==1) pt_pass_Map_Run2_["pt_pass_Et_22"]->Fill(eleProbeSclEt);
    
    if(hasL1Emu_looseiso26==1) pt_pass_Map_Run2_["pt_pass_Et_26_loose"]->Fill(eleProbeSclEt);
    if(hasL1Emu_looseiso24==1) pt_pass_Map_Run2_["pt_pass_Et_24_loose"]->Fill(eleProbeSclEt);
    if(hasL1Emu_looseiso22==1) pt_pass_Map_Run2_["pt_pass_Et_22_loose"]->Fill(eleProbeSclEt);
    
    if(hasL1Emu_tightiso26==1) pt_pass_Map_Run2_["pt_pass_Et_26_tight"]->Fill(eleProbeSclEt);
    if(hasL1Emu_tightiso24==1) pt_pass_Map_Run2_["pt_pass_Et_24_tight"]->Fill(eleProbeSclEt);
    if(hasL1Emu_tightiso22==1) pt_pass_Map_Run2_["pt_pass_Et_22_tight"]->Fill(eleProbeSclEt);
  } //End of Event loop for turnon

  //Turnons/Nvtx_Eff/Eta_Eff  for various Ets
  for(UInt_t e = 5 ; e <= 55 ; e += 1) {
    for (auto it :lutProgOptVec_) {
      TString PtPassName_= "pt_pass_option";
      TString turnOn_Option_="turnOn_Option";
      
      PtPassName_ += it;
      PtPassName_ +="_Et_";
      PtPassName_ += std::to_string(e);
      
      turnOn_Option_ += it;
      turnOn_Option_ +="_Et_";
      turnOn_Option_ += std::to_string(e);
      outputFile_->cd();
      
      td3->cd();
      turnOn_Map_[turnOn_Option_]=new TGraphAsymmErrors( pt_pass_Map_[PtPassName_],pT_all,"cp"); 
      turnOn_Map_[turnOn_Option_]->Write();

      double acceptance = (pt_pass_Map_[PtPassName_]->GetEntries())/(pT_all->GetEntries()) ;
      std::cout<<"accc"<<acceptance<<"pt_all"<<pT_all->GetEntries()<<"pass"<<pt_pass_Map_[PtPassName_]->GetEntries()<<std::endl;
      
      double rate = FindRate(it, e);
      double rateD = FindRateD(it, e);
      if(e==24) std::cout<<turnOn_Option_ <<"  24 :"<<rate<<std::endl;
      if(e==10) std::cout<<turnOn_Option_ <<"  10 :"<<rateD<<std::endl;

      th1fStore["EtForTurnons"]->Fill(turnOn_Map_[turnOn_Option_]->GetName(),e);
      th1fStore["AcceptanceForTurnons"]->Fill(turnOn_Map_[turnOn_Option_]->GetName(),acceptance);
      th1fStore["RateForTurnons"]->Fill(turnOn_Map_[turnOn_Option_]->GetName(),rate);
      th1fStore["RateDForTurnons"]->Fill(turnOn_Map_[turnOn_Option_]->GetName(),rateD);

      //Nvtx/Eta Iso
      TString NvtxPassIso_ ="Nvtx_Pass_Iso_";
      NvtxPassIso_+= it;
      NvtxPassIso_+="_Et_";
      NvtxPassIso_+=std::to_string(e);

      TString EtaPassIso_ ="Eta_Pass_Iso_";
      EtaPassIso_+= it;
      EtaPassIso_+="_Et_";
      EtaPassIso_+=std::to_string(e);
      TString NvtxFailIso_ ="Nvtx_Fail_Iso_";
      NvtxFailIso_+= it;
      NvtxFailIso_+="_Et_";
      NvtxFailIso_+=std::to_string(e);

      TString EtaFailIso_ ="Eta_Fail_Iso_";
      EtaFailIso_+= it;
      EtaFailIso_+="_Et_";
      EtaFailIso_+=std::to_string(e);

      TString EtaIsoEff_ ="Eta_Eff_Iso_";
      EtaIsoEff_+= it;
      EtaIsoEff_+="_Et_";
      EtaIsoEff_+=std::to_string(e);

      TString NvtxIsoEff_ ="Nvtx_Eff_Iso_";
      NvtxIsoEff_+= it;
      NvtxIsoEff_+="_Et_";
      NvtxIsoEff_+=std::to_string(e);

      createEfficiencyHistograms(Nvtx_pass_Map_[NvtxPassIso_],Nvtx_fail_Map_[NvtxFailIso_],Nvtx_Eff_Map_[NvtxIsoEff_]);
      createEfficiencyHistograms(eta_pass_Map_[EtaPassIso_],eta_fail_Map_[EtaFailIso_],eta_Eff_Map_[EtaIsoEff_]);
    
    }
    if(e==24){turnon_inclusive=new TGraphAsymmErrors(PtPass_inclusive,pT_all,"cp");
      turnon_inclusive->Write();
    }
    //Nvtx/Eta without Iso
      TString NvtxPass_ ="Nvtx_Pass_";
      NvtxPass_+="_Et_";
      NvtxPass_+=std::to_string(e);

      TString EtaPass_ ="Eta_Pass_";
      EtaPass_+="_Et_";
      EtaPass_+=std::to_string(e);

      TString NvtxFail_ ="Nvtx_Fail_";
      NvtxFail_+="_Et_";
      NvtxFail_+=std::to_string(e);

      TString EtaFail_ ="Eta_Fail_";
      EtaFail_+="_Et_";
      EtaFail_+=std::to_string(e);

      TString NvtxEff_ ="Nvtx_Eff_";
      NvtxEff_+="_Et_";
      NvtxEff_+=std::to_string(e);

      TString EtaEff_ ="Eta_Eff_";
      EtaEff_+="_Et_";
      EtaEff_+=std::to_string(e);

      createEfficiencyHistograms(Nvtx_pass_Map_[NvtxPass_],Nvtx_fail_Map_[NvtxFail_],Nvtx_Eff_Map_[NvtxEff_]);
      createEfficiencyHistograms(eta_pass_Map_[EtaPass_],eta_fail_Map_[EtaFail_],eta_Eff_Map_[EtaEff_]);

  }

  
  // using Run2 LUT                                                                                                                                      
  td4->cd();
  std::vector<int> v = {22, 24, 26};
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

    double acceptance = (pt_pass_Map_Run2_[PtPassName_]->GetEntries())/(pT_all->GetEntries()) ;
    double acceptance_loose = (pt_pass_Map_Run2_[PtPassLooseIsoName_]->GetEntries())/(pT_all->GetEntries()) ;
    double acceptance_tight = (pt_pass_Map_Run2_[PtPassTightIsoName_]->GetEntries())/(pT_all->GetEntries()) ;

    th1fStore["Run2Turnons_Acceptance"]->Fill(turnOn_Map_Run2_[turnOnName_]->GetName(),acceptance);
    th1fStore["Run2Turnons_LooseIso_Acceptance"]->Fill(turnOn_Map_Run2_[turnOnLooseIsoName_]->GetName(),acceptance_loose);
    th1fStore["Run2Turnons_TightIso_Acceptance"]->Fill(turnOn_Map_Run2_[turnOnTightIsoName_]->GetName(),acceptance_tight);
  }
}

void ApplyIsolation::bookHistograms() {
  outputFile_ = new TFile(outputFileName_.c_str(), "RECREATE");
  outputFile_->cd();
  
  for (auto it : lutProgOptVec_) {    
    TString  CurrentNameHisto = "pt_Progression";
    TString CurrentNameHisto1= "rate_Progression";
    CurrentNameHisto += it;
    CurrentNameHisto1 += it;
    
    if(!check_pt_rate_dir) { 
      td = outputFile_->mkdir("pt_progression_for_rate");
      check_pt_rate_dir = true;
    }
    td->cd();
    TH1F* pt_Progression=  new TH1F(CurrentNameHisto, CurrentNameHisto ,ET_MAX ,  0.0 , ET_MAX);
    ptMap_.insert(std::make_pair(CurrentNameHisto,pt_Progression));
    

    TString  CurrentNameHistoD = "pt_Progression_double";
    TString CurrentNameHisto1D= "rate_Progression_double";
    CurrentNameHistoD += it;
    CurrentNameHisto1D += it;
    if(!check_pt_rate_double_dir) {
      tdD = outputFile_->mkdir("pt_progression_for_rate_double");
      check_pt_rate_double_dir = true;
    }
    tdD->cd();
    TH1F* pt_ProgressionD=  new TH1F(CurrentNameHistoD, CurrentNameHistoD ,ET_MAX ,  0.0 , ET_MAX);
    ptMap_.insert(std::make_pair(CurrentNameHistoD,pt_ProgressionD));


    if(!check_rate_dir) {
      td1 = outputFile_->mkdir("rate_progression");
      check_rate_dir = true;
    }
    td1->cd();   
    TH1F* rate_Progression = new TH1F(CurrentNameHisto1, CurrentNameHisto1, ET_MAX ,  0.0 , ET_MAX );
    rateMap_.insert(std::make_pair(CurrentNameHisto1,rate_Progression));

    if(!check_rate_double_dir) {
      td1D = outputFile_->mkdir("rate_progression_double");
      check_rate_double_dir = true;
    }
    td1D->cd();
    TH1F* rate_ProgressionD = new TH1F(CurrentNameHisto1D, CurrentNameHisto1D, ET_MAX ,  0.0 , ET_MAX );
    rateMap_.insert(std::make_pair(CurrentNameHisto1D,rate_ProgressionD));

    
    for(UInt_t e = 5 ; e <= 55 ; e += 1) {
      TString PtPassName_= "pt_pass_option";
      TString turnOn_Option_="turnOn_Option";
      
      PtPassName_ += it;
      PtPassName_ +="_Et_";
      PtPassName_ += std::to_string(e);  
      
      if(!check_pt_turn_on_dir) {
	td2 = outputFile_->mkdir("pt_progression_for_turnon");
	check_pt_turn_on_dir = true;
      }
      td2->cd();
      TH1F*PtPass= new TH1F(PtPassName_, PtPassName_, nBins_fine,xEdges_fine);
      pt_pass_Map_.insert(std::make_pair(PtPassName_, PtPass));
      
      turnOn_Option_ += it;
      turnOn_Option_ +="_Et_";
      turnOn_Option_ += std::to_string(e);
      
      TGraphAsymmErrors* turnOn_Option;
      turnOn_Map_.insert(std::make_pair(turnOn_Option_,turnOn_Option));
    }
  }
  td2->cd();
  pT_all = new TH1F("pT_all","pT_all",nBins_fine,xEdges_fine);
  outputFile_->cd();

  th1fStore["EtForTurnons"]   = new TH1F("EtForTurnons","",3,0.0,3.0);
  th1fStore["EtForTurnons"]->SetCanExtend(TH1::kAllAxes);
  th1fStore["AcceptanceForTurnons"]   = new TH1F("AcceptanceForTurnons","",3,0.0,3.0);
  th1fStore["AcceptanceForTurnons"]->SetCanExtend(TH1::kAllAxes);
  th1fStore["RateForTurnons"]   = new TH1F("RateForTurnons","",3,0.0,3.0);
  th1fStore["RateForTurnons"]->SetCanExtend(TH1::kAllAxes);
  th1fStore["RateDForTurnons"]   = new TH1F("RateDForTurnons","",3,0.0,3.0);
  th1fStore["RateDForTurnons"]->SetCanExtend(TH1::kAllAxes);


  if(!check_turn_on_dir) {
    td3 = outputFile_->mkdir("turn_on_progression");
    check_turn_on_dir = true;
  }

  //using Run2 LUT
  std::vector<int> v = {22, 24, 26};
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
    
    if(!check_turn_on_dir_Run2) {
      td4 = outputFile_->mkdir("turnon_progression_Run2");
      check_turn_on_dir_Run2 = true;
    }
    td4->cd();
    TH1F*PtPass= new TH1F(PtPassName_, PtPassName_, nBins_fine,xEdges_fine);
    TH1F*PtPassLoose= new TH1F(PtPassLooseIsoName_, PtPassLooseIsoName_, nBins_fine,xEdges_fine);
    TH1F*PtPassTight= new TH1F(PtPassTightIsoName_, PtPassTightIsoName_, nBins_fine,xEdges_fine);
    
    pt_pass_Map_Run2_.insert(std::make_pair(PtPassName_, PtPass));
    pt_pass_Map_Run2_.insert(std::make_pair(PtPassLooseIsoName_, PtPassLoose));
    pt_pass_Map_Run2_.insert(std::make_pair(PtPassTightIsoName_, PtPassTight));
    
    
    TGraphAsymmErrors* turnOn;
    TGraphAsymmErrors* turnOnLoose;
    TGraphAsymmErrors* turnOnTight;

    turnOn_Map_Run2_.insert(std::make_pair(turnOnName_,turnOn));
    turnOn_Map_Run2_.insert(std::make_pair(turnOnLooseIsoName_,turnOnLoose));
    turnOn_Map_Run2_.insert(std::make_pair(turnOnTightIsoName_,turnOnTight));
  }

  PtPass_inclusive= new TH1F("PtPass_inclusive"," PtPass_inclusive", nBins_fine,xEdges_fine);
 
  th1fStore["Run2Turnons_Acceptance"]   = new TH1F("Run2Turnons_Acceptance","",3,0.0,3.0);
  th1fStore["Run2Turnons_Acceptance"]->SetCanExtend(TH1::kAllAxes);
  th1fStore["Run2Turnons_LooseIso_Acceptance"]   = new TH1F("Run2Turnons_LooseIso_Acceptance","",3,0.0,3.0);
  th1fStore["Run2Turnons_LooseIso_Acceptance"]->SetCanExtend(TH1::kAllAxes);
  th1fStore["Run2Turnons_TightIso_Acceptance"]   = new TH1F("Run2Turnons_TightIso_Acceptance","",3,0.0,3.0);
  th1fStore["Run2Turnons_TightIso_Acceptance"]->SetCanExtend(TH1::kAllAxes);
  

  //For Vertex
  for(UInt_t e = 5 ; e <= 55 ; e += 1) {

      TString NvtxPass_ ="Nvtx_Pass_";
      NvtxPass_+="_Et_";
      NvtxPass_+=std::to_string(e);

      TString NvtxFail_ ="Nvtx_Fail_";
      NvtxFail_+="_Et_";
      NvtxFail_+=std::to_string(e);

      TString NvtxEff_ ="Nvtx_Eff_";
      NvtxEff_+="_Et_";
      NvtxEff_+=std::to_string(e);

    if(!check_Nvtx) {
      td5 = outputFile_->mkdir("Nvtx");
      check_Nvtx = true;
    }
    td5->cd();
    TH1D*thNvtxP = new TH1D(NvtxPass_, NvtxPass_, 100, -0.5, 99.5);
    TH1D*thNvtxF = new TH1D(NvtxFail_, NvtxFail_, 100, -0.5, 99.5);
    TH1D*thNvtxE = new TH1D(NvtxEff_, NvtxEff_, 100, -0.5, 99.5);
    Nvtx_pass_Map_.insert(std::make_pair(NvtxPass_, thNvtxP));
    Nvtx_fail_Map_.insert(std::make_pair(NvtxFail_,thNvtxF));
    Nvtx_Eff_Map_.insert(std::make_pair(NvtxEff_,thNvtxE));
  }


//For Vertex Iso            
  for(UInt_t e = 5 ; e <= 55 ; e += 1) {
    for (auto it :lutProgOptVec_) {

      TString NvtxPassIso_ ="Nvtx_Pass_Iso_";
      NvtxPassIso_+= it;
      NvtxPassIso_+="_Et_";
      NvtxPassIso_+=std::to_string(e);

      TString NvtxFailIso_ ="Nvtx_Fail_Iso_";
      NvtxFailIso_+= it;
      NvtxFailIso_+="_Et_";
      NvtxFailIso_+=std::to_string(e);

      TString NvtxIsoEff_ ="Nvtx_Eff_Iso_";
      NvtxIsoEff_+= it;
      NvtxIsoEff_+="_Et_";
      NvtxIsoEff_+=std::to_string(e);

      if(!check_Nvtx_Iso) {
	td6 = outputFile_->mkdir("Nvtx_Iso");
	check_Nvtx_Iso = true;
      }
      td6->cd();
      TH1D*thNvtxPI = new TH1D(NvtxPassIso_, NvtxPassIso_, 100, -0.5, 99.5);
      TH1D*thNvtxFI = new TH1D(NvtxFailIso_, NvtxFailIso_, 100, -0.5, 99.5);
      TH1D*thNvtxEI = new TH1D(NvtxIsoEff_, NvtxIsoEff_, 100, -0.5, 99.5);
      Nvtx_pass_Map_.insert(std::make_pair(NvtxPassIso_,thNvtxPI));
      Nvtx_fail_Map_.insert(std::make_pair(NvtxFailIso_,thNvtxFI));
      Nvtx_Eff_Map_.insert(std::make_pair(NvtxIsoEff_,thNvtxEI));
    }
  }
  
  //For Eta
  for(UInt_t e = 5 ; e <= 55 ; e += 1) {
      
      TString EtaPass_ ="Eta_Pass_";
      EtaPass_+="_Et_";
      EtaPass_+=std::to_string(e);
      
      TString EtaFail_ ="Eta_Fail_";
      EtaFail_+="_Et_";
      EtaFail_+=std::to_string(e);

      TString EtaEff_ ="Eta_Eff_";
      EtaEff_+="_Et_";
      EtaEff_+=std::to_string(e);

      if(!check_eta) {
	td7 = outputFile_->mkdir("Eta");
	check_eta = true;
      }
      td7->cd();

      TH1D*thetaP = new TH1D(EtaPass_, EtaPass_, 10, -2.5, 2.5);
      TH1D*thetaF = new TH1D(EtaFail_, EtaFail_, 10, -2.5, 2.5);
      TH1D*thetaE = new TH1D(EtaEff_, EtaEff_, 10, -2.5, 2.5);
      
      eta_pass_Map_.insert(std::make_pair(EtaPass_, thetaP));
      eta_fail_Map_.insert(std::make_pair(EtaFail_,thetaF));
      eta_Eff_Map_.insert(std::make_pair(EtaEff_,thetaE));
  }
  
  
  //For Eta Iso                                                                                                    
  for(UInt_t e = 5 ; e <= 55 ; e += 1) {
    for (auto it :lutProgOptVec_) {

      TString EtaPassIso_ ="Eta_Pass_Iso_";
      EtaPassIso_+= it;
      EtaPassIso_+="_Et_";
      EtaPassIso_+=std::to_string(e);

      TString EtaFailIso_ ="Eta_Fail_Iso_";
      EtaFailIso_+= it;
      EtaFailIso_+="_Et_";
      EtaFailIso_+=std::to_string(e);

      TString EtaIsoEff_ ="Eta_Eff_Iso_";
      EtaIsoEff_+= it;
      EtaIsoEff_+="_Et_";
      EtaIsoEff_+=std::to_string(e);

      if(!check_eta_Iso) {
        td8 = outputFile_->mkdir("Eta_Iso");
        check_eta_Iso = true;
      }
      td8->cd();

      TH1D*thetaPI = new TH1D(EtaPassIso_, EtaPassIso_, 10, -2.5, 2.5);
      TH1D*thetaFI = new TH1D(EtaFailIso_, EtaFailIso_, 10, -2.5, 2.5);
      TH1D*thetaEI = new TH1D(EtaIsoEff_, EtaIsoEff_, 10, -2.5, 2.5);

      eta_pass_Map_.insert(std::make_pair(EtaPassIso_,thetaPI));
      eta_fail_Map_.insert(std::make_pair(EtaFailIso_,thetaFI));
      eta_Eff_Map_.insert(std::make_pair(EtaIsoEff_,thetaEI));
    }
  }


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
      else if(key=="WeightFileName") weightFileName_ = value.c_str();
      else if (key=="OptionsFileName") optionsFileName_ = value.c_str();
      else if (key=="OutputFileName")  outputFileName_ = value.c_str();
      else if (key=="EtLUTFileName")    readLUTTable(value,nBinsIEt, lutMapEt);
      else if (key=="EtaLUTFileName")   readLUTTable(value,nBinsIEta, lutMapEta);
      else if (key=="NTTLUTFileName")   readLUTTable(value,nBinsNTT, lutMapNTT);
      else if (key=="nTTRange") nTTRange_ = std::stoi(value.c_str());
      else if (key=="ReportEvery")    reportEvery_ = atoi(value.c_str());
      else if (key=="MaxEntriesForEfficency") maxEntriesForEfficiency_ = atoi(value.c_str());
      else if (key=="MaxEntriesForRate")   maxEntriesForRate_ = atoi(value.c_str());
      else if (key=="FixedRate") frate_ = std::stof(value.c_str());
      else if (key=="FixedRateD") frateD_ = std::stof(value.c_str());

      else if (key=="LUTProgressionOptions")
        {
	  std::string tmp_string = value;
	  std::vector<std::string> tmp_vec;
          tokenize(tmp_string,tmp_vec,",");
          for (auto it : tmp_vec) {
            lutProgOptVec_.push_back(it);
          }
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
  fChain_1->SetMakeClass(1);
  fChain_1->SetBranchAddress("nPV_True", &nPV_True);

  //for turnon
  fChain1->SetBranchAddress("l1tEmuPt", &l1tEmuPt);
  fChain1->SetBranchAddress("l1tEmuNTT",&l1tEmuNTT);
  fChain1->SetBranchAddress("l1tEmuRawEt",&l1tEmuRawEt);
  fChain1->SetBranchAddress("l1tEmuTowerIEta",&l1tEmuTowerIEta);
  fChain1->SetBranchAddress("eleProbeSclEt",&eleProbeSclEt);
  fChain1->SetBranchAddress("l1tEmuIsoEt",&l1tEmuIsoEt);

  fChain1->SetBranchAddress("eleProbeEta",&eleProbeEta);
  fChain1->SetBranchAddress("eleProbePhi",&eleProbePhi);
  fChain1->SetBranchAddress("eleTagEta",&eleTagEta);
  fChain1->SetBranchAddress("eleTagPhi",&eleTagPhi);
  fChain1->SetBranchAddress("isProbeLoose",&isProbeLoose);
  fChain1->SetBranchAddress("Nvtx",&Nvtx);

  
  fChain1->SetBranchAddress("hasL1Emu_26",&hasL1Emu_26);
  fChain1->SetBranchAddress("hasL1Emu_24",&hasL1Emu_24);
  fChain1->SetBranchAddress("hasL1Emu_22",&hasL1Emu_22);

  fChain1->SetBranchAddress("hasL1Emu_looseiso26",&hasL1Emu_looseiso26);
  fChain1->SetBranchAddress("hasL1Emu_looseiso24",&hasL1Emu_looseiso24);
  fChain1->SetBranchAddress("hasL1Emu_looseiso22",&hasL1Emu_looseiso22);

  fChain1->SetBranchAddress("hasL1Emu_tightiso26",&hasL1Emu_tightiso26);
  fChain1->SetBranchAddress("hasL1Emu_tightiso24",&hasL1Emu_tightiso24);
  fChain1->SetBranchAddress("hasL1Emu_tightiso22",&hasL1Emu_tightiso22);



}

double ApplyIsolation::BetaInverse(double x,double p, double q)
{
  double result(0.0);
  double dy = 0.001;
  double eMin = 100;
  for(int i=0;i<1000;i++){
    double y = i*dy;
    double e = fabs(TMath::BetaIncomplete(y,p,q)-x);
    if (e<eMin)
      {
        eMin = e;
        result = y;
      }
  }
  return result;
}

double ApplyIsolation::FindRate(std::string it, UInt_t e) { 
  TString CurrentNameHisto1= "rate_Progression";
  CurrentNameHisto1 += it;
  double rate = rateMap_[CurrentNameHisto1]->GetBinContent(e);
  return rate;
}

double ApplyIsolation::FindRateD(std::string it, UInt_t e) {
  TString CurrentNameHisto1D= "rate_Progression_double";
  CurrentNameHisto1D += it;
  double rateD =  rateMap_[CurrentNameHisto1D]->GetBinContent(e);
  return rateD;
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


