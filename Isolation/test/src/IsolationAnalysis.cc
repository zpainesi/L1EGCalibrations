#include <string>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TObjString.h>
#include <iostream>
#include <fstream>
#include "IsolationAnalysis.h"

IsolationAnalysis::IsolationAnalysis(const std::string& inputFileName){
  reportEvery=5000;
  maxEntries=-1;
  readParameters(inputFileName);
  if (ntupleFileName_.size() == 0) {
    std::cout << " Inputfile list missing !!!" << ntupleFileName_ << std::endl;
    return;
  }

  accessTree(ntupleFileName_);
  
  assert(fChain);
  
  bookedHistograms_ = false;
 
  short last = -1;
  for (auto ieta : lutMapIEta_)
    {
      short val = ieta.second;
      if (val != last) {
	lutIEtaVec_.push_back(ieta.first);
      }
      last = val; 
    }
  lutIEtaVec_.push_back(31);
  std::cout << "LUTIEtaBins: ";
  for (auto ieta : lutIEtaVec_) std::cout << ieta << " , " ;
  std::cout << " (Size " << lutIEtaVec_.size() << ")"<< std::endl;
  last = -1;
  for (auto iet : lutMapIEt_)
    {
      short val = iet.second;
      if (val != last) 
	{
	  lutIEtVec_.push_back(iet.first);
	}
      last = val; 
    }
  lutIEtVec_.push_back(255);
  std::cout << "LUTIEtBins: ";
  for (auto iet : lutIEtVec_) std::cout << iet << " , " ;
  std::cout << " (Size " << lutIEtVec_.size()<< ")"<< std::endl;
  last = -1;
  for (auto intt : lutMapnTT_)
    {
      short val = intt.second;
      if (val != last) {
	lutnTTVec_.push_back(intt.first);
      }
      last = val; 
    }
  lutnTTVec_.push_back(255);
  std::cout << "LUTnTTBins: ";
  for (auto intt : lutnTTVec_) std::cout << intt << " , " ;
  std::cout << " (Size " << lutnTTVec_.size() <<  ")" << std::endl;
}


IsolationAnalysis::~IsolationAnalysis(){
  /*  for (auto it : Histos_PerBin) {
      if (it.second) it.second->Delete();
      }
      Histos_PerBin.clear();*/
}
void IsolationAnalysis::accessTree(std::string & input_filelist) {
  std::ifstream myFile;
  myFile.open(input_filelist.c_str(), std::ios::in);
  if (!myFile) {
    std::cout << "Input File: " << input_filelist << " could not be opened!" << std::endl;
    fChain = 0;
  } else {
    fChain = new TChain("Ntuplizer/TagAndProbe");
    static constexpr int BUF_SIZE = 256;
    char buf[BUF_SIZE];
    while (myFile.getline(buf, BUF_SIZE, '\n')) {  // Pops off the newline character
      std::string line(buf); 
      fChain->AddFile(line.c_str(),-1);
      std::cout << "Adding file " << line << " Entries " << fChain->GetEntries() <<  std::endl;
    }
  }
}

void IsolationAnalysis::analyse() {
  
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  std::cout<<"Total Available Entries : "<<nentries<<std::endl;
  Long64_t nbytes = 0, nb = 0;
  if(maxEntries > -1) nentries = maxEntries < nentries ? maxEntries :nentries ;
  std::cout<<"Processing a total of "<<nentries<<" Entries \n";

  auto t_start = std::chrono::high_resolution_clock::now();
  auto t_end = std::chrono::high_resolution_clock::now();

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    if (jentry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if(jentry%reportEvery == 0 )
       {
             t_end = std::chrono::high_resolution_clock::now();
             std::cout<<"Loop 1/3 : Processing Entry in event loop : "<<jentry<<" / "<<nentries<<"  [ "<<100.0*jentry/nentries<<"  % ]  "
                      << " Elapsed time : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()/1000.0
                      <<"  Estimated time left : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()*( nentries - jentry)/(1e-9 + jentry)* 0.001
                      <<std::endl;

       } 

    if( isProbeLoose==0 ) continue;
    if( fabs(eleProbeEta) >= 2.5) continue;
    if( sqrt(pow(eleProbeEta-eleTagEta,2)+pow(eleProbePhi-eleTagPhi,2)) < 0.6 ) continue;
    hprof_IEt->Fill(et, iso);
    hprof_IEta->Fill(eta, iso);
    hprof_nTT->Fill(ntt, iso);
    
    if(et < 0 )  continue;

    TString Name_Histo = getHistoName(eta,et, ntt);
    std::map<TString, TH1F*>::iterator iPos = Histos_PerBin.find(Name_Histo);
    if (iPos != Histos_PerBin.end()) iPos->second->Fill(iso);
    if (eta > 31.5) {
      pt_large_eta->Fill(et);
      eta_large_eta->Fill(eta);
      nTT_large_eta->Fill(ntt);
      iso_large_eta->Fill(iso);
    }
  }

  std::cout << "  Total Number of Histograms " << Histos_PerBin.size() << std::endl;
  Int_t NumberOfHistosWithLowStats = 0;

  Int_t c=0;
  t_start = std::chrono::high_resolution_clock::now();
  t_end = std::chrono::high_resolution_clock::now();

  for (auto it : Histos_PerBin) {
    TString Name_Histo = it.first;
    TH1F* th = it.second;
    
    short ibin_eta;
    short ibin_et;
    short ibin_ntt;    
    bool name_flag = getHistoBin(Name_Histo, ibin_eta, ibin_et, ibin_ntt);
    if (!name_flag) {
      std::cout << " ===> " <<  Name_Histo.Data() << ibin_eta << " " << ibin_et<<  " " << ibin_ntt << std::endl;
      continue;
    }
    if(th->GetEntries()<40){
      std::cout<<"Low stat Bin : "<<ibin_et<<" , "<<ibin_eta<<" , "<< ibin_ntt <<" , "<< th->GetEntries()<<"\n";
      NumberOfHistosWithLowStats++;
    }
    for(UInt_t iEff = 1 ; iEff < 101 ; ++iEff){
      Float_t Efficiency = 0.01*iEff;
      
      for(UInt_t iIso = 0 ; iIso < 100 ; ++iIso) {
	if(th->Integral(1,iIso+1)/th->Integral(1,100+1)>=Efficiency){
	  if (IsoCut_PerEfficiency_PerBin[iEff][Name_Histo] == -1) {
	    IsoCut_PerEfficiency_PerBin[iEff][Name_Histo] =  iIso;
	    IsoCut_PerBin[iEff]->SetBinContent(ibin_eta+1,ibin_et+1,ibin_ntt+1,iIso);
	  }
	}
      }
    }
    c++;
    if(c%200 == 0)
       {
             t_end = std::chrono::high_resolution_clock::now();
             std::cout<<"Loop 2/3 : Processing Hist in event loop : "<<c<<" / "<<Histos_PerBin.size()<<"  [ "<<100.0*c/Histos_PerBin.size()<<"  % ]  "
                      << " Elapsed time : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()/1000.0
                      <<"  Estimated time left : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()*( Histos_PerBin.size() - c)/(1e-9 + c )* 0.001
                      <<std::endl;
       }
  }
  
  t_start = std::chrono::high_resolution_clock::now();
  for (Long64_t jentry=0; jentry<nentries;jentry++) {

      if(jentry%reportEvery == 0 )
       {
             t_end = std::chrono::high_resolution_clock::now();
             std::cout<<"Loop 3/3 : Processing Entry in event loop : "<<jentry<<" / "<<nentries<<"  [ "<<100.0*jentry/nentries<<"  % ]  "
                      << " Elapsed time : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()/1000.0
                      <<"  Estimated time left : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()*( nentries - jentry)/(1e-9 + jentry)* 0.001
                      <<std::endl;
       }
    if (jentry < 0) break;

    nb = fChain->GetEntry(jentry);   nbytes += nb;
    
    if( isProbeLoose==0 ) continue;
    if( fabs(eleProbeEta) >= 2.5) continue;
    if( sqrt(pow(eleProbeEta-eleTagEta,2)+pow(eleProbePhi-eleTagPhi,2)) < 0.6 ) continue;
    
    pt_all->Fill(et);
    if ( et < 0) continue;     
    
    eta_all->Fill(eta);
    nTT_all->Fill(ntt);
    
    TString Name_Histo = getHistoName(eta,et, ntt);  
    
    short ibin_eta;
    short ibin_et;
    short ibin_ntt;      
    
    getHistoBin(Name_Histo, ibin_eta,  ibin_et, ibin_ntt);
    
    std::map<TString,TH1F*>::iterator iPos = Histos_PerBin.find(Name_Histo);
    if (iPos == Histos_PerBin.end()) continue;
    TH1F* th = iPos->second;
    if (!th->GetEntries()) continue;
    float eff_at_isolation = th->Integral(1,iso)/th->Integral(1,100+1);
    
    for(UInt_t iEff = 1 ; iEff < 101 ; ++iEff) 	{
      std::map<Int_t,std::map<TString,Int_t>>::iterator jPos = IsoCut_PerEfficiency_PerBin.find(iEff);
      
      if (jPos == IsoCut_PerEfficiency_PerBin.end()) continue;
      std::map<TString,Int_t>::iterator kPos = jPos->second.find(Name_Histo);
      if (kPos == jPos->second.end())  continue;
      if(iso<=kPos->second){
	eta_pass_efficiency[iEff]->Fill(eta);
	pt_pass_efficiency[iEff]->Fill(et);
	nTT_pass_efficiency[iEff]->Fill(ntt);
      }
      if(iso<=IsoCut_PerBin[iEff]->GetBinContent(ibin_eta+1, ibin_et+1, ibin_ntt+1)) pt_pass_efficiency_TH3[iEff]->Fill(et);
      //      if (eff_at_isolation <= iEff*0.01  ){
      //	eta_pass_efficiency[iEff]->Fill(eta);
      //	pt_pass_efficiency[iEff]->Fill(et);
      //	nTT_pass_efficiency[iEff]->Fill(ntt);
      //	
      //	IsoCut_PerBin[iEff]->SetBinContent(ibin_eta+1, ibin_et+1, ibin_ntt+1, isolation);
      //	}
    }
  }
  outputFile_->cd();  
  outputFile_->cd("Step1Histos");
  for(UInt_t iEff = 0 ; iEff < 101 ; ++iEff) {
    TGraphAsymmErrors* temp_histo1 = new TGraphAsymmErrors(pt_pass_efficiency[iEff],pt_all,"cp");
    TGraphAsymmErrors* temp_histo_TH3 = new TGraphAsymmErrors(pt_pass_efficiency_TH3[iEff],pt_all,"cp");
    
    TGraphAsymmErrors* temp_histo2 = new TGraphAsymmErrors(eta_pass_efficiency[iEff],eta_all,"cp");
    
    TGraphAsymmErrors* temp_histo3 = new TGraphAsymmErrors(nTT_pass_efficiency[iEff],nTT_all,"cp");
    
    temp_histo1->Write();
    temp_histo_TH3->Write();
    temp_histo2->Write();
    temp_histo3->Write();    
  }
 
  std::cout<<"NumberOfHistosWithLowStats/Tot = "<<NumberOfHistosWithLowStats<<"/"
	   << Histos_PerBin.size()<<std::endl;
}

void IsolationAnalysis::bookHistograms(std::string option) 
{  

  outputFile_ = new TFile(outputFileName_.c_str(), "RECREATE");
  outputFile_->cd();

  if(option == "do_1" || option == "do_all") {
    TDirectory* td_step1 =   outputFile_->mkdir("Step1Histos");
    if (td_step1) td_step1->cd();
    for(UInt_t i = 0 ; i < nBinsIEta ; ++i) {
      for(UInt_t j = 0 ; j < nBinsIEt ; ++j) {
	for(UInt_t k = 0 ; k < nBinsnTT ; ++k) {
	  TString Name_Histo = "Hist_";
	  Name_Histo +=  i;
	  Name_Histo += "_";
	  Name_Histo += j;
	  Name_Histo += "_";
	  Name_Histo += k;
	  
	  TH1F* temp_histo = new TH1F(Name_Histo.Data(),Name_Histo.Data(),100,0.,100.);
	  Histos_PerBin.insert({Name_Histo,temp_histo});
	}
      }
    }
    hprof_IEt  = new TProfile("hprof_IEt","Profile L1_Iso vs. L1_IEt",100,0.,200.,0,20);
    hprof_IEta  = new TProfile("hprof_IEta","Profile L1_Iso vs. L1_IEta",32,0.,32.,0,20);
    hprof_nTT  = new TProfile("hprof_nTT","Profile L1_Iso vs. L1_IEta",180,20.,200.,0,20);
    std::map<TString, Int_t> tempMap;
    for(UInt_t iEff = 0 ; iEff < 101 ; ++iEff){
      
      for (auto it : Histos_PerBin) {
	tempMap.insert({it.first, -1});
      }      
      IsoCut_PerEfficiency_PerBin.insert({iEff, tempMap});
      
      TString NameEff = "Eff_";
      NameEff += iEff;
      TH3F* temp = new TH3F(NameEff,NameEff,nBinsIEta,0,nBinsIEta,nBinsIEt,0,nBinsIEt,nBinsnTT,0,nBinsnTT);
      IsoCut_PerBin.insert({iEff,temp});
      
      TString nameHisto1 = "pt_pass_efficiency_";
      nameHisto1 += iEff;
      
      TString nameHisto_TH3 = "pt_pass_efficiency_TH3_";
      nameHisto_TH3 += iEff;
      
      TH1F* temp_histo = new TH1F(nameHisto1,nameHisto1,100,0,200);
      TH1F* temp_histo_TH3 = new TH1F(nameHisto_TH3,nameHisto_TH3,100,0,200);
      
      pt_pass_efficiency.insert({iEff,temp_histo});
      pt_pass_efficiency_TH3.insert({iEff,temp_histo_TH3});
      
      TString nameHisto2 = "eta_pass_efficiency_";
      nameHisto2 += iEff;
      
      TH1F* temp_histo2 = new TH1F(nameHisto2,nameHisto2,100,0,100);
      eta_pass_efficiency.insert({iEff,temp_histo2});
      
      TString nameHisto3 = "nTT_pass_efficiency_";
      nameHisto3 += iEff;
      
      TH1F* temp_histo3 = new TH1F(nameHisto3,nameHisto3,180,20,200);
      nTT_pass_efficiency.insert({iEff,temp_histo3});
    }
    
    pt_large_eta = new  TH1F("Pt_large_eta",  "Pt_large_eta",  100, 0., 200.);
    eta_large_eta = new TH1F("Eta_large_eta", "Eta_large_eta", 100, 0., 100.);
    nTT_large_eta = new TH1F("NTT_large_eta", "NTT_large_eta", 180, 20., 200.);  
    iso_large_eta = new TH1F("Isolation_large_eta", "Isolation_large_eta", 100, 0., 100.);  
    
    pt_all = new  TH1F("Pt_all",  "Pt_all",  100, 0., 200.);
    eta_all = new TH1F("Eta_all", "Eta_all", 100, 0., 100.);
    nTT_all = new TH1F("NTT_all", "NTT_all", 180, 20., 200.);
    
  }
  if(option == "do_2" || option == "do_all") {
    TDirectory* td_step2 =   outputFile_->mkdir("Step2Histos");
    if (td_step2) td_step2->cd();
    int i=0;
    for (auto it : lutProgOptVec_) {
      std::string tmp_string = it;
      std::vector<std::string> options;
      tokenize(tmp_string,options,":");
      
      TString hname = "LUT_Progression_";
      hname += options[0];
      TH3F* th3 = new TH3F(hname,hname,
			   lutIEtaVec_.size()-1, 0, lutIEtaVec_.size()-1,
			   lutIEtVec_.size()-1, 0, lutIEtVec_.size()-1,
			   lutnTTVec_.size()-1, 0, lutnTTVec_.size()-1);
      lutProgHistoMap_.insert({it, th3});
    }  
    for(UInt_t iEff = 0 ; iEff <= 100 ; ++iEff)
      {
	TString NameHisto = "LUT_WP";
	NameHisto += iEff;
	TH3F* LUT_temp = new TH3F(NameHisto.Data(),NameHisto.Data(),
				  lutIEtaVec_.size()-1, 0, lutIEtaVec_.size()-1,
				  lutIEtVec_.size()-1, 0, lutIEtVec_.size()-1,
				  lutnTTVec_.size()-1, 0, lutnTTVec_.size()-1);
	LUT_WP.push_back(LUT_temp);
      }
  }
  bookedHistograms_ = true;
}
void IsolationAnalysis::fillLUTProgression(std::string option){ 
  Double_t AvEt;
  TFile*WPFile;
  if(option == "do_2") {
    WPFile = TFile::Open(outputWPFileName_.c_str(), "READ");
    if (!WPFile) return;
    std::cout<<"Filling LUT only"<<std::endl;
  }
  std::cout<<"Filling LUT"<<std::endl;
  Int_t maxEntriesForLUT=(lutIEtaVec_.size()-1)*(lutIEtVec_.size()-1)*(lutnTTVec_.size()-1);
  Int_t count=0;
  auto  t_start = std::chrono::high_resolution_clock::now();
  auto  t_end = std::chrono::high_resolution_clock::now();
  for(Int_t i = 0 ; i < lutIEtaVec_.size()-1; i++)
    {
      for(Int_t j = 0 ; j < lutIEtVec_.size()-1 ; j++)
	{
	  for(Int_t k = 0 ; k < lutnTTVec_.size()-1; k++)
	    {
          count++;
          if(count%200==0)
          {
             t_end = std::chrono::high_resolution_clock::now();
             std::cout<<"LUT Filling : Processing  loop id : "<<count<<" / "<<maxEntriesForLUT<<"  [ "<<100.0*count/maxEntriesForLUT<<"  % ]  "
                      << " Elapsed time : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()/1000.0
                      <<"  Estimated time left : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()*( maxEntriesForLUT - count)/(1e-9 + count)* 0.001
                      <<std::endl;

          }
	      for (auto it : lutProgHistoMap_)
		{
             
		  std::string tmp_string = it.first;
		  std::vector<std::string> options;
		  tokenize(tmp_string,options,":");
                  Double_t minPt =  std::stod(options[1]);
                  Double_t effLowMinPt = std::stod(options[2]);
		  Double_t reach100pc= std::stod(options[3]);
                  
		  Double_t Efficiency_Progression = findEfficiencyProgression((lutIEtVec_[j]+lutIEtVec_[j+1])/2.0, minPt,effLowMinPt,reach100pc);
		  if(Efficiency_Progression >= 0.9999) Efficiency_Progression = 1.0001;
		  Int_t Int_Efficiency_Progression = int(Efficiency_Progression*100);
		  TH3F* eff_histo;		  
		  if(option == "do_2") {
		    TString WPName = "Eff_";
		    WPName +=std::to_string(Int_Efficiency_Progression); 
		    eff_histo = dynamic_cast<TH3F*>(WPFile->Get("Step1Histos/"+WPName));
		    //std::cout<<WPName<<std::endl;		  
          }
		  else 
		    eff_histo = IsoCut_PerBin[Int_Efficiency_Progression];
          Int_t IsoCut_Progression = eff_histo->GetBinContent(i+1,j+1, k+1);
		  if(Int_Efficiency_Progression==100) IsoCut_Progression = 1000;
		  it.second->SetBinContent(i+1,j+1,k+1,IsoCut_Progression);
		}

	      for(UInt_t iEff = 0 ; iEff < 101 ; ++iEff)
		{
		  TH3F* th;
                  if(option == "do_2") {
		    TString WPName1 = "Eff_";
		    WPName1 +=TString::Itoa(iEff,10);
		    th = dynamic_cast<TH3F*>(WPFile->Get("Step1Histos/"+WPName1));
		  }
		  else
		    th= IsoCut_PerBin[iEff];
		  Int_t IsoCut = th->GetBinContent(i+1,j+1, k+1);

		  if(iEff==100) IsoCut = 1000;
		  LUT_WP.at(iEff)->SetBinContent(i+1,j+1,k+1,IsoCut);
		}		 
	      
	    }
	}
    }
}
Double_t IsolationAnalysis::findEfficiencyProgression(Double_t IEt, Double_t MinPt,
						      Double_t Efficiency_low_MinPt, Double_t Reaching_100pc_at) {
  Double_t Efficiency = 0; 
  Double_t Pt = IEt/2.;
  if(Pt>=Reaching_100pc_at) Efficiency = 1.;
  else if(Pt<MinPt) Efficiency = Efficiency_low_MinPt;
  else
    {
      Double_t  Slope = (1.-Efficiency_low_MinPt)/(Reaching_100pc_at-MinPt);
      Efficiency = Slope*Pt + (1. - Slope*Reaching_100pc_at);
    }
  
  if(Efficiency<0) Efficiency = 0.;
  if(Efficiency>=1) Efficiency = 1.;

  return Efficiency ;
}

void IsolationAnalysis::readLUTTable(std::string& file_name, unsigned int& nbin,
				     std::map<short, short>& lut_map){
  
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
      // enable '#' and '//' style comments
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
void IsolationAnalysis::saveHistograms() {
  if (outputFile_ && bookedHistograms_) {
    /*    for(UInt_t iEff = 0 ; iEff < 101 ; ++iEff) {
	  pt_pass_efficiency[iEff]->Delete();
	  pt_pass_efficiency_TH3[iEff]->Delete();
	  eta_pass_efficiency[iEff]->Delete();
	  nTT_pass_efficiency[iEff]->Delete();
	  
	  }*/
    /*    for (auto it : Histos_PerBin) {
      if (it.second) it.second->Delete();
    }*/
    
    outputFile_->Write();
    outputFile_->Close();
  }
}
void IsolationAnalysis::readParameters(const std::string jfile) {
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
      if(key=="NtupleFileName")        ntupleFileName_= value;
      else if (key=="OutputWPStepFileName") outputWPFileName_ = value.c_str();      
      else if (key=="OutputFileName")  outputFileName_ = value.c_str();
      else if (key=="EtLUTFileName")	readLUTTable(value,nBinsIEt, lutMapIEt_);
      else if (key=="EtaLUTFileName")	readLUTTable(value,nBinsIEta, lutMapIEta_);
      else if (key=="NTTLUTFileName")	readLUTTable(value,nBinsnTT, lutMapnTT_);
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
      else if (key=="MaxEntries")	    {
                maxEntries= atoi(value.c_str());
      }
      else 
	std::cout << " unknown option " << " key " << key << std::endl;
    }
  }
  jobcardFile.close();
  
}
void IsolationAnalysis::readTree()
{
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("l1tEmuNTT",&ntt);
  fChain->SetBranchAddress("l1tEmuRawEt",&et);
  fChain->SetBranchAddress("l1tEmuTowerIEta",&eta);
  fChain->SetBranchAddress("l1tEmuIsoEt",&iso);

  fChain->SetBranchAddress("eleProbeEta"  ,&eleProbeEta			);
  fChain->SetBranchAddress("eleProbePhi"  ,&eleProbePhi			);
  fChain->SetBranchAddress("eleTagEta"    ,&eleTagEta		    	);
  fChain->SetBranchAddress("eleTagPhi"    ,&eleTagPhi			    );
  fChain->SetBranchAddress("isProbeLoose" ,&isProbeLoose			);
}

bool IsolationAnalysis::getHistoBin(TString str, short& eta_bin, short& et_bin, short& ntt_bin){
  bool rval = false;
  if (str.Length() == 0) rval;
  
  str.Remove(0,5);
  TObjArray* toa = str.Tokenize("_");
  if (toa->GetEntries() >= 3) {
    eta_bin =  (dynamic_cast<TObjString*>(toa->At(0)))->GetString().Atoi();
    et_bin  =  (dynamic_cast<TObjString*>(toa->At(1)))->GetString().Atoi();
    ntt_bin =  (dynamic_cast<TObjString*>(toa->At(2)))->GetString().Atoi();
    rval = true;
  } else rval = false;
  return rval;
}
TString IsolationAnalysis::getHistoName(short eta, short et, short ntt){
  
  TString Name_Histo = "Hist_";
  
  unsigned int jEta, jEt, jnTT; 
  std::map<short, short >::iterator jEtaPos = lutMapIEta_.find(abs(eta)); 
  if (jEtaPos != lutMapIEta_.end()) jEta = jEtaPos->second;
  else                          jEta = nBinsIEta-1;
  
  std::map<short, short >::iterator jEtPos = lutMapIEt_.find(et); 
  if (jEtPos != lutMapIEt_.end()) jEt = jEtPos->second;
  else                          jEt = nBinsIEt-1;
  
  std::map<short, short >::iterator jnTTPos = lutMapnTT_.find(ntt); 
  if (jnTTPos != lutMapnTT_.end()) jnTT = jnTTPos->second;
  else                          jnTT = nBinsnTT-1;
  
  Name_Histo += jEta;
  Name_Histo += "_";
  Name_Histo += jEt;
  Name_Histo += "_";
  Name_Histo += jnTT;
  return Name_Histo;
}
 void IsolationAnalysis::tokenize(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters) {
  
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
  std::string to_do = argv[2];
  
  IsolationAnalysis treeReader(data_file);
  treeReader.readTree();
  
  treeReader.bookHistograms(to_do);
  
  if(to_do == "do_1" || to_do == "do_all") {
    std::cout << " Calling Analyse" << std::endl;
    treeReader.analyse();
  }
  if(to_do =="do_2" || to_do == "do_all") {
    std::cout << " Calling LUT Options" << std::endl;
    treeReader.fillLUTProgression(to_do);
  }  
  treeReader.saveHistograms();
  return 0;
}
