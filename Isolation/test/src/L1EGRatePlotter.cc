#include "L1EGRatePlotter.h"
#include "chrono"
#include <iostream>
#include <fstream>

L1EGRatePlotter::L1EGRatePlotter(std::string& inputFileName){
  readParameters(inputFileName);
  if (ntupleFileName_.size() == 0) {
    std::cout << " Inputfile list missing !!!" << ntupleFileName_ << std::endl;
    return;
  }
  accessTree(ntupleFileName_);

  assert(fChain);
  reportEvery_=10000;
  nttProfileTrue=nullptr;
  puProfileTrue=nullptr;
  nttProfile=nullptr;
  puProfile=nullptr;
  bookedHistograms_ = false;
  sumEventWeight=0.0;

}

void L1EGRatePlotter::accessTree(std::string & input_filelist) {
  std::ifstream myFile;
  myFile.open(input_filelist.c_str(), std::ios::in);
  if (!myFile) {
    std::cout << "Input File: " << input_filelist << " could not be opened!" << std::endl;
    fChain = 0;
    fChainForPU=0;
  } else {
    //    fChain = new TChain("l1UpgradeEmuTree/L1UpgradeTree");
    //fChain = new TChain("l1UpgradeEmuTree/L1UpgradeTree");
    //fChainForPU= new TChain("l1EventTree/L1EventTree");
    //fChain = new TChain("L1UpgradeEmuTree");
    fChain = new TChain("L1UpgradeTree");
    fChainForPU= new TChain("L1EventTree");
    static constexpr int BUF_SIZE = 256;
    char buf[BUF_SIZE];
    while (myFile.getline(buf, BUF_SIZE, '\n')) {  // Pops off the newline character
      std::string line(buf); 
      fChain->AddFile(line.c_str(),-1);
      fChainForPU->AddFile(line.c_str(),-1);
      std::cout << "Adding file " << line << " Entries " << fChain->GetEntries() << "( "<<fChainForPU->GetEntries()<<" )"<<  std::endl;
    }
    std::cout<<"Total : "<<fChain->GetEntries()<<"\n";
   // fChain->Print();
  }
}

void L1EGRatePlotter::loop()
{
   if (fChain == 0) return;
   nEntries_ = fChain->GetEntriesFast()*0 + 1e6;
   std::cout << " Entries " << nEntries_ << std::endl;
   Long64_t nbytes = 0, nb = 0;
   Float_t weight(0.0);
   Float_t puSum(0.0);
  auto t_start = std::chrono::high_resolution_clock::now();
  auto t_end = std::chrono::high_resolution_clock::now();
  Long64_t evtCount(0);

   for (Long64_t jentry=0; jentry < nEntries_; jentry++) {
      Long64_t ientry = fChain->LoadTree(jentry);
      if (ientry < 0) break;
      Long64_t iientry = fChainForPU->LoadTree(jentry);
      if (iientry < 0) break;
      
      
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      nb = fChainForPU->GetEntry(jentry);
      
      puProfileTrue->Fill(nPV_True);
      nttProfileTrue->Fill(nTT);
      
      weight=0.0;
//    if( nPV_True <= 50 and nPV_True>=44 ) { weight=1.0; evtCount++;}
      if( nPV_True <= 52 and nPV_True>=44 ) { weight=1.0; evtCount++;}
      else continue;
//     if (jentry > 20) break;
      int Et_Threshold=0; 
      int Bunches = 4;     

      number_EGsBranch->GetEntry(jentry);
      EG_ETBranch->GetEntry(jentry);
      EG_ETABranch->GetEntry(jentry);
      EG_BXBranch->GetEntry(jentry);
      EG_ISOBranch->GetEntry(jentry);
      
      //std::cout << " jentry " << jentry << " nEGs " << nEGs << " " << egEt.size() <<std::endl;

   // if(run != 352912) continue;
   // if( lumi < 100 or lumi > 124) continue;
    if(jentry%reportEvery_ == 0 )
      {
        	t_end = std::chrono::high_resolution_clock::now();
        	std::cout<<"Processing Entry in event loop (Rate) : "<<jentry<<" / "<<nEntries_<<"  [ "<<100.0*jentry/nEntries_<<"  % ]  "
        		 << " Elapsed time : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()/1000.0
        		 <<"  Estimated time left : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()*( nEntries_ - jentry)/(1e-9 + jentry)* 0.001
                 <<" evts  : "<<evtCount
        		 <<std::endl;
            std::cout<<"\t\t Run , Lumi : "<<run<<" , "<<lumi<<"\n";
      }

    nTT=0;
    if(nEGs >0) nTT=egNTT[0];

    for(Et_Threshold=0; Et_Threshold<=etMax_; Et_Threshold++)
	{    
	  Int_t nPassed_Single = 0;
	  Int_t nPassed_Double = 0;
	  Int_t nPassed_TightIso_Single = 0;	  	  
      Int_t nPassed_LooseIso_Single = 0;
	  Int_t nPassed_LooseIso_Double = 0;
	  bool flag_Double = false;
	  bool flag_LooseIso_Double = false;	  
	  for (UShort_t iEG=0; iEG < nEGs; ++iEG)
	    {
	      if (egBx[iEG]!=0)	continue;
	      float EG_Et  = egEt[iEG];
	      float EG_Eta  = egEta[iEG];
	      short EG_Iso = egIso[iEG];

          if( abs(EG_Eta) > 2.5 ) continue;
          //std::cout<<jentry<<" , "<<iEG<<" Eta : "<<EG_Eta<<"\n";

	      if ( EG_Et >= Et_Threshold) {
		nPassed_Single++;

		if (EG_Et >= Et_Threshold+10) flag_Double = true;
		nPassed_Double++;
		std::cout<<EG_Iso<<"\n";
		if (EG_Iso==1 || EG_Iso==3) {
            nPassed_TightIso_Single++;
                std::cout<<" isTight \n";
		}
        if (EG_Iso==2 || EG_Iso==3) {
                std::cout<<" isLoose \n";
		  nPassed_LooseIso_Single++;
		  if (EG_Et >= Et_Threshold+10) flag_LooseIso_Double = true;
		  nPassed_LooseIso_Double++;	      
		}
	      }
	      if (EG_Et >= Et_Threshold) {
	      }		
	    }
	  for (std::map<std::string, TH1F*>::iterator it =  histoMap_.begin(); it != histoMap_.end(); it++)
	    {
	      bool flag = false;
	      
	      if      ( (it->first == "SingleEG_rate_inclusive") && nPassed_Single >= 1) flag = true;
	      else if ( (it->first == "DoubleEG_rate_inclusive") && nPassed_Double >= 2 &&  flag_Double) flag = true;
	      else if ( (it->first == "SingleEG_rate_TightIso" ) && nPassed_TightIso_Single >= 1) flag = true;
          else if ( (it->first == "SingleEG_rate_LooseIso" ) && nPassed_LooseIso_Single >= 1) flag = true;
          else if ( (it->first == "DoubleEG_rate_LooseIso" ) && nPassed_Double >= 2 && flag_LooseIso_Double) flag = true;
	      else flag = false;
	      if (flag) it->second->Fill(Et_Threshold);
	    }
	}
     sumEventWeight+=weight;
     puProfile->Fill(nPV_True,weight);
     nttProfile->Fill(nTT,weight);
     puSum+=weight*nPV_True;
   }

   std::cout<<" Average PU = "<<puSum/sumEventWeight<<"\n";
}

void L1EGRatePlotter::bookHistograms() 
{

  outputFile_ = new TFile(outputFileName_.c_str(), "RECREATE");
  outputFile_->cd();

  for (auto it : histoLabelVec_) {
    std::string hname = it;
    TH1F* th = new TH1F(hname.c_str() ,hname.c_str(), etMax_, 0., etMax_*1.0);
    th->GetYaxis()->SetTitle("Rate [kHz]");
    th->GetXaxis()->SetTitle("E_{T}^{L1} Threshold [GeV]");
    th->GetXaxis()->SetTitleOffset(1.3);
    th->Sumw2();
    std::cout << " Booking Histogram " << hname << std::endl;
    histoMap_.insert(std::make_pair(hname, th));
    
  }
 
  puProfileTrue  = new TH1F("puProfileTrue" ,"PU Profile for the Dataset ", 100, -0.5,  99.5);
  nttProfileTrue = new TH1F("nttProfileTrue","nTT Profile for the Dataset ", 200, -0.5, 199.5);
  puProfile  = new TH1F("puProfile" ,"PU Profile for the Dataset ", 100, -0.5,  99.5);
  nttProfile = new TH1F("nttProfile","nTT Profile for the Dataset ", 200, -0.5, 199.5);
  
  bookedHistograms_ = true;
}
void L1EGRatePlotter::scaleHistograms(){
  int Bunches = 2544;     
  double scaleFactor = 11.245*Bunches;
  if (sumEventWeight>0) scaleFactor /= sumEventWeight;
  else scaleFactor = 1.0;
  int xbin;
  float vl,bl;
  float frate=3.9;
  for (std::map<std::string, TH1F*>::iterator it =  histoMap_.begin(); it != histoMap_.end(); it++) {
    TH1F* th = it->second;
    for (Int_t i = 1; i < th->GetNbinsX()+1; ++i) {
      Double_t cont = th->GetBinContent(i);
      Double_t err = th->GetBinError(i);
      th->SetBinContent(i, cont*scaleFactor);
      th->SetBinError(i, err*scaleFactor);
    }

   if      ( (it->first == "SingleEG_rate_inclusive") )  frate=16.7;
   else if ( (it->first == "DoubleEG_rate_inclusive") )  frate=3.9 ;
   else if ( (it->first == "SingleEG_rate_TightIso" ) )  frate=16.7;
   else if ( (it->first == "SingleEG_rate_LooseIso" ) )  frate=16.7;
   else if ( (it->first == "DoubleEG_rate_LooseIso" ) )  frate=3.9;
   else frate=16.7;

  for(UInt_t i=0;i< etMax_;i++) {
    float rate = th->GetBinContent(i+1);
    if(rate > frate) {
      bl =i;
      vl=rate;
    }
    
    if(rate <=frate){
      xbin = (rate - frate) < abs(rate-vl) ? i : bl;
      break;
    }
  }
   std::cout<< it->first<<" , "<<"et: "<<xbin<<" @ "<<frate<<" +/- "<<th->GetBinError(xbin) <<" kHz "<<std::endl;
  
  }
}

void L1EGRatePlotter::saveHistograms() {
  if (outputFile_ && bookedHistograms_) {
    scaleHistograms();
    puProfileTrue->Write();
    nttProfileTrue->Write();
    puProfile->Write();
    nttProfile->Write();
    outputFile_->Write();
    outputFile_->Close();
  }
}
void L1EGRatePlotter::readParameters(const std::string jfile) {
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
      std::cout<<line<<"\n";
      tokenize(line,tokens,"=");
      std::cout << tokens[0] << ":" << tokens[1] << std::endl;
      std::string key   = tokens.at(0);
      std::string value = tokens.at(1); 
      if(key=="NtupleFileName")        ntupleFileName_= value;
      else if (key=="ReportEvery")     reportEvery_ = std::stoi(value.c_str());
      else if (key=="OutputFileName")  outputFileName_ = value.c_str();
      else if (key=="EtMaxCutoff")     etMax_ = std::atoi(value.c_str());      
      else if (key=="RateHistoNames")      
	{
	  std::string tmp_string = value;
	  tokenize(tmp_string,histoLabelVec_,",");
	}   
      
      else 
	std::cout << " unknown option " << " key " << key << std::endl;
    }
  }
  jobcardFile.close();

}
L1EGRatePlotter::~L1EGRatePlotter()
{
  //  if (!fChain) return;
  //  delete fChain->GetCurrentFile();
}

void L1EGRatePlotter::readTree()
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).
  fChain->SetMakeClass(1);
  
  fChain->SetBranchAddress("nEGs",  &nEGs, &number_EGsBranch);
  fChain->SetBranchAddress("egEt",  &egEt, &EG_ETBranch);
  fChain->SetBranchAddress("egEta", &egEta, &EG_ETABranch);
  fChain->SetBranchAddress("egPhi", &egPhi, &EG_PHIBranch);
  fChain->SetBranchAddress("egIso", &egIso, &EG_ISOBranch);
  fChain->SetBranchAddress("egBx",  &egBx, &EG_BXBranch);
  fChain->SetBranchAddress("egNTT", &egNTT);
  
  fChainForPU->SetMakeClass(1);
  fChainForPU->SetBranchAddress("nPV_True", &nPV_True);
  fChainForPU->SetBranchAddress("run", &run);
  fChainForPU->SetBranchAddress("lumi", &lumi);
}

void L1EGRatePlotter::tokenize(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters) {

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
  L1EGRatePlotter treeReader(data_file);
  treeReader.readTree();
  treeReader.bookHistograms();
  std::cout << " Calling Loop" << std::endl;
  treeReader.loop();
  treeReader.saveHistograms();
  return 0;
}
