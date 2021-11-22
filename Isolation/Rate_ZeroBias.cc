#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <sstream>
#include <string>
#include <ctime>
#include <cstdlib>

#include <TROOT.h>
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLine.h"
#include "TString.h"
#include "TMath.h"
#include "TString.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TTree.h"
#include "TFile.h"
#include "TLatex.h"
#include "TMarker.h"
#include "TAxis.h"
#include "L1UpgradeTree.h"

#define N_OPTIONS 23 // 1 [ no Iso ] + 22
#define ET_MAX 255

//using namespace std;
TCanvas* CreateCanvas(TString CanvasName = "myPlot", bool LogY = false, bool Grid = true);
void DrawPrelimLabel(TCanvas* c);
void DrawLumiLabel(TCanvas* c, TString Lumi = "35.9");
void SaveCanvas(TCanvas* c, TString PlotName = "myPlotName");
void tokenize(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters);
void readLUTTable(std::string& file_name, unsigned int& nbin, std::map<unsigned int, unsigned int>& lut_map);

void Rate_ZeroBias(TString optionFile="",TString zeroBiasFile="", TString outFile="", float Bunches=1.0)
{
   // "calibFiles/data2018_EGIsoV2_Options.root"
   TFile f_Isolation(optionFile,"READ");
  

   std::map<TString,TH3F*> histosIsolation;
   TChain *treeChain = new TChain("L1UpgradeTree");
   
   //"/grid_mnt/t3storage3/athachay/l1egamma/isolation/CMSSW_7_6_0/src/Isolation/slimL1EmulatorTrees/L1Ntuple_ZeroBias_325170_reEmulation_11_2_0_run3Recalib.root";
   treeChain->Add(zeroBiasFile);

   cout << "Total Number of Events Available : " << treeChain->GetEntriesFast()  << endl;
   L1UpgradeTree ntupleRawTree(treeChain);

 
  std::map<unsigned int, unsigned int> lutMapEt;
  std::map<unsigned int, unsigned int> lutMapEta;
  std::map<unsigned int, unsigned int> lutMapNTT;
  
  for(UInt_t i = 1 ; i < N_OPTIONS ; ++i)
    {
      if(i == 0) continue ;
      TString CurrentNameHisto = "LUT_Progression_";
      ostringstream convert;
      convert << i;
      CurrentNameHisto += convert.str();
      std::cout<<"Reading Option  :  "<<CurrentNameHisto<<"\n";
      TH3F* current_Histo = (TH3F*)f_Isolation.Get(CurrentNameHisto.Data());
      histosIsolation.insert(make_pair(CurrentNameHisto,current_Histo));
    } 

   int optionsPtIsoMatrix[100][ET_MAX+1];
 
   // Integrating the event pass distribution
    for(int i=0;i< N_OPTIONS ; i++)
    {
        for(int j= 0 ; j < ET_MAX  ;j++)
        {
            optionsPtIsoMatrix[i][j]= 0 ;
       }
    }
   UShort_t nEGs ;
   std::vector<short>   egIEt;
   std::vector<short>   egIEta;
   std::vector<short>   egIsoEt;
   std::vector<short>   egNTT;

   ntupleRawTree.fChain->SetBranchStatus("*",0);
   ntupleRawTree.fChain->SetBranchStatus("nEGs",1);
   ntupleRawTree.fChain->SetBranchStatus("egIEt",1);
   ntupleRawTree.fChain->SetBranchStatus("egIEta",1);
   ntupleRawTree.fChain->SetBranchStatus("egIsoEt",1);
   ntupleRawTree.fChain->SetBranchStatus("egNTT",1);

  
   short in_compressedieta;
   short in_compressedE;
   short in_compressednTT;
   short clippedEta;
   short clippedEt;

   int eTMax[N_OPTIONS];

   //inTree->Print();
   //inTree->Scan("nEGs");
   
   Long64_t nentries = ntupleRawTree.fChain->GetEntries();
   Long64_t nbytes = 0, nb = 0;
  
    auto t_start = std::chrono::high_resolution_clock::now();
    auto t_end = std::chrono::high_resolution_clock::now();

   
 cout << "Total Number of Events Available : " <<ntupleRawTree.fChain->GetEntries()  << endl;
 cout << "Total Number of Events to process : " << nentries  << endl;
 Long64_t nEventsConsidered=0;
 for (Long64_t jentry=0; jentry<nentries;jentry++) {

      Long64_t ientry = ntupleRawTree.LoadTree(jentry);
      if (ientry < 0) break;

         nb = ntupleRawTree.fChain->GetEntry(jentry);   nbytes += nb;
         
         nEventsConsidered++;

        if(jentry%12000 == 0) 
        {
             t_end = std::chrono::high_resolution_clock::now();
             std::cout<<"Processing Entry in event loop : "<<jentry<<" / "<<nentries<<"  [ "<<100.0*jentry/nentries<<"  % ]  "
             << " Elapsed time : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()/1000.0
             <<"  Estimated time left : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()*( nentries - jentry)/jentry* 0.001
             <<std::endl;
        }

         //std::cout<<"Entry #"<<jentry<<endl; 
         //std::cout<<"\tnEGs = "<<ntupleRawTree.nEGs<<"\n";


         // if(PU_per_LS.find(in_lumi)==PU_per_LS.end()) continue;
         // Float_t weight = PU_per_LS[48]/PU_per_LS[in_lumi];
         // if(weight<0.5) cout<<"lumi = "<<in_lumi<<", weight = "<<weight<<endl;
         // if(weight>2) cout<<"lumi = "<<in_lumi<<", weight = "<<weight<<endl;

         Float_t weight = 1.;
         for(UInt_t i = 0 ; i < N_OPTIONS ; ++i)
         {
               eTMax[i] = -1; 
         }
         for(auto  eg_idx=0 ; eg_idx < ntupleRawTree.nEGs ; eg_idx++)
         {
             clippedEt  = ntupleRawTree.egIEt.at(eg_idx) < ET_MAX ? ntupleRawTree.egIEt.at(eg_idx) : ET_MAX ;
             //std::cout<<"\t eg_idx "<<eg_idx<<" = "<<ntupleRawTree.egIEt.at(eg_idx)<<" [ "<<clippedEt<<" ] isoEt : "<<ntupleRawTree.egIsoEt.at(eg_idx)<<" \n";
             
             if( ntupleRawTree.egIEta.at(eg_idx) > 31 )
                clippedEta = 31 ;
             else if ( ntupleRawTree.egIEta.at(eg_idx) < -31 )
                clippedEta= -31 ;
             else
             {
                clippedEta = ntupleRawTree.egIEta.at(eg_idx);
             }

             in_compressedieta =  lutMapEt[clippedEta]; 
             in_compressedE    =  lutMapEta[clippedEt];
             in_compressednTT  =  lutMapNTT[ntupleRawTree.egNTT.at(eg_idx)];

             if(clippedEt > eTMax[0]) eTMax[0] = clippedEt;
            

             for(UInt_t i = 1 ; i < N_OPTIONS ; ++i)
               {
                 TString CurrentNameHisto = "LUT_Progression_";
                 ostringstream convert;
                 convert << i;
                 CurrentNameHisto += convert.str();
                 Int_t Cut_L1EG_Iso   = histosIsolation[CurrentNameHisto] ->GetBinContent(in_compressedieta+1,in_compressedE+1,in_compressednTT+1);
                 
                 //std::cout<<" option : "<<i<<" : cut : "<<Cut_L1EG_Iso<<"\n";

                 if ( ntupleRawTree.egIsoEt.at(eg_idx) > Cut_L1EG_Iso ) continue ;
                 if(clippedEt > eTMax[i]) eTMax[i] = clippedEt;
               } 
          }
          
          for(UInt_t i=0 ;i < N_OPTIONS ; i++)
          {
               if( eTMax[i] > 0  )
               {
                   //std::cout<<" Pass for  "<<i<<" for "<<eTMax[i]<<"\n";
                   optionsPtIsoMatrix[i][eTMax[i]] +=1;
               }
          }

    }

    // Integrating the event pass distribution
    for(int i=0;i< N_OPTIONS ; i++)
    {
        for(int j= ET_MAX  ; j > 0  ;j--)
        {
            optionsPtIsoMatrix[i][j-1]+=optionsPtIsoMatrix[i][j] ;
        }
    }

  // Filling the hist
  TFile *file =new TFile(outFile,"RECREATE");
  file->cd();

  double scaleFactor = 11.2456 * Bunches /nEventsConsidered  ;
  // No Isolation
      TString CurrentNameHisto = "pass_noIsolation";
      TH1F* pass_Histo = new TH1F(CurrentNameHisto, CurrentNameHisto , ET_MAX , 0.0 - 0.5 , ET_MAX -0.5 );
      CurrentNameHisto = "rate_noIsolation";
      TH1F* rate_Histo = new TH1F(CurrentNameHisto, CurrentNameHisto ,ET_MAX , 0.0 - 0.5 , ET_MAX -0.5 );
      
      for(int j=1;j<=ET_MAX;j++)
        {
            pass_Histo->SetBinContent(j,optionsPtIsoMatrix[0][j-1])  ;
            rate_Histo->SetBinContent(j,optionsPtIsoMatrix[0][j-1]*scaleFactor) ;
        }
  pass_Histo->Write();
  rate_Histo->Write();

  for(UInt_t i = 1 ; i < N_OPTIONS ; ++i )
    {
      CurrentNameHisto = "pass_LUT_Progression_";
      ostringstream convert;
      convert.clear();  convert << i;
      CurrentNameHisto += convert.str();
      pass_Histo = new TH1F(CurrentNameHisto, CurrentNameHisto , ET_MAX , 0.0 - 0.5 , ET_MAX -0.5 );
      CurrentNameHisto = "rate_LUT_Progression_" + convert.str();
      rate_Histo = new TH1F(CurrentNameHisto, CurrentNameHisto ,ET_MAX , 0.0 - 0.5 , ET_MAX -0.5 );

      for(int j=1;j<=ET_MAX;j++)
        {
            pass_Histo->SetBinContent(j,optionsPtIsoMatrix[i][j])  ;
            rate_Histo->SetBinContent(j,optionsPtIsoMatrix[i][j]*scaleFactor) ;
        }

      pass_Histo->Write();
      rate_Histo->Write();
    } 

    file->Write();
    file->Close();

}

void tokenize(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters) {

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

void readLUTTable(std::string& file_name, unsigned int& nbin,
				     std::map<unsigned int, unsigned int>& lut_map){

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



