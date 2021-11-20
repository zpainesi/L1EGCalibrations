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


#define N_OPTIONS 23 // 1 [ no Iso ] + 22
#define ET_MAX 255

//using namespace std;
TCanvas* CreateCanvas(TString CanvasName = "myPlot", bool LogY = false, bool Grid = true);
void DrawPrelimLabel(TCanvas* c);
void DrawLumiLabel(TCanvas* c, TString Lumi = "35.9");
void SaveCanvas(TCanvas* c, TString PlotName = "myPlotName");

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


void rate_calculation()
{


    std::map<TString,TH3F*> histosIsolation;
   TFile f_Isolation("/home/llr/cms/mperez/EGTrigger/CMSSW_9_2_10/src/EGTagAndProbe/EGTagAndProbe/test/Trees_meanparam/Iso_LUTs_Options_MC_92X_mean_newnTT.root","READ");
   TString FileName_in = "/home/llr/cms/mperez/EGTrigger/CMSSW_9_2_10/src/EGTagAndProbe/EGTagAndProbe/test/Trees_meanparam/calibratedOutput_Ntuple_ZeroBias_Run305310_92X_mean_newnTT.root";
  
  std::map<unsigned int, unsigned int> lutMapEt;
  std::map<unsigned int, unsigned int> lutMapEta;
  std::map<unsigned int, unsigned int> lutMapNTT;
  
  for(UInt_t i = 0 ; i < N_OPTIONS ; ++i)
    {
      if(i == 0) continue ;

      TString CurrentNameHisto = "LUT_Progression_";
      ostringstream convert;
      convert << i;
      CurrentNameHisto += convert.str();
      TH3F* current_Histo = (TH3F*)f_Isolation.Get(CurrentNameHisto.Data());
      histosIsolation.insert(make_pair(CurrentNameHisto,current_Histo));
    } 


   int optionsPtIsoMatrix[100][ET_MAX+1];

   TFile f_in(FileName_in.Data(),"READ");
   TTree* inTree = (TTree*)f_in.Get("outTreeForCalibration");
 
   UShort_t nEGs ;
   std::vector<short>   egIEt;
   std::vector<short>   egIEta;
   std::vector<short>   egIsoEt;
   std::vector<short>   egNTT;
  

   inTree->SetBranchAddress("nEGs"   , &nEGs    );
   inTree->SetBranchAddress("egIEt"  , &egIEt   );
   inTree->SetBranchAddress("egIEta" , &egIEta  );
   inTree->SetBranchAddress("egNTT"  , &egNTT   );
   inTree->SetBranchAddress("isoEt"  , &egIsoEt );
   
   int eTMin[N_OPTIONS];

   for(UInt_t i = 0 ; i < inTree->GetEntries() ; ++i)
    {
         inTree->GetEntry(i);
         if(i%10000==0) cout<<"Entry #"<<i<<endl; 

         // if(PU_per_LS.find(in_lumi)==PU_per_LS.end()) continue;
         // Float_t weight = PU_per_LS[48]/PU_per_LS[in_lumi];
         // if(weight<0.5) cout<<"lumi = "<<in_lumi<<", weight = "<<weight<<endl;
         // if(weight>2) cout<<"lumi = "<<in_lumi<<", weight = "<<weight<<endl;

         Float_t weight = 1.;

         ++Denominator;
         for(UInt_t i = 0 ; i < N_OPTIONS ; ++i)
         {
               eTMin = 10 +  ET_MAX; 
         }

         for(auto  eg_idx=0 ; eg_idx < nEGs ; eg_idx++)
         {
             clippedEt  = L1_EG_Et[eg_idx] < ET_MAX ? L1_EG_Et[eg_idx] : ET_MAX ;

             if( egEta[eg_idx] > 31 )
                clippedEta = 31 ;
             else if ( egEta[eg_idx] < -31 )
                clippedEta= -31 ;
             else
             {
                clippedEta = egEta[eg_idx];
             }

             in_compressedieta =  lutMapEt[clippedEta]; 
             in_compressedE    =  lutMapEta[clippedEt];
             in_compressednTT  =  lutMapNTT[egNTT[eg_idx]];

             for(UInt_t i = 1 ; i < N_OPTIONS ; ++i)
               {
                 TString CurrentNameHisto = "LUT_Progression_";
                 ostringstream convert;
                 convert << i;
                 CurrentNameHisto += convert.str();
                 Int_t Cut_L1EG_Iso   = histosIsolation[CurrentNameHisto] ->GetBinContent(in_compressedieta+1,in_compressedE+1,in_compressednTT+1);

                 if ( egIsoEt[eg_idx] > Cut_L1EG_Iso ) continue ;
                 if(clippedEt < eTMin[i]) eTMin[i] = clippedEt;
               } 
             
          }
          
          for(UInt_t i=0 ;i < N_OPTIONS ; i++)
          {
               if( eTMin[i] <= ET_MAX )
               {
                   optionsPtIsoMatrix[i][eTMin[i]] +=1;
               }
          }

    }

    // Integrating the event pass distribution
    for(int i=0;i< N_OPTIONS ; i++)
    {
        for(int j= ET_MAX-1 ; j >=0 ;j--)
        {
            optionsPtIsoMatrix[i][j]+=optionsPtIsoMatrix[i][j+1] ;
        }
    }

  // Filling the hist
  ostringstream convert;
  TFile *file =new TFile("L1EG_Rate_2018.root","RECREATE");
  file->cd();

  for(UInt_t i = 1 ; i < 23 ; ++i)
    {
      TString CurrentNameHisto = "pass_LUT_Progression_";
      convert.clea();  convert << i;
      CurrentNameHisto += convert.str();
      TH3F* pass_Histo = new TH3F(CurrentNameHisto.Data(), ET_MAX , 0.0 - 0.5 , ET_MAX -0.5 );
      CurrentNameHisto = "rate_LUT_Progression_" + convert.c_str();
      TH3F* rate_Histo = new TH3F(CurrentNameHisto.Data(), ET_MAX , 0.0 - 0.5 , ET_MAX -0.5 );

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



TCanvas* CreateCanvas(TString CanvasName = "myPlot", bool LogY = false, bool Grid = true)
{
    TCanvas* c = new TCanvas(CanvasName.Data(),CanvasName.Data(),800,800);
    c->SetLeftMargin(0.11);
    if(Grid)
        c->SetGrid();
    if(LogY)
        c->SetLogy();
    return c;
}

void DrawPrelimLabel(TCanvas* c)
{
    c->cd();

    TLatex tex;
    tex.SetTextSize(0.03);
    tex.DrawLatexNDC(0.11,0.91,"#scale[1.5]{CMS}");
    tex.Draw("same");

    return;
}

void DrawLumiLabel(TCanvas* c, TString Lumi = "35.9")
{
    c->cd();

    TLatex tex;
    tex.SetTextSize(0.035);
    TString toDisplay = Lumi + " pb^{-1} (13 TeV)";
    tex.DrawLatexNDC(0.66,0.91,toDisplay.Data());
    tex.Draw("same");

    return;
}

void SaveCanvas(TCanvas* c, TString PlotName = "myPlotName")
{
    c->cd();
    c->SaveAs(PlotName + ".pdf");
    c->SaveAs(PlotName + ".root");

    return;
}


