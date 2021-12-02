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
#include "ApplyCalibration.C"

using namespace std;

Double_t FindEfficiency_Progression(Double_t IEt, Double_t MinPt, Double_t Efficiency_low_MinPt, Double_t Reaching_100pc_at)
{
  Double_t Efficiency = 0; 
  Double_t Pt = IEt/2.;

  if(Pt>=Reaching_100pc_at) Efficiency = 1.;
  else if(Pt<MinPt) Efficiency = Efficiency_low_MinPt;
  else
    {
      Double_t Slope = (1.-Efficiency_low_MinPt)/(Reaching_100pc_at-MinPt);
      Efficiency = Slope*Pt + (1. - Slope*Reaching_100pc_at);
      // Efficiency = (Effiency_low_MinPt-(1.-Effiency_low_MinPt)) + (1.-Effiency_low_MinPt)/(Reaching_100pc_at-MinPt)*Pt;
    }

  if(Efficiency<0) Efficiency = 0.;
  if(Efficiency>=1) Efficiency = 1.;

  return Efficiency ;
}


void Fill_Isolation(TString inFile,TString oFileName="Iso_LUTs_Options.root")
{
  std::map<TString,TH3F*> histosIsolation;
  TFile  f_Isolation( inFile     ,  "READ"    );
  TFile LUTs_Options( oFileName  ,  "RECREATE");

  // TFile f_Isolation("Iso_LUTs_Distributions_ggH_DY_NewCompression_NewLayer1_Thomas_nTT_OlivierFlatWP_SK1616.root");
  // TFile f_Isolation("isolationLUTS/Iso_LUTs_Distributions_ggH_DY_NewCompression_Layer1_Thomas_nTT.root");
  // TFile f_Isolation("isolationLUTS/Iso_LUTs_Distributions_ggH_DY_NewCompression_Layer1.root");
  // TFile f_Isolation("isolationLUTS/Iso_LUTs_Distributions_DownTo0pcEff.root");

  for(UInt_t i = 0 ; i < 101 ; ++i)
    {
      TString CurrentNameHisto = "Eff_";
      ostringstream convert;
      convert << i;
      CurrentNameHisto += TString(convert.str());
      TH3F* current_Histo = (TH3F*)f_Isolation.Get(CurrentNameHisto.Data());
      histosIsolation.insert(make_pair(TString(convert.str()),current_Histo));
    }  
  TH3F* LUT_Progression_1  = new TH3F("LUT_Progression_1","LUT_Progression_1",NbinsIEta-1,0,NbinsIEta-1,NbinsIEt2-1,0,NbinsIEt2-1,NbinsnTT2-1,0,NbinsnTT2-1);
  TH3F* LUT_Progression_2  = new TH3F("LUT_Progression_2","LUT_Progression_2",NbinsIEta-1,0,NbinsIEta-1,NbinsIEt2-1,0,NbinsIEt2-1,NbinsnTT2-1,0,NbinsnTT2-1);
  TH3F* LUT_Progression_3  = new TH3F("LUT_Progression_3","LUT_Progression_3",NbinsIEta-1,0,NbinsIEta-1,NbinsIEt2-1,0,NbinsIEt2-1,NbinsnTT2-1,0,NbinsnTT2-1);
  TH3F* LUT_Progression_4  = new TH3F("LUT_Progression_4","LUT_Progression_4",NbinsIEta-1,0,NbinsIEta-1,NbinsIEt2-1,0,NbinsIEt2-1,NbinsnTT2-1,0,NbinsnTT2-1);
  TH3F* LUT_Progression_5  = new TH3F("LUT_Progression_5","LUT_Progression_5",NbinsIEta-1,0,NbinsIEta-1,NbinsIEt2-1,0,NbinsIEt2-1,NbinsnTT2-1,0,NbinsnTT2-1);
  TH3F* LUT_Progression_6  = new TH3F("LUT_Progression_6","LUT_Progression_6",NbinsIEta-1,0,NbinsIEta-1,NbinsIEt2-1,0,NbinsIEt2-1,NbinsnTT2-1,0,NbinsnTT2-1);
  TH3F* LUT_Progression_7  = new TH3F("LUT_Progression_7","LUT_Progression_7",NbinsIEta-1,0,NbinsIEta-1,NbinsIEt2-1,0,NbinsIEt2-1,NbinsnTT2-1,0,NbinsnTT2-1);
  TH3F* LUT_Progression_8  = new TH3F("LUT_Progression_8","LUT_Progression_8",NbinsIEta-1,0,NbinsIEta-1,NbinsIEt2-1,0,NbinsIEt2-1,NbinsnTT2-1,0,NbinsnTT2-1);
  TH3F* LUT_Progression_9  = new TH3F("LUT_Progression_9","LUT_Progression_9",NbinsIEta-1,0,NbinsIEta-1,NbinsIEt2-1,0,NbinsIEt2-1,NbinsnTT2-1,0,NbinsnTT2-1);
  TH3F* LUT_Progression_10 = new TH3F("LUT_Progression_10","LUT_Progression_10",NbinsIEta-1,0,NbinsIEta-1,NbinsIEt2-1,0,NbinsIEt2-1,NbinsnTT2-1,0,NbinsnTT2-1);

  TH3F* LUT_Progression_11 = new TH3F("LUT_Progression_11","LUT_Progression_11",NbinsIEta-1,0,NbinsIEta-1,NbinsIEt2-1,0,NbinsIEt2-1,NbinsnTT2-1,0,NbinsnTT2-1);
  TH3F* LUT_Progression_12 = new TH3F("LUT_Progression_12","LUT_Progression_12",NbinsIEta-1,0,NbinsIEta-1,NbinsIEt2-1,0,NbinsIEt2-1,NbinsnTT2-1,0,NbinsnTT2-1);
  TH3F* LUT_Progression_13 = new TH3F("LUT_Progression_13","LUT_Progression_13",NbinsIEta-1,0,NbinsIEta-1,NbinsIEt2-1,0,NbinsIEt2-1,NbinsnTT2-1,0,NbinsnTT2-1);
  TH3F* LUT_Progression_14 = new TH3F("LUT_Progression_14","LUT_Progression_14",NbinsIEta-1,0,NbinsIEta-1,NbinsIEt2-1,0,NbinsIEt2-1,NbinsnTT2-1,0,NbinsnTT2-1);
  TH3F* LUT_Progression_15 = new TH3F("LUT_Progression_15","LUT_Progression_15",NbinsIEta-1,0,NbinsIEta-1,NbinsIEt2-1,0,NbinsIEt2-1,NbinsnTT2-1,0,NbinsnTT2-1);
  TH3F* LUT_Progression_16 = new TH3F("LUT_Progression_16","LUT_Progression_16",NbinsIEta-1,0,NbinsIEta-1,NbinsIEt2-1,0,NbinsIEt2-1,NbinsnTT2-1,0,NbinsnTT2-1);
  TH3F* LUT_Progression_17 = new TH3F("LUT_Progression_17","LUT_Progression_17",NbinsIEta-1,0,NbinsIEta-1,NbinsIEt2-1,0,NbinsIEt2-1,NbinsnTT2-1,0,NbinsnTT2-1);
  TH3F* LUT_Progression_18 = new TH3F("LUT_Progression_18","LUT_Progression_18",NbinsIEta-1,0,NbinsIEta-1,NbinsIEt2-1,0,NbinsIEt2-1,NbinsnTT2-1,0,NbinsnTT2-1);
  TH3F* LUT_Progression_19 = new TH3F("LUT_Progression_19","LUT_Progression_19",NbinsIEta-1,0,NbinsIEta-1,NbinsIEt2-1,0,NbinsIEt2-1,NbinsnTT2-1,0,NbinsnTT2-1);
  TH3F* LUT_Progression_20 = new TH3F("LUT_Progression_20","LUT_Progression_20",NbinsIEta-1,0,NbinsIEta-1,NbinsIEt2-1,0,NbinsIEt2-1,NbinsnTT2-1,0,NbinsnTT2-1);
  TH3F* LUT_Progression_21 = new TH3F("LUT_Progression_21","LUT_Progression_21",NbinsIEta-1,0,NbinsIEta-1,NbinsIEt2-1,0,NbinsIEt2-1,NbinsnTT2-1,0,NbinsnTT2-1);
  TH3F* LUT_Progression_22 = new TH3F("LUT_Progression_22","LUT_Progression_22",NbinsIEta-1,0,NbinsIEta-1,NbinsIEt2-1,0,NbinsIEt2-1,NbinsnTT2-1,0,NbinsnTT2-1);

  std::vector<TH3F*> LUT_WP ;
  for(UInt_t iEff = 0 ; iEff <= 100 ; ++iEff)
    {
      stringstream ss_i;
      ss_i << iEff;
      TString Appendix_i = TString(ss_i.str());

      TString NameHisto = "LUT_WP";
      NameHisto += Appendix_i ;
      TH3F* LUT_temp = new TH3F(NameHisto.Data(),NameHisto.Data(),NbinsIEta-1,0,NbinsIEta-1,NbinsIEt2-1,0,NbinsIEt2-1,NbinsnTT2-1,0,NbinsnTT2-1);
      LUT_WP.push_back(LUT_temp);
    }

  // TH3F* LUT_WP70 = new TH3F("LUT_WP70","LUT_WP70",NbinsIEta-1,0,NbinsIEta-1,NbinsIEt2-1,0,NbinsIEt2-1,NbinsnTT2-1,0,NbinsnTT2-1);
  // TH3F* LUT_WP80 = new TH3F("LUT_WP80","LUT_WP80",NbinsIEta-1,0,NbinsIEta-1,NbinsIEt2-1,0,NbinsIEt2-1,NbinsnTT2-1,0,NbinsnTT2-1);
  // TH3F* LUT_WP90 = new TH3F("LUT_WP90","LUT_WP90",NbinsIEta-1,0,NbinsIEta-1,NbinsIEt2-1,0,NbinsIEt2-1,NbinsnTT2-1,0,NbinsnTT2-1);
  int progressionCount=0;
  for(Int_t i = 0 ; i < NbinsIEta-1 ; ++i)
    {
      for(Int_t j = 0 ; j < NbinsIEt2-1 ; ++j)
	{
	  for(Int_t k = 0 ; k < NbinsnTT2-1 ; ++k)
	  // for(Int_t k = 0 ; k < NbinsnTT-1 ; ++k)
	    {
	      //Progression_1
          progressionCount++;
	      //cout<<"Doing progression "<<progressionCount<<endl;
	      Double_t Efficiency_Progression_1 = FindEfficiency_Progression((hardcodedIetBins2[j]+hardcodedIetBins2[j+1])/2., 25., 0.1, 50.);
	      if(Efficiency_Progression_1>=0.9999) Efficiency_Progression_1 = 1.0001;
	      Int_t Int_Efficiency_Progression_1 = int(Efficiency_Progression_1*100);
	      ostringstream convert_Progression_1;
	      convert_Progression_1 << Int_Efficiency_Progression_1 ;
	      TString Result_Progression_1 = TString(convert_Progression_1.str());
	      Int_t IsoCut_Progression_1 = histosIsolation[Result_Progression_1]->GetBinContent(i+1,FindBinCorrespondenceIEt(hardcodedIetBins2[j])+1,FindBinCorrespondencenTT(hardcodednTTBins2[k])+1);
	      if(Int_Efficiency_Progression_1==100) IsoCut_Progression_1 = 1000;
	      LUT_Progression_1->SetBinContent(i+1,j+1,k+1,IsoCut_Progression_1);


	      //Progression_2
          progressionCount++;
	      //cout<<"Doing progression "<<progressionCount<<endl;
	      Double_t Efficiency_Progression_2 = FindEfficiency_Progression((hardcodedIetBins2[j]+hardcodedIetBins2[j+1])/2., 25., 0.4, 50.);
	      if(Efficiency_Progression_2>=0.9999) Efficiency_Progression_2 = 1.0001;
	      Int_t Int_Efficiency_Progression_2 = int(Efficiency_Progression_2*100);
	      ostringstream convert_Progression_2;
	      convert_Progression_2 << Int_Efficiency_Progression_2 ;
	      TString Result_Progression_2 = TString(convert_Progression_2.str());
	      Int_t IsoCut_Progression_2 = histosIsolation[Result_Progression_2]->GetBinContent(i+1,FindBinCorrespondenceIEt(hardcodedIetBins2[j])+1,FindBinCorrespondencenTT(hardcodednTTBins2[k])+1);
	      if(Int_Efficiency_Progression_2==100) IsoCut_Progression_2 = 1000;
	      LUT_Progression_2->SetBinContent(i+1,j+1,k+1,IsoCut_Progression_2);

	      //	      cout<<"after progression 2"<<endl;

	      //Progression_3
          progressionCount++;
	      //cout<<"Doing progression "<<progressionCount<<endl;
	      Double_t Efficiency_Progression_3 = FindEfficiency_Progression((hardcodedIetBins2[j]+hardcodedIetBins2[j+1])/2., 25., 0.5, 50.);
	      if(Efficiency_Progression_3>=0.9999) Efficiency_Progression_3 = 1.0001;
	      Int_t Int_Efficiency_Progression_3 = int(Efficiency_Progression_3*100);
	      ostringstream convert_Progression_3;
	      convert_Progression_3 << Int_Efficiency_Progression_3 ;
	      TString Result_Progression_3 = TString(convert_Progression_3.str());
	      Int_t IsoCut_Progression_3 = histosIsolation[Result_Progression_3]->GetBinContent(i+1,FindBinCorrespondenceIEt(hardcodedIetBins2[j])+1,FindBinCorrespondencenTT(hardcodednTTBins2[k])+1);
	      if(Int_Efficiency_Progression_3==100) IsoCut_Progression_3 = 1000;
	      LUT_Progression_3->SetBinContent(i+1,j+1,k+1,IsoCut_Progression_3);


	      //Progression_4
          progressionCount++;
	      //cout<<"Doing progression "<<progressionCount<<endl;
	      Double_t Efficiency_Progression_4 = FindEfficiency_Progression((hardcodedIetBins2[j]+hardcodedIetBins2[j+1])/2., 25., 0.6, 50.);
	      if(Efficiency_Progression_4>=0.9999) Efficiency_Progression_4 = 1.0001;
	      Int_t Int_Efficiency_Progression_4 = int(Efficiency_Progression_4*100);
	      ostringstream convert_Progression_4;
	      convert_Progression_4 << Int_Efficiency_Progression_4 ;
	      TString Result_Progression_4 = TString(convert_Progression_4.str());
	      Int_t IsoCut_Progression_4 = histosIsolation[Result_Progression_4]->GetBinContent(i+1,FindBinCorrespondenceIEt(hardcodedIetBins2[j])+1,FindBinCorrespondencenTT(hardcodednTTBins2[k])+1);
	      if(Int_Efficiency_Progression_4==100) IsoCut_Progression_4 = 1000;
	      LUT_Progression_4->SetBinContent(i+1,j+1,k+1,IsoCut_Progression_4);

	      // if(j==3)
	      // 	{
	      // 	  cout<<"iet bin (4 bits) = "<<FindBinCorrespondenceIEt(hardcodedIetBins2[j])<<endl;
	      // 	  cout<<"iet bin (6 bits) = "<<j<<endl;
	      // 	  cout<<"ieta bin = "<<i<<endl;
	      // 	  cout<<"nTT bin = "<<k<<endl;
	      // 	  cout<<"iso = "<<IsoCut_Progression_4<<endl;
	      // 	}

	      //Progression_5
          progressionCount++;
	      //cout<<"Doing progression "<<progressionCount<<endl;
	      Double_t Efficiency_Progression_5 = FindEfficiency_Progression((hardcodedIetBins2[j]+hardcodedIetBins2[j+1])/2., 25., 0.7, 50.);
	      if(Efficiency_Progression_5>=0.9999) Efficiency_Progression_5 = 1.0001;
	      Int_t Int_Efficiency_Progression_5 = int(Efficiency_Progression_5*100);
	      ostringstream convert_Progression_5;
	      convert_Progression_5 << Int_Efficiency_Progression_5 ;
	      TString Result_Progression_5 = TString(convert_Progression_5.str());
	      Int_t IsoCut_Progression_5 = histosIsolation[Result_Progression_5]->GetBinContent(i+1,FindBinCorrespondenceIEt(hardcodedIetBins2[j])+1,FindBinCorrespondencenTT(hardcodednTTBins2[k])+1);
	      if(Int_Efficiency_Progression_5==100) IsoCut_Progression_5 = 1000;
	      LUT_Progression_5->SetBinContent(i+1,j+1,k+1,IsoCut_Progression_5);

	      //Progression_6
          progressionCount++;
	      //cout<<"Doing progression "<<progressionCount<<endl;
	      Double_t Efficiency_Progression_6 = FindEfficiency_Progression((hardcodedIetBins2[j]+hardcodedIetBins2[j+1])/2., 25., 0.5, 40.);
	      if(Efficiency_Progression_6>=0.9999) Efficiency_Progression_6 = 1.0001;
	      Int_t Int_Efficiency_Progression_6 = int(Efficiency_Progression_6*100);
	      ostringstream convert_Progression_6;
	      convert_Progression_6 << Int_Efficiency_Progression_6 ;
	      TString Result_Progression_6 = TString(convert_Progression_6.str());
	      Int_t IsoCut_Progression_6 = histosIsolation[Result_Progression_6]->GetBinContent(i+1,FindBinCorrespondenceIEt(hardcodedIetBins2[j])+1,FindBinCorrespondencenTT(hardcodednTTBins2[k])+1);
	      if(Int_Efficiency_Progression_6==100) IsoCut_Progression_6 = 1000;
	      LUT_Progression_6->SetBinContent(i+1,j+1,k+1,IsoCut_Progression_6);

	      //	      cout<<"after progression 6"<<endl;

	      //Progression_7
          progressionCount++;
	      //cout<<"Doing progression "<<progressionCount<<endl;
	      Double_t Efficiency_Progression_7 = FindEfficiency_Progression((hardcodedIetBins2[j]+hardcodedIetBins2[j+1])/2., 25., 0.5, 40.);
	      if(Efficiency_Progression_7>=0.9999) Efficiency_Progression_7 = 1.0001;
	      Int_t Int_Efficiency_Progression_7 = int(Efficiency_Progression_7*100);
	      ostringstream convert_Progression_7;
	      convert_Progression_7 << Int_Efficiency_Progression_7 ;
	      TString Result_Progression_7 = TString(convert_Progression_7.str());
	      Int_t IsoCut_Progression_7 = histosIsolation[Result_Progression_7]->GetBinContent(i+1,FindBinCorrespondenceIEt(hardcodedIetBins2[j])+1,FindBinCorrespondencenTT(hardcodednTTBins2[k])+1);
	      if(Int_Efficiency_Progression_7==100) IsoCut_Progression_7 = 1000;
	      LUT_Progression_7->SetBinContent(i+1,j+1,k+1,IsoCut_Progression_7);

	      //Progression_8
          progressionCount++;
	      //cout<<"Doing progression "<<progressionCount<<endl;
	      Double_t Efficiency_Progression_8 = FindEfficiency_Progression((hardcodedIetBins2[j]+hardcodedIetBins2[j+1])/2., 25., 0.5, 60.);
	      if(Efficiency_Progression_8>=0.9999) Efficiency_Progression_8 = 1.0001;
	      Int_t Int_Efficiency_Progression_8 = int(Efficiency_Progression_8*100);
	      ostringstream convert_Progression_8;
	      convert_Progression_8 << Int_Efficiency_Progression_8 ;
	      TString Result_Progression_8 = TString(convert_Progression_8.str());
	      Int_t IsoCut_Progression_8 = histosIsolation[Result_Progression_8]->GetBinContent(i+1,FindBinCorrespondenceIEt(hardcodedIetBins2[j])+1,FindBinCorrespondencenTT(hardcodednTTBins2[k])+1);
	      if(Int_Efficiency_Progression_8==100) IsoCut_Progression_8 = 1000;
	      LUT_Progression_8->SetBinContent(i+1,j+1,k+1,IsoCut_Progression_8);

	      //Progression_9
          progressionCount++;
	      //cout<<"Doing progression "<<progressionCount<<endl;
	      Double_t Efficiency_Progression_9 = FindEfficiency_Progression((hardcodedIetBins2[j]+hardcodedIetBins2[j+1])/2., 25., 0.6, 40.);
	      if(Efficiency_Progression_9>=0.9999) Efficiency_Progression_9 = 1.0001;
	      Int_t Int_Efficiency_Progression_9 = int(Efficiency_Progression_9*100);
	      ostringstream convert_Progression_9;
	      convert_Progression_9 << Int_Efficiency_Progression_9 ;
	      TString Result_Progression_9 = TString(convert_Progression_9.str());
	      Int_t IsoCut_Progression_9 = histosIsolation[Result_Progression_9]->GetBinContent(i+1,FindBinCorrespondenceIEt(hardcodedIetBins2[j])+1,FindBinCorrespondencenTT(hardcodednTTBins2[k])+1);
	      if(Int_Efficiency_Progression_9==100) IsoCut_Progression_9 = 1000;
	      LUT_Progression_9->SetBinContent(i+1,j+1,k+1,IsoCut_Progression_9);

	      //Progression_10
          progressionCount++;
	      //cout<<"Doing progression "<<progressionCount<<endl;
	      Double_t Efficiency_Progression_10 = FindEfficiency_Progression((hardcodedIetBins2[j]+hardcodedIetBins2[j+1])/2., 25., 0.7, 40.);
	      if(Efficiency_Progression_10>=0.9999) Efficiency_Progression_10 = 1.0001;
	      Int_t Int_Efficiency_Progression_10 = int(Efficiency_Progression_10*100);
	      ostringstream convert_Progression_10;
	      convert_Progression_10 << Int_Efficiency_Progression_10 ;
	      TString Result_Progression_10 = TString(convert_Progression_10.str());
	      Int_t IsoCut_Progression_10 = histosIsolation[Result_Progression_10]->GetBinContent(i+1,FindBinCorrespondenceIEt(hardcodedIetBins2[j])+1,FindBinCorrespondencenTT(hardcodednTTBins2[k])+1);
	      if(Int_Efficiency_Progression_10==100) IsoCut_Progression_10 = 1000;
	      LUT_Progression_10->SetBinContent(i+1,j+1,k+1,IsoCut_Progression_10);

	      //Progression_11
          progressionCount++;
	      //cout<<"Doing progression "<<progressionCount<<endl;
	      Double_t Efficiency_Progression_11 = FindEfficiency_Progression((hardcodedIetBins2[j]+hardcodedIetBins2[j+1])/2., 20., 0.7, 40.);
	      if(Efficiency_Progression_11>=0.9999) Efficiency_Progression_11 = 1.0001;
	      Int_t Int_Efficiency_Progression_11 = int(Efficiency_Progression_11*100);
	      ostringstream convert_Progression_11;
	      convert_Progression_11 << Int_Efficiency_Progression_11 ;
	      TString Result_Progression_11 = TString(convert_Progression_11.str());
	      Int_t IsoCut_Progression_11 = histosIsolation[Result_Progression_11]->GetBinContent(i+1,FindBinCorrespondenceIEt(hardcodedIetBins2[j])+1,FindBinCorrespondencenTT(hardcodednTTBins2[k])+1);
	      if(Int_Efficiency_Progression_11==100) IsoCut_Progression_11 = 1000;
	      LUT_Progression_11->SetBinContent(i+1,j+1,k+1,IsoCut_Progression_11);

	      //	      cout<<"after progression 11"<<endl;

	      //Progression_12
          progressionCount++;
	      //cout<<"Doing progression "<<progressionCount<<endl;
	      Double_t Efficiency_Progression_12 = FindEfficiency_Progression((hardcodedIetBins2[j]+hardcodedIetBins2[j+1])/2., 20., 0.7, 35.);
	      if(Efficiency_Progression_12>=0.9999) Efficiency_Progression_12 = 1.0001;
	      Int_t Int_Efficiency_Progression_12 = int(Efficiency_Progression_12*100);
	      ostringstream convert_Progression_12;
	      convert_Progression_12 << Int_Efficiency_Progression_12 ;
	      TString Result_Progression_12 = TString(convert_Progression_12.str());
	      Int_t IsoCut_Progression_12 = histosIsolation[Result_Progression_12]->GetBinContent(i+1,FindBinCorrespondenceIEt(hardcodedIetBins2[j])+1,FindBinCorrespondencenTT(hardcodednTTBins2[k])+1);
	      if(Int_Efficiency_Progression_12==100) IsoCut_Progression_12 = 1000;
	      LUT_Progression_12->SetBinContent(i+1,j+1,k+1,IsoCut_Progression_12);

	      //Progression_13
          progressionCount++;
	      //cout<<"Doing progression "<<progressionCount<<endl;
	      Double_t Efficiency_Progression_13 = FindEfficiency_Progression((hardcodedIetBins2[j]+hardcodedIetBins2[j+1])/2., 20., 0.6, 35.);
	      if(Efficiency_Progression_13>=0.9999) Efficiency_Progression_13 = 1.0001;
	      Int_t Int_Efficiency_Progression_13 = int(Efficiency_Progression_13*100);
	      ostringstream convert_Progression_13;
	      convert_Progression_13 << Int_Efficiency_Progression_13 ;
	      TString Result_Progression_13 = TString(convert_Progression_13.str());
	      Int_t IsoCut_Progression_13 = histosIsolation[Result_Progression_13]->GetBinContent(i+1,FindBinCorrespondenceIEt(hardcodedIetBins2[j])+1,FindBinCorrespondencenTT(hardcodednTTBins2[k])+1);
	      if(Int_Efficiency_Progression_13==100) IsoCut_Progression_13 = 1000;
	      LUT_Progression_13->SetBinContent(i+1,j+1,k+1,IsoCut_Progression_13);

	      //Progression_14
          progressionCount++;
	      //cout<<"Doing progression "<<progressionCount<<endl;
	      Double_t Efficiency_Progression_14 = FindEfficiency_Progression((hardcodedIetBins2[j]+hardcodedIetBins2[j+1])/2., 20., 0.5, 35.);
	      if(Efficiency_Progression_14>=0.9999) Efficiency_Progression_14 = 1.0001;
	      Int_t Int_Efficiency_Progression_14 = int(Efficiency_Progression_14*100);
	      ostringstream convert_Progression_14;
	      convert_Progression_14 << Int_Efficiency_Progression_14 ;
	      TString Result_Progression_14 = TString(convert_Progression_14.str());
	      Int_t IsoCut_Progression_14 = histosIsolation[Result_Progression_14]->GetBinContent(i+1,FindBinCorrespondenceIEt(hardcodedIetBins2[j])+1,FindBinCorrespondencenTT(hardcodednTTBins2[k])+1);
	      if(Int_Efficiency_Progression_14==100) IsoCut_Progression_14 = 1000;
	      LUT_Progression_14->SetBinContent(i+1,j+1,k+1,IsoCut_Progression_14);

	      //Progression_15
          progressionCount++;
	      //cout<<"Doing progression "<<progressionCount<<endl;
	      Double_t Efficiency_Progression_15 = FindEfficiency_Progression((hardcodedIetBins2[j]+hardcodedIetBins2[j+1])/2., 20., 0.7, 50.);
	      if(Efficiency_Progression_15>=0.9999) Efficiency_Progression_15 = 1.0001;
	      Int_t Int_Efficiency_Progression_15 = int(Efficiency_Progression_15*100);
	      ostringstream convert_Progression_15;
	      convert_Progression_15 << Int_Efficiency_Progression_15 ;
	      TString Result_Progression_15 = TString(convert_Progression_15.str());
	      Int_t IsoCut_Progression_15 = histosIsolation[Result_Progression_15]->GetBinContent(i+1,FindBinCorrespondenceIEt(hardcodedIetBins2[j])+1,FindBinCorrespondencenTT(hardcodednTTBins2[k])+1);
	      if(Int_Efficiency_Progression_15==100) IsoCut_Progression_15 = 1000;
	      LUT_Progression_15->SetBinContent(i+1,j+1,k+1,IsoCut_Progression_15);

	      //	      cout<<"after progression 15"<<endl;

	      //Progression_16
          progressionCount++;
	      //cout<<"Doing progression "<<progressionCount<<endl;
	      Double_t Efficiency_Progression_16 = FindEfficiency_Progression((hardcodedIetBins2[j]+hardcodedIetBins2[j+1])/2., 23., 0., 30.);
	      // Double_t Efficiency_Progression_16 = FindEfficiency_Progression((hardcodedIetBins2[j]+hardcodedIetBins2[j+1])/2., 20., 0., 30.);
	      //cout<<"Efficiency_Progression_16 = "<<Efficiency_Progression_16<<endl;
	      if(Efficiency_Progression_16>=0.9999) Efficiency_Progression_16 = 1.0001;
	      Int_t Int_Efficiency_Progression_16 = int(Efficiency_Progression_16*100);
	      ostringstream convert_Progression_16;
	      convert_Progression_16 << Int_Efficiency_Progression_16 ;
	      TString Result_Progression_16 = TString(convert_Progression_16.str());
	      Int_t IsoCut_Progression_16 = histosIsolation[Result_Progression_16]->GetBinContent(i+1,FindBinCorrespondenceIEt(hardcodedIetBins2[j])+1,FindBinCorrespondencenTT(hardcodednTTBins2[k])+1);
	      if(Int_Efficiency_Progression_16==100) IsoCut_Progression_16 = 1000;
	      LUT_Progression_16->SetBinContent(i+1,j+1,k+1,IsoCut_Progression_16);

	      //Progression_17
          progressionCount++;
	      //cout<<"Doing progression "<<progressionCount<<endl;
	      Double_t Efficiency_Progression_17 = FindEfficiency_Progression((hardcodedIetBins2[j]+hardcodedIetBins2[j+1])/2., 15., 0., 30.);
	      if(Efficiency_Progression_17>=0.9999) Efficiency_Progression_17 = 1.0001;
	      Int_t Int_Efficiency_Progression_17 = int(Efficiency_Progression_17*100);
	      ostringstream convert_Progression_17;
	      convert_Progression_17 << Int_Efficiency_Progression_17 ;
	      TString Result_Progression_17 = TString(convert_Progression_17.str());
	      Int_t IsoCut_Progression_17 = histosIsolation[Result_Progression_17]->GetBinContent(i+1,FindBinCorrespondenceIEt(hardcodedIetBins2[j])+1,FindBinCorrespondencenTT(hardcodednTTBins2[k])+1);
	      if(Int_Efficiency_Progression_17==100) IsoCut_Progression_17 = 1000;
	      LUT_Progression_17->SetBinContent(i+1,j+1,k+1,IsoCut_Progression_17);

	      //Progression_18
          progressionCount++;
	      //cout<<"Doing progression "<<progressionCount<<endl;
	      Double_t Efficiency_Progression_18 = FindEfficiency_Progression((hardcodedIetBins2[j]+hardcodedIetBins2[j+1])/2., 20., 0., 35.);
	      if(Efficiency_Progression_18>=0.9999) Efficiency_Progression_18 = 1.0001;
	      Int_t Int_Efficiency_Progression_18 = int(Efficiency_Progression_18*100);
	      ostringstream convert_Progression_18;
	      convert_Progression_18 << Int_Efficiency_Progression_18 ;
	      TString Result_Progression_18 = TString(convert_Progression_18.str());
	      Int_t IsoCut_Progression_18 = histosIsolation[Result_Progression_18]->GetBinContent(i+1,FindBinCorrespondenceIEt(hardcodedIetBins2[j])+1,FindBinCorrespondencenTT(hardcodednTTBins2[k])+1);
	      if(Int_Efficiency_Progression_18==100) IsoCut_Progression_18 = 1000;
	      LUT_Progression_18->SetBinContent(i+1,j+1,k+1,IsoCut_Progression_18);
	      
	      
	      //Progression_19
          progressionCount++;
	      //cout<<"Doing progression "<<progressionCount<<endl;
	      Double_t Efficiency_Progression_19 = FindEfficiency_Progression((hardcodedIetBins2[j]+hardcodedIetBins2[j+1])/2., 15., 0., 35.);
	      //cout<<"Efficiency_Progression_19 = "<<Efficiency_Progression_19<<endl;
	      if(Efficiency_Progression_19>=0.9999) Efficiency_Progression_19 = 1.0001;
	      Int_t Int_Efficiency_Progression_19 = int(Efficiency_Progression_19*100);
	      ostringstream convert_Progression_19;
	      convert_Progression_19 << Int_Efficiency_Progression_19 ;
	      TString Result_Progression_19 = TString(convert_Progression_19.str());
	      Int_t IsoCut_Progression_19 = histosIsolation[Result_Progression_19]->GetBinContent(i+1,FindBinCorrespondenceIEt(hardcodedIetBins2[j])+1,FindBinCorrespondencenTT(hardcodednTTBins2[k])+1);
	      if(Int_Efficiency_Progression_19==100) IsoCut_Progression_19 = 1000;
	      LUT_Progression_19->SetBinContent(i+1,j+1,k+1,IsoCut_Progression_19);

	      //Progression_20
	      Double_t Efficiency_Progression_20 = FindEfficiency_Progression((hardcodedIetBins2[j]+hardcodedIetBins2[j+1])/2., 10., 0., 35.);
	      if(Efficiency_Progression_20>=0.9999) Efficiency_Progression_20 = 1.0001;
	      Int_t Int_Efficiency_Progression_20 = int(Efficiency_Progression_20*100);
	      ostringstream convert_Progression_20;
	      convert_Progression_20 << Int_Efficiency_Progression_20 ;
	      TString Result_Progression_20 = TString(convert_Progression_20.str());
	      Int_t IsoCut_Progression_20 = histosIsolation[Result_Progression_20]->GetBinContent(i+1,FindBinCorrespondenceIEt(hardcodedIetBins2[j])+1,FindBinCorrespondencenTT(hardcodednTTBins2[k])+1);
	      if(Int_Efficiency_Progression_20==100) IsoCut_Progression_20 = 1000;
	      LUT_Progression_20->SetBinContent(i+1,j+1,k+1,IsoCut_Progression_20);

	      //	      cout<<"after progression 20"<<endl;

	      //Progression_21
	      //cout<<"i = "<<i<<endl;
	      //	      cout<<"j = "<<j<<endl;
	      //	      cout<<"k = "<<k<<endl;
          progressionCount++;
	      //cout<<"Doing progression "<<progressionCount<<endl;
	      Double_t Efficiency_Progression_21 = FindEfficiency_Progression((hardcodedIetBins2[j]+hardcodedIetBins2[j+1])/2., 20., 0.80, 50.);
	      if(Efficiency_Progression_21>=0.9999) Efficiency_Progression_21 = 1.0001;
	      Int_t Int_Efficiency_Progression_21 = int(Efficiency_Progression_21*100);
	      //	      cout<<"Int_Efficiency_Progression_21 = "<<Int_Efficiency_Progression_21<<endl;
	      ostringstream convert_Progression_21;
	      convert_Progression_21 << Int_Efficiency_Progression_21 ;
	      TString Result_Progression_21 = TString(convert_Progression_21.str());
	      Int_t IsoCut_Progression_21 = histosIsolation[Result_Progression_21]->GetBinContent(i+1,FindBinCorrespondenceIEt(hardcodedIetBins2[j])+1,FindBinCorrespondencenTT(hardcodednTTBins2[k])+1);
	      if(Int_Efficiency_Progression_21==100) IsoCut_Progression_21 = 1000;
	      //	      cout<<"IsoCut_Progression_21 = "<<IsoCut_Progression_21<<endl;

	     /* if(hardcodedIetBins2[j]/2.>30.&&hardcodedIetBins2[j]/2.<55. && k==10 && i==0)
		{
		  cout<<"ieta bin = "<<i<<endl;
		  cout<<"pT bin = "<<j<<endl;
		  cout<<"nTT bin = "<<k<<endl;

		  cout<<"looking for..."<<endl;
		  cout<<"ieta bin = "<<i<<endl;
		  cout<<"pT bin = "<<j<<endl;
		  cout<<"nTT bin = "<<k<<endl;

		  cout<<"pT = "<<hardcodedIetBins2[j]/2.<<endl;
		  // cout<<"pT = "<<hardcodedIetBins[FindBinCorrespondenceIEt((hardcodedIetBins2[j]+hardcodedIetBins2[j+1])/2.)]/2.<<endl;
		  cout<<"Target efficiency = "<<Efficiency_Progression_5<<endl;
		  cout<<"Isolation cut = "<<IsoCut_Progression_5<<endl;
		  cout<<"--"<<endl;
		}*/
	      
	//      if(IsoCut_Progression_21==0)
	//	{
	//	  cout<<"j compressed = "<<FindBinCorrespondenceIEt(hardcodedIetBins2[j])<<endl;
	//	  cout<<"k compressed = "<<FindBinCorrespondencenTT(hardcodednTTBins2[k])<<endl;
	//	}
	      LUT_Progression_21->SetBinContent(i+1,j+1,k+1,IsoCut_Progression_21);

          progressionCount++;
	      //cout<<"Doing progression "<<progressionCount<<endl;
	      Double_t Efficiency_Progression_22 = FindEfficiency_Progression((hardcodedIetBins2[j]+hardcodedIetBins2[j+1])/2., 25., 0.70, 70.);
	      //cout<<"Bin j = "<<j<<", TargetEfficiency = "<<Efficiency_Progression_22<<endl;
	      // Double_t Efficiency_Progression_22 = FindEfficiency_Progression((hardcodedIetBins2[j]+hardcodedIetBins2[j+1])/2., 20., 0.60, 80.);
	      if(Efficiency_Progression_22>=0.9999) Efficiency_Progression_22 = 1.0001;
	      Int_t Int_Efficiency_Progression_22 = int(Efficiency_Progression_22*100);
	      ostringstream convert_Progression_22;
	      convert_Progression_22 << Int_Efficiency_Progression_22 ;
	      TString Result_Progression_22 = TString(convert_Progression_22.str());
	      Int_t IsoCut_Progression_22 = histosIsolation[Result_Progression_22]->GetBinContent(i+1,FindBinCorrespondenceIEt(hardcodedIetBins2[j])+1,FindBinCorrespondencenTT(hardcodednTTBins2[k])+1);
	      if(Int_Efficiency_Progression_22==100) IsoCut_Progression_22 = 1000;

	      LUT_Progression_22->SetBinContent(i+1,j+1,k+1,IsoCut_Progression_22);

	      //	      cout<<"after progression 21"<<endl;

	      // LUT_WP

	      for(UInt_t iEff = 0 ; iEff < 101 ; ++iEff)
		{
		  //		  cout<<"iEff = "<<iEff<<endl;
		  Double_t Efficiency = iEff/100.;
		  if(Efficiency>=0.9999) Efficiency = 1.0001;
		  Int_t Int_Efficiency = int(Efficiency*100);
		  ostringstream convert;
		  convert << Int_Efficiency ;
		  TString Result = TString(convert.str());
		  //		  cout<<"Result = "<<Result<<endl;
		  Int_t IsoCut = histosIsolation[Result]->GetBinContent(i+1,FindBinCorrespondenceIEt(hardcodedIetBins2[j])+1,FindBinCorrespondencenTT(hardcodednTTBins2[k])+1);
		  if(Int_Efficiency==100) IsoCut = 1000;
		  //		  cout<<"IsoCut = "<<IsoCut<<endl;
		  LUT_WP.at(iEff)->SetBinContent(i+1,j+1,k+1,IsoCut);
		  //		  cout<<"after at.iEff"<<endl;
		}

	      /*
	      //WP90
	      Double_t Efficiency_WP90 = 0.9;
	      if(Efficiency_WP90>=0.9999) Efficiency_WP90 = 1.0001;
	      Int_t Int_Efficiency_WP90 = int(Efficiency_WP90*100);
	      ostringstream convert_WP90;
	      convert_WP90 << Int_Efficiency_WP90 ;
	      TString Result_WP90 = TString(convert_WP90.str());
	      Int_t IsoCut_WP90 = histosIsolation[Result_WP90]->GetBinContent(i+1,FindBinCorrespondenceIEt(hardcodedIetBins2[j])+1,k+1);
	      if(Int_Efficiency_WP90==100) IsoCut_WP90 = 1000;
	      LUT_WP90->SetBinContent(i+1,j+1,k+1,IsoCut_WP90);

	      //WP80
	      Double_t Efficiency_WP80 = 0.8;
	      if(Efficiency_WP80>=0.9999) Efficiency_WP80 = 1.0001;
	      Int_t Int_Efficiency_WP80 = int(Efficiency_WP80*100);
	      ostringstream convert_WP80;
	      convert_WP80 << Int_Efficiency_WP80 ;
	      TString Result_WP80 = TString(convert_WP80.str());
	      Int_t IsoCut_WP80 = histosIsolation[Result_WP80]->GetBinContent(i+1,FindBinCorrespondenceIEt(hardcodedIetBins2[j])+1,k+1);
	      if(Int_Efficiency_WP80==100) IsoCut_WP80 = 1000;
	      LUT_WP80->SetBinContent(i+1,j+1,k+1,IsoCut_WP80);

	      //WP70
	      Double_t Efficiency_WP70 = 0.7;
	      if(Efficiency_WP70>=0.9999) Efficiency_WP70 = 1.0001;
	      Int_t Int_Efficiency_WP70 = int(Efficiency_WP70*100);
	      ostringstream convert_WP70;
	      convert_WP70 << Int_Efficiency_WP70 ;
	      TString Result_WP70 = TString(convert_WP70.str());
	      Int_t IsoCut_WP70 = histosIsolation[Result_WP70]->GetBinContent(i+1,FindBinCorrespondenceIEt(hardcodedIetBins2[j])+1,k+1);
	      if(Int_Efficiency_WP70==100) IsoCut_WP70 = 1000;
	      LUT_WP70->SetBinContent(i+1,j+1,k+1,IsoCut_WP70);
	      */
	    }
	}
    }

  LUTs_Options.cd();
 // TFile LUTs_Options("Iso_LUTs_Options_ggH_DY_NewCompression_NewLayer1_Thomas_nTT_18_07_16_OlivierFlatWP_SK1616.root","RECREATE");
  LUT_Progression_1->Write();
  LUT_Progression_2->Write();
  LUT_Progression_3->Write();
  LUT_Progression_4->Write();
  LUT_Progression_5->Write();
  LUT_Progression_6->Write();
  LUT_Progression_7->Write();
  LUT_Progression_8->Write();
  LUT_Progression_9->Write();
  LUT_Progression_10->Write();
  LUT_Progression_11->Write();
  LUT_Progression_12->Write();
  LUT_Progression_13->Write();
  LUT_Progression_14->Write();
  LUT_Progression_15->Write();
  LUT_Progression_16->Write();
  LUT_Progression_17->Write();
  LUT_Progression_18->Write();
  LUT_Progression_19->Write();
  LUT_Progression_20->Write();
  LUT_Progression_21->Write();
  LUT_Progression_22->Write();

  for(UInt_t iEff = 0 ; iEff < 101 ; ++iEff)
    {
      LUT_WP.at(iEff)->Write();
    }

  // LUT_WP90->Write();
  // LUT_WP80->Write();
  // LUT_WP70->Write();
}
