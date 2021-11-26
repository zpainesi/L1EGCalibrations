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
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TProfile.h>
#include <sstream>
#include <TBranchElement.h>
#include <TGraphAsymmErrors.h>
#include <fstream>

// Regular LUT COPMRESSION SCHEME
// const Int_t NbinsIEta = 32+1; const Int_t IEtaOffset=15;
// const Int_t NbinsIEt  = 16+1; const Int_t IEtOffset = 0;
// const Int_t NbinsnTT  = 32+1; const Int_t INTTOffset= 0;//Acceptable LUT Thomas

// MODIFIED LUT COPMRESSION SCHEME
// Eta :  4 bit
// E   :  4 bit
// nTT :  5 bit, max 12

// const Int_t NbinsIEta = 4+1;  const Int_t IEtaOffset=1;
// const Int_t NbinsIEt  = 16+1; const Int_t IEtOffset = 0;
// const Int_t NbinsnTT  = 12+1; const Int_t INTTOffset= 0;//Acceptable LUT Thomas
// 
// 
// const Int_t NbinsIEt2 = 32+1;
// const Int_t NbinsnTT2 = 32+1;

 
// TAU LUT COMPRESSION SCHEME

const Int_t NbinsIEta = 4+1;  const Int_t IEtaOffset= 0;                        
const Int_t NbinsIEt = 32+1;  const Int_t IEtOffset = 0;
const Int_t NbinsnTT = 32+1;  const Int_t INTTOffset= 0;//Acceptable LUT Thomas//Acceptable LUT Thomas

const Int_t NbinsIEt2 = 32+1;
const Int_t NbinsnTT2 = 32+1;


using namespace std;

void Build_Isolation_WPs(TString inFile,TString oFile="Iso_LUTs_Distributions.root",Long64_t maxEvents=1000)
{
  TChain data("EGIsoCalibration");
  data.Add(inFile);
  TFile* f_out = new TFile(oFile,"RECREATE");
  

  Long64_t nentries=-1;
  int   L1Tau_IEt      = -99;
  int   L1Tau_compressed_IEt      = -99;
  int   L1Tau_IEta     = -99;
  int   L1Tau_compressed_IEta     = -99;
  int   L1Tau_compressed_IEta_offseted     = -99 + IEtaOffset;
  int   L1Tau_hasEM    = -99;
  int   L1Tau_isMerged = -99;
  int   L1Tau_nTT      = -99;
  int   L1Tau_Iso      = -99;
  // float L1Tau_Calibrated_pt = -99.;
  // float L1Tau_UnCalibrated_pt = -99.;
  // float L1Tau_CalibrationConstant = -99.;
  float OfflineTau_pt = -99.;
  bool OfflineTau_isMatched = false;
  int supercompressedE = -99;
  int supercompressednTT = -99;

   auto t_start = std::chrono::high_resolution_clock::now();
   auto t_end = std::chrono::high_resolution_clock::now();

  data.SetBranchAddress("l1RawE", &L1Tau_IEt);
  data.SetBranchAddress("ieta", &L1Tau_IEta);
  data.SetBranchAddress("compressedE", &supercompressedE);
  data.SetBranchAddress("compressedNTT", &supercompressednTT);
  data.SetBranchAddress("compressedieta", &L1Tau_compressed_IEta);
  data.SetBranchAddress("nTT", &L1Tau_nTT);
  data.SetBranchAddress("isoEt", &L1Tau_Iso);
  data.SetBranchAddress("offlineEta",&OfflineTau_pt);
  OfflineTau_isMatched = true;
  
  // data.SetBranchAddress("compressedE", &L1Tau_compressed_IEt);
  // data.SetBranchAddress("L1Tau_isMerged", &L1Tau_isMerged);
  // data.SetBranchAddress("L1Tau_hasEM", &L1Tau_hasEM);
  // data.SetBranchAddress("OfflineTau_isMatched",&OfflineTau_isMatched);
  // data.SetBranchAddress("L1Tau_Calibrated_pt", &L1Tau_Calibrated_pt);
  // data.SetBranchAddress("L1Tau_CalibrationConstant",&L1Tau_CalibrationConstant);

  std::map<TString,TH1F*> Histos_PerBin ;
  std::map<Int_t,TH3F*> IsoCut_PerBin ;
  std::map<Int_t,std::map<TString,Int_t>> IsoCut_PerEfficiency_PerBin;
  // std::map<TString,Int_t> IsoCut_PerBin_90pc ;
  // std::map<TString,Int_t> IsoCut_PerBin_80pc ;

  for(UInt_t i = 0 ; i < NbinsIEta-1 ; ++i)
    {
      for(UInt_t j = 0 ; j < NbinsIEt-1 ; ++j)
	    {
	      for(UInt_t k = 0 ; k < NbinsnTT-1 ; ++k)
	        {
	          TString Name_Histo = "Hist_";

	          stringstream ss_i;
	          ss_i << i;
	          TString Appendix_i = TString(ss_i.str());
	          Name_Histo += Appendix_i;
	          Name_Histo += "_";

	          stringstream ss_j;
	          ss_j << j;
	          TString Appendix_j = TString(ss_j.str());
	          Name_Histo += Appendix_j;
	          Name_Histo += "_";

	          stringstream ss_k;
	          ss_k << k;
	          TString Appendix_k = TString(ss_k.str());
	          Name_Histo += Appendix_k;

	          TH1F* temp_histo = new TH1F(Name_Histo.Data(),Name_Histo.Data(),100,0.,100.);
	          Histos_PerBin.insert(make_pair(Name_Histo,temp_histo));

	        //   cout<<"Name_Histo = "<<Name_Histo<<endl;
	        }
	   }
    }
  
  cout<<"entering loop"<<endl;

  TProfile* hprof_IEt  = new TProfile("hprof_IEt","Profile L1_Iso vs. L1_IEt",100,0.,200.,0,20);
  TProfile* hprof_IEta = new TProfile("hprof_IEta","Profile L1_Iso vs. L1_IEta",28,0.,28.,0,20);
  TProfile* hprof_nTT  = new TProfile("hprof_nTT","Profile L1_Iso vs. L1_IEta",70,0.,70.,0,20);

  nentries = data.GetEntries();
  if(maxEvents >0 ) nentries = nentries>maxEvents ? maxEvents : nentries;
  
  t_start = std::chrono::high_resolution_clock::now();
  t_end = std::chrono::high_resolution_clock::now();

  for(UInt_t i = 0 ; i < nentries ; ++i)
    {
      data.GetEntry(i);
      if(!OfflineTau_isMatched) continue;
      L1Tau_compressed_IEta_offseted= L1Tau_compressed_IEta + IEtaOffset ;
      if(i%10000 == 0) {
             t_end = std::chrono::high_resolution_clock::now();
             cout<<"Processing Entry in event loop : "<<i<<" / "<<nentries<<"  [ "<<100.0*i/nentries<<"  % ]  "
             << " Elapsed time : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()/1000.0
             <<"  Estimated time left : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()*( nentries - i)/i* 0.001
             <<endl;
            //cout<<"     supercompressedE = "<<supercompressedE<<endl;
            //cout<<"     L1Tau_compressed_IEta = "<<L1Tau_compressed_IEta_offseted<<endl;
            //cout<<"     supercompressednTT = "<<supercompressednTT<<endl;
            //cout<<"     L1Tau_nTT = "<<L1Tau_nTT<<endl;
            //
            //cout<<endl;
      }

      hprof_IEt->Fill(L1Tau_IEt,L1Tau_Iso,1);
      hprof_IEta->Fill(L1Tau_IEta,L1Tau_Iso,1);
      hprof_nTT->Fill(L1Tau_nTT,L1Tau_Iso,1);


      std::vector<Int_t> binForIsolation ;
      binForIsolation.clear();
      binForIsolation.push_back(L1Tau_compressed_IEta_offseted);
      binForIsolation.push_back(supercompressedE);
      binForIsolation.push_back(supercompressednTT);

      
      TString Name_Histo = "Hist_";

      stringstream ss_i;
      ss_i << binForIsolation.at(0);
      TString Appendix_i = TString(ss_i.str());
      Name_Histo += Appendix_i;
      Name_Histo += "_";

      stringstream ss_j;
      ss_j << binForIsolation.at(1);
      TString Appendix_j = TString(ss_j.str());
      Name_Histo += Appendix_j;
      Name_Histo += "_";

      stringstream ss_k;
      ss_k << binForIsolation.at(2);
      TString Appendix_k = TString(ss_k.str());
      Name_Histo += Appendix_k;


     // cout<<"Name_Histo = "<<Name_Histo<<endl;
      Histos_PerBin[Name_Histo]->Fill(L1Tau_Iso);
    }
  

  std::cout<<"Event Loop completed "<<endl;


  hprof_IEt->Write();
  hprof_IEta->Write();
  hprof_nTT->Write();

  Int_t NumberOfHistosWithLowStats = 0;

  std::cout<<"Allocating IsoCutForEfficiency Maps ";
  for(UInt_t iEff = 0 ; iEff < 101 ; ++iEff)
    {
      if(iEff%20==0) {
      
      std::cout<<"\nAllocating IsoCutForEfficiency Maps for  eff = "<<iEff;
      
      }
      std::map<TString,Int_t> temp;

      for(UInt_t i = 0 ; i < NbinsIEta-1 ; ++i)
	{
	  for(UInt_t j = 0 ; j < NbinsIEt-1 ; ++j)
	    {
	      for(UInt_t k = 0 ; k < NbinsnTT-1 ; ++k)
		{
		  TString Name_Histo = "Hist_";
		  
		  stringstream ss_i;
		  ss_i << i;
		  TString Appendix_i = TString(ss_i.str());
		  Name_Histo += Appendix_i;
		  Name_Histo += "_";
		  
		  stringstream ss_j;
		  ss_j << j;
		  TString Appendix_j = TString(ss_j.str());
		  Name_Histo += Appendix_j;
		  Name_Histo += "_";
		  
		  stringstream ss_k;
		  ss_k << k;
		  TString Appendix_k = TString(ss_k.str());
		  Name_Histo += Appendix_k;
		  
		  temp.insert(make_pair(Name_Histo,-1));
		}
	    }
	}
      IsoCut_PerEfficiency_PerBin.insert(make_pair(iEff,temp));
      temp.clear();
    }
  
   std::cout<<"\n";
  std::cout<<"Allocating 3-D Histograms for  IsoCut_PerBin \n";
  for(UInt_t iEff = 0 ; iEff < 101 ; ++iEff)
    {
      ostringstream convert;
      convert << iEff ;
      TString NameEff = "Eff_";
      TString TempStr(convert.str());
      NameEff += TempStr;
      TH3F* temp = new TH3F(NameEff.Data(),NameEff.Data(),NbinsIEta-1,0,NbinsIEta-1,NbinsIEt-1,0,NbinsIEt-1,NbinsnTT-1,0,NbinsnTT-1);
      IsoCut_PerBin.insert(make_pair(iEff,temp));
    }

  std::cout<<"Filling 3-D Histograms for  IsoCut_PerBin \n";

  int nMax = (NbinsIEta -1 )*(NbinsIEt -1)*(NbinsnTT -1) ;
  int currN=0;
  
  t_start = std::chrono::high_resolution_clock::now();
  t_end = std::chrono::high_resolution_clock::now();

  for(UInt_t i = 0 ; i < NbinsIEta-1 ; ++i)
    {
      for(UInt_t j = 0 ; j < NbinsIEt-1 ; ++j)
	{
	  for(UInt_t k = 0 ; k < NbinsnTT-1 ; ++k)
	    {
            if(currN%256 == 0) 
            {
                 t_end = std::chrono::high_resolution_clock::now();
                 std::cout<<"\t for eta,iEt,nTT : "<<i<<","<<j<<","<<k  <<" ( "<<currN<<" / "<<nMax <<" ) "
                   << " Elapsed time : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()/1000.0
                   <<"  Estimated time left : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()*( nMax - currN)/currN *  0.001
                   <<endl;
            } 
          currN++;

	      TString Name_Histo = "Hist_";

	      stringstream ss_i;
	      ss_i << i;
	      TString Appendix_i = TString(ss_i.str());
	      Name_Histo += Appendix_i;
	      Name_Histo += "_";

	      stringstream ss_j;
	      ss_j << j;
	      TString Appendix_j = TString(ss_j.str());
	      Name_Histo += Appendix_j;
	      Name_Histo += "_";

	      stringstream ss_k;
	      ss_k << k;
	      TString Appendix_k = TString(ss_k.str());
	      Name_Histo += Appendix_k;


	      for(UInt_t iEff = 0 ; iEff < 101 ; ++iEff)
		{
		  Float_t Efficiency = 0.01*iEff;

		  for(UInt_t iIso = 0 ; iIso < 100 ; ++iIso)
		    {
		      if(Histos_PerBin[Name_Histo]->Integral(1,iIso+1)/Histos_PerBin[Name_Histo]->Integral(1,100+1)>=Efficiency)
			{
			  if(IsoCut_PerEfficiency_PerBin[iEff][Name_Histo]==-1)
			    {
			      //cout<<"Efficiency = "<<Efficiency<<", bin = "<<Name_Histo<<", cut = "<<iIso<<endl;
			      //if(iEff==97) cout<<"IsoCut = "<<iIso<<endl;
			      IsoCut_PerEfficiency_PerBin[iEff][Name_Histo]=iIso;
			      IsoCut_PerBin[iEff]->SetBinContent(i+1,j+1,k+1,iIso);
			      //cout<<"IsoCut_PerEfficiency_PerBin[iEff][Name_Histo] = "<<IsoCut_PerEfficiency_PerBin[iEff][Name_Histo]<<endl;
			    }
			}
		    }
		}

	      // IsoCut_PerBin_90pc.insert(make_pair(Name_Histo,-1));
	      // IsoCut_PerBin_80pc.insert(make_pair(Name_Histo,-1));

	      // for(UInt_t iIso = 0 ; iIso < 100 ; ++iIso)
	      // 	{
	      // 	  // cout<<"testing iIso = "<<iIso<<endl;
	      // 	  // cout<<"Histos_PerBin[Name_Histo]->Integral(1,iIso+1) = "<<Histos_PerBin[Name_Histo]->Integral(1,iIso+1)<<endl;
	      // 	  // cout<<"Histos_PerBin[Name_Histo]->Integral(1,100+1) = "<<Histos_PerBin[Name_Histo]->Integral(1,100+1)<<endl;
	      // 	  if(Histos_PerBin[Name_Histo]->Integral(1,iIso+1)/Histos_PerBin[Name_Histo]->Integral(1,100+1)>=0.9)
	      // 	    {
	      // 	      if(IsoCut_PerBin_90pc[Name_Histo]==-1) IsoCut_PerBin_90pc[Name_Histo]=Histos_PerBin[Name_Histo]->GetBinLowEdge(iIso+1);
	      // 	    }
	      // 	  if(Histos_PerBin[Name_Histo]->Integral(1,iIso+1)/Histos_PerBin[Name_Histo]->Integral(1,100+1)>=0.8)
	      // 	    {
	      // 	      if(IsoCut_PerBin_80pc[Name_Histo]==-1) IsoCut_PerBin_80pc[Name_Histo]=Histos_PerBin[Name_Histo]->GetBinLowEdge(iIso+1);
	      // 	    }
	      // 	}

	      // cout<<"50% iso cut for "<<Name_Histo<<" is = "<<IsoCut_PerBin_90pc[Name_Histo]<<endl;

	      //if(Histos_PerBin[Name_Histo]->GetEntries()<20)
	      if(Histos_PerBin[Name_Histo]->GetEntries()<40)
		{
		  NumberOfHistosWithLowStats++;
	//	  cout<<"Histo "<<Name_Histo<<" has low stat.: "<<Histos_PerBin[Name_Histo]->GetEntries()<<endl;
		}

	      Histos_PerBin[Name_Histo]->Write();
	    }
	}

    }

  //Efficiency as function of pT
  std::cout<<"Making Efficiency Histograms as a function of pT"<<endl;

  std::map<Int_t,TH1F*> pt_pass_efficiency ;
  std::map<Int_t,TH1F*> pt_pass_efficiency_TH3 ;
  for(UInt_t iEff = 0 ; iEff < 101 ; ++iEff)
    {
      TString nameHisto = "pt_pass_efficiency_";
      ostringstream convert;
      convert << iEff;
      TString temp(convert.str());
      nameHisto += temp;

      TString nameHisto_TH3 = "pt_pass_efficiency_TH3_";
      nameHisto_TH3 += temp;
      
      TH1F* temp_histo = new TH1F(nameHisto.Data(),nameHisto.Data(),100,0,200);
      TH1F* temp_histo_TH3 = new TH1F(nameHisto_TH3.Data(),nameHisto_TH3.Data(),100,0,200);

      pt_pass_efficiency.insert(make_pair(iEff,temp_histo));
      pt_pass_efficiency_TH3.insert(make_pair(iEff,temp_histo_TH3));
    }

  std::cout<<"Making Efficiency Histograms as a function of eta"<<endl;
  std::map<Int_t,TH1F*> eta_pass_efficiency ;
  for(UInt_t iEff = 0 ; iEff < 101 ; ++iEff)
    {
      TString nameHisto = "eta_pass_efficiency_";
      ostringstream convert;
      convert << iEff;
      TString temp(convert.str());
      nameHisto += temp;
      
      TH1F* temp_histo = new TH1F(nameHisto.Data(),nameHisto.Data(),100,0,100);
      eta_pass_efficiency.insert(make_pair(iEff,temp_histo));
    }

  std::cout<<"Making Efficiency Histograms as a function of nTT"<<endl;
  std::map<Int_t,TH1F*> nTT_pass_efficiency ;
  for(UInt_t iEff = 0 ; iEff < 101 ; ++iEff)
    {
      TString nameHisto = "nTT_pass_efficiency_";
      ostringstream convert;
      convert << iEff;
      TString temp(convert.str());
      nameHisto += temp;
      
      TH1F* temp_histo = new TH1F(nameHisto.Data(),nameHisto.Data(),100,0,100);
      nTT_pass_efficiency.insert(make_pair(iEff,temp_histo));

    }

  TH1F* pt = new TH1F("pt","pt",100,0,200);
  TH1F* eta = new TH1F("eta","eta",100,0,100);
  TH1F* nTT = new TH1F("nTT","nTT",100,0,100);
  
   t_start = std::chrono::high_resolution_clock::now();
   t_end = std::chrono::high_resolution_clock::now();

  std::cout<<"Looing over events to fill up the Efficiecy Histograms based on the obtaied cuts "<<endl;
  for(UInt_t i = 0 ; i < nentries ; ++i)
    {
      if(i%10000 == 0) 
      {
           t_end = std::chrono::high_resolution_clock::now();
           cout<<"Processing Entry in event loop : "<<i<<" / "<<nentries<<"  [ "<<100.0*i/nentries<<"  % ]  "
             << " Elapsed time : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()/1000.0
             <<"  Estimated time left : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()*( nentries - i)/i* 0.001
             <<endl;
      }
      data.GetEntry(i);
      if(!OfflineTau_isMatched) continue;

      //cout<<"L1_IEta = "<<L1_IEta<<endl;
      // if(i>10) continue;

      std::vector<Int_t> binForIsolation ;
      binForIsolation.clear();
      binForIsolation.push_back(L1Tau_compressed_IEta_offseted);
      binForIsolation.push_back(supercompressedE);
      binForIsolation.push_back(supercompressednTT);
      //cout<<"nTT = "<<L1_nTT<<", bin = "<<binForIsolation.at(2)<<endl;

      TString Name_Histo = "Hist_";

      stringstream ss_i;
      ss_i << binForIsolation.at(0);
      TString Appendix_i = TString(ss_i.str());
      Name_Histo += Appendix_i;
      Name_Histo += "_";
      
      stringstream ss_j;
      ss_j << binForIsolation.at(1);
      TString Appendix_j = TString(ss_j.str());
      Name_Histo += Appendix_j;
      Name_Histo += "_";
      
      stringstream ss_k;
      ss_k << binForIsolation.at(2);
      TString Appendix_k = TString(ss_k.str());
      Name_Histo += Appendix_k;

      for(UInt_t iEff = 0 ; iEff < 101 ; ++iEff)
	{    
	  if(L1Tau_Iso<=IsoCut_PerEfficiency_PerBin[iEff][Name_Histo])
	    {
	      eta_pass_efficiency[iEff]->Fill(L1Tau_IEta);
	      pt_pass_efficiency[iEff]->Fill(L1Tau_IEt);
	      nTT_pass_efficiency[iEff]->Fill(L1Tau_nTT);
	    }
	  //cout<<"L1Tau_Iso = "<<L1Tau_Iso<<", cut from map = "<<IsoCut_PerEfficiency_PerBin[iEff][Name_Histo]<<", cut from histo = "<<IsoCut_PerBin[iEff]->GetBinContent(binForIsolation.at(0)+1,binForIsolation.at(1)+1,binForIsolation.at(2)+1)<<endl;
	  if(L1Tau_Iso<=IsoCut_PerBin[iEff]->GetBinContent(binForIsolation.at(0)+1,binForIsolation.at(1)+1,binForIsolation.at(2)+1)) pt_pass_efficiency_TH3[iEff]->Fill(L1Tau_IEt);
	  // if(L1Tau_Iso<=IsoCut_PerBin[iEff]->GetBinContent(binForIsolation.at(0)+1,binForIsolation.at(1)+1,binForIsolation.at(2)+1)) pt_pass_efficiency_TH3[iEff]->Fill(L1Tau_IEt);
	}

      // cout<<"IsoCut_PerEfficiency_PerBin[90][Name_Histo] = "<<IsoCut_PerEfficiency_PerBin[90][Name_Histo]<<endl;

      // if(L1Tau_Iso<=IsoCut_PerEfficiency_PerBin[90][Name_Histo]) pt_pass_90pc->Fill(L1Tau_IEt);
      // if(L1Tau_Iso<=IsoCut_PerEfficiency_PerBin[80][Name_Histo]) pt_pass_80pc->Fill(L1Tau_IEt);
      pt->Fill(L1Tau_IEt);
      eta->Fill(L1Tau_IEta);
      nTT->Fill(L1Tau_nTT);
    }
    std::cout<<"Event loop completed ! "<<endl;
    
    std::cout<<"Making the Efficiency graphs !  "<<endl;


  for(UInt_t iEff = 0 ; iEff < 101 ; ++iEff)
    {
      TString nameHisto = "pt_pass_efficiency_";
      ostringstream convert;
      convert << iEff;
      TString temp(convert.str());
      nameHisto += temp;
      
  //    std::cout<<nameHisto<<" <- "<<endl;
      TGraphAsymmErrors* temp_histo = new TGraphAsymmErrors(pt_pass_efficiency[iEff],pt,"cp");
      TGraphAsymmErrors* temp_histo_TH3 = new TGraphAsymmErrors(pt_pass_efficiency_TH3[iEff],pt,"cp");
      temp_histo->Write();
      temp_histo_TH3->Write();
    }

  for(UInt_t iEff = 0 ; iEff < 101 ; ++iEff)
    {
      TString nameHisto = "eta_pass_efficiency_";
      ostringstream convert;
      convert << iEff;
      TString temp(convert.str());
      nameHisto += temp;
      
//      std::cout<<nameHisto<<" <- "<<endl;
      TGraphAsymmErrors* temp_histo = new TGraphAsymmErrors(eta_pass_efficiency[iEff],eta,"cp");
      temp_histo->Write();
    }

  for(UInt_t iEff = 0 ; iEff < 101 ; ++iEff)
    {
      TString nameHisto = "nTT_pass_efficiency_";
      ostringstream convert;
      convert << iEff;
      TString temp(convert.str());
      nameHisto += temp;
      
 //     std::cout<<nameHisto<<" <- "<<endl;
      TGraphAsymmErrors* temp_histo = new TGraphAsymmErrors(nTT_pass_efficiency[iEff],nTT,"cp");
      temp_histo->Write();

    }

  for(UInt_t iEff = 0 ; iEff < 101 ; ++iEff)
    {
      IsoCut_PerBin[iEff]->Write();
    }

  pt->Write();
  eta->Write();
  nTT->Write();

  cout<<"NumberOfHistosWithLowStats/Tot = "<<NumberOfHistosWithLowStats<<"/"<<NbinsIEta*NbinsIEt*NbinsnTT<<endl;

  cout<<"Purging the file out now !! "<<endl;
  f_out->Purge();
  cout<<"Writing the file out now !! "<<endl;
  f_out->Write();
  cout<<"Closing the output file : "<<f_out->GetName()<<endl;
  f_out->Close();

}
