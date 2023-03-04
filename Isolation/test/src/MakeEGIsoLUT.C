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

void MakeEGIsoLUT(string option, bool isOption, bool includeCompression)
{
  TFile* fLUTS = new TFile("HistgramFile_step2.root");
  TH3F* LUT ;
 
  fLUTS->cd("Step2Histos");
  TString Option="LUT_Progression_"+option;
  if(isOption) LUT = (TH3F*) gDirectory->Get(Option);
  else { 
    int Option_=std::atoi(option.c_str());
    TString Option="LUT_WP"+to_string(Option_);
    LUT = (TH3F*) gDirectory->Get(Option);
  }
  if(!LUT) std::cout<<"leave";
  cout << "LUT name: " << LUT->GetName() << endl;
  cout << "Max LUT CuT: " << LUT->GetMaximum() << endl;

  //  const Float_t hardcodedIetaBins[] = {0, 6, 12, 18, 31}; //for tau
  const Float_t hardcodedIetaBins[] = {0, 5, 6, 9, 10, 12, 13, 14, 17, 18, 19, 20, 23, 24, 25, 26, 31};
    //  const Float_t hardcodedIetBins[]  = {0, 15, 18, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47, 50, 53, 56, 59, 62, 65, 69, 73, 77, 82, 88, 95, 105, 120, 157, 255}; //for Tau
  const Float_t hardcodedIetBins[]  = {0, 18, 20, 22, 28, 32, 37, 42, 52, 63, 73, 81, 87, 91, 111, 151, 255};
  const Float_t hardcodednTTBins[] = {0, 6, 11, 16, 21, 26, 31, 36, 41, 46, 51, 56, 61, 66, 71, 76, 81, 86, 91, 96, 101, 106, 111, 116, 121, 126, 131, 136, 141, 146, 161, 161, 255};


  int neta = sizeof(hardcodedIetaBins) / sizeof(int) -1 ;
  int nEt  = sizeof(hardcodedIetBins) / sizeof(int)-1;
  int nnTT = sizeof(hardcodednTTBins) / sizeof(int)-1;

  TH3F* binning = new TH3F ("binning", "binning", neta, hardcodedIetaBins, nEt, hardcodedIetBins, nnTT, hardcodednTTBins);

  float maxeta = binning->GetXaxis()->GetBinLowEdge(LUT->GetNbinsX()+1);
  float maxEt  = binning->GetYaxis()->GetBinLowEdge(LUT->GetNbinsY()+1);
  float maxnTT = binning->GetZaxis()->GetBinLowEdge(LUT->GetNbinsZ()+1);

  cout << maxeta << " " << maxEt << " " << maxnTT << endl;
  cout << "Output cmpr bins: " << neta << " " << nEt << " " << nnTT << endl;
  

  TString outFile ;
  if(isOption && includeCompression) outFile = "EG_Iso_LUT_Option_" + option + ".txt";
  else if(!isOption && includeCompression) outFile = "EG_Iso_LUT_Flat_WP_" + option + ".txt";
  else if(isOption && !includeCompression) outFile = "EG_Iso_LUT_Option_" + option + ".txt";
  else if(!isOption && !includeCompression) outFile = "EG_Iso_LUT_Flat_WP_" + option + ".txt";

  std::ofstream LUTfile (outFile.Data());

  ////////////////////////////////////////////////////////////////////////////
  //// print compression block on top of the LUT /////////////////////////////

  int etaBits = 5;
  int etBits  = 8;
  int nTTBits = 10;

  int TotalIDX = 0;

  ////////////////////////////////////////////////////////////////////////////
  ///// print header
  ////////////////////////////////////////////////////////////////////////////
  //# comments...
  //#<header> V1 12 7 </header>
  // eta compr bits: 2
  // Et compr bits:  3
  // nTT compr bits: 3
  // input bits: 2**6 + 2**13 + 2**10 + 2**(2+3+3) = 271k ==> stays on 14 bits
    
  int cmprEta = 4;
  int cmprEt = 4;
  int cmprnTT = 5;

  int totInBits  = 13;//23 ;
  int totOutBits = 9;//12 ;  // i.e. cut max 127, but it is largely overestimated: max in LUT is 27 --> 5 bits are enough
  
  if(isOption) LUTfile << "# Tau isolation LUT with Option = " << option << "" << endl;
  else LUTfile << "# Tau isolation LUT with WP = " << option << "%" << endl;
  if(includeCompression) LUTfile << "# on top of the table there is the in/out compression scheme, followed by the LUT" << endl;
  LUTfile << "# iso LUT structure is ieta --> iEt -->  nTT" << endl;
  if(includeCompression) LUTfile << "# Input bits: ieta: " << etaBits << " iEt " << etBits << " nTT: " << nTTBits << endl;
  LUTfile << "# Compr bits: ieta: " << cmprEta << " iEt: " << cmprEt << " nTT: " << cmprnTT  << endl;

  LUTfile << "#<header> V10.0 " << totInBits << " " << totOutBits << " </header>"<<endl;
  LUTfile << "\n";
  std::cout<<"hi:"<<includeCompression<<std::endl;
  if(includeCompression)
    {

      ////////// iEt compr
      for (int i = 0; i < pow(2, etBits); i++, TotalIDX++)
	{
	  int BinEt = binning->GetYaxis()->FindBin (i);
	  if (BinEt < 1) BinEt = 1;
	  if (BinEt > nEt) BinEt = nEt;
	  LUTfile << TotalIDX << " " << BinEt-1;
	  if (i == 0) LUTfile << " # begin of ET compr block";
	  LUTfile << endl;
	}

      ////////// inTT compr
      for (int i = 0; i < pow(2, nTTBits); i++, TotalIDX++)
	{
	  int BinnTT = binning->GetZaxis()->FindBin (i);
	  if (BinnTT < 1) BinnTT = 1;
	  if (BinnTT > nnTT) BinnTT = nnTT;
	  LUTfile << TotalIDX << " " << BinnTT-1;
	  if (i == 0) LUTfile << " # begin of nTT compr block";
	  LUTfile << endl;
	}

      ////////// iEta compr
      for (int i = 0; i < pow(2, etaBits); i++, TotalIDX++)
	{
	  int BinEta = binning->GetXaxis()->FindBin (i);
	  if (BinEta < 1) BinEta = 1;
	  if (BinEta > neta) BinEta = neta;
	  LUTfile << TotalIDX << " " << BinEta-1;
	  if (i == 0) LUTfile << " # begin of eta compr block";
	  LUTfile << endl;
	}

    }

  ////////////////////////////////////////////////////////////////////////////
  //// print isolation LUT ///////////////////////////////////////////////////

  // structure: eta -- Et -- nTT
  int onlyIsoAddr = 0;
  for (int ieta = 0; ieta < pow(2, cmprEta); ieta++)
    {
      for (int iEt = 0; iEt < pow(2, cmprEt); iEt++)
	{
	  for (int inTT = 0; inTT < pow(2, cmprnTT); inTT++)
	    {
	      int binEt = iEt;
	      int bineta = ieta;
	      int binnTT = inTT;
                
	      if (binEt >= nEt)   binEt  = nEt -1;   // these are never used -- values that are outside the compression
	      if (bineta >= neta) bineta = neta -1;
	      if (binnTT >= nnTT) binnTT = nnTT -1;

	      int thr = LUT->GetBinContent(bineta+1,binEt+1,binnTT+1)+1;
	      // if(ieta==3 && ((iEt==2) || (iEt==3)) && inTT>=11) thr = LUT->GetBinContent(bineta,binEt+1,binnTT+1)+1;
	      if(thr>=pow(2,totOutBits)-1) thr=pow(2,totOutBits)-1;
	      LUTfile << TotalIDX << " " << thr;
	      if (iEt == 0 && ieta == 0 && inTT == 0) LUTfile << " # start of isolation LUT -- ieta : iEt : nTT = " << ieta << " : " << iEt << " : " << inTT; 
	      else LUTfile << " # ieta : iEt : inTT = " << ieta << " : " << iEt << " : " << inTT; 
	      LUTfile << endl;
	      cout << onlyIsoAddr << " # ieta : iEt : inTT = " << ieta << " : " << iEt << " : " << inTT <<  " --> thr: " << thr << endl; 
	      
	      TotalIDX++;
	      onlyIsoAddr++;
            }
        }
    }

  int totIDXexpected = pow(2,etaBits) + pow(2,etBits) + pow (2, nTTBits) + pow (2, cmprEta+cmprEt+cmprnTT);
  cout << "TOTAL IDX " << TotalIDX << " expected: " << totIDXexpected << endl;
  cout << "Max allowed by input LUT bits: " << pow(2, totInBits) << endl;
  if (TotalIDX >= pow(2, totInBits))
    {
      cout << " **** WARNING!!! too few input bits " << totInBits << " , please increase value in header" << endl;
    } 
  

}
