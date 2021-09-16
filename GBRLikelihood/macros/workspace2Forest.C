#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooExponential.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "RooConstVar.h"
#include "RooDataSet.h"
#include "RooHybridBDTAutoPdf.h"
#include "RooFormulaVar.h"
#include "RooProdPdf.h"
#include "RooUniform.h"
#include "TRandom.h"
#include "TGraph.h"
#include "RooAddPdf.h"
#include "RooNDKeysPdf.h"
#include "RooExtendPdf.h"
#include "RooMinimizer.h"
#include "TFile.h"
#include "TNtuple.h"
#include "HybridGBRForestFlex.h"
#include "RooProduct.h"
#include "RooChebychev.h"
#include "RooBernstein.h"
#include "RooPolynomial.h"
#include "RooGenericPdf.h"
//#include "HZZ2L2QRooPdfs.h"
#include "RooDoubleCBFast.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooCBShape.h"
#include "RooWorkspace.h"
#include "TH1D.h"
#include "TChain.h"
#include "TCut.h"
#include "TLine.h"
#include "TLegend.h"
#include "RooRandom.h"
#include "RooAddition.h"
#include "TSystem.h"
#include "RooLinearVar.h"

#include <iostream>

//#include "Cintex/Cintex.h"

//#include "HybridGBRForestFlex.h"

#ifndef __CINT__
  #include "CondFormats/EgammaObjects/interface/GBRForestD.h"
#endif

using namespace RooFit;

void workspace2Forest() {

    //gSystem->Load("../../../../lib/slc6_amd64_gcc491/libHiggsAnalysisGBRLikelihood.so");
  
  //output dir
//   TString dirname = "/data/bendavid/eregexampletest/eregexampletest_test/"; 
//   gSystem->mkdir(dirname,true);
//   gSystem->cd(dirname);    
  
  //read workspace from training
  //TFile *fws = TFile::Open("/home/llr/cms/sauvan/Regression/RegressionTraining/GBR_Clustering_746_bx25_Electrons_GedGsfElectron_RegPs1_results.root"); 
  TFile *fws = TFile::Open("wereg_ele_eb.root"); 
  RooWorkspace *ws = (RooWorkspace*)fws->Get("wereg");
  
  //read variables from workspace
  RooGBRTargetFlex *meantgt = static_cast<RooGBRTargetFlex*>(ws->arg("sigmeant"));  
  RooRealVar *tgtvar = ws->var("tgtvar");
  
  
  RooArgList vars;
  vars.add(meantgt->FuncVars());
  vars.add(*tgtvar);
  
  //ROOT::Cintex::Cintex::Enable();
  
  std::cout<<"Initial response = "<<meantgt->Forest()->InitialResponse()<<"\n";

  TFile *fout = new TFile("wereg_test.root","RECREATE");
  GBRForestD *forestout = new GBRForestD(*meantgt->Forest());
  //HybridGBRForestFlex* forestout = new HybridGBRForestFlex(*meantgt->Forest());
  fout->WriteObject(forestout,"regmean");
  fout->Close();
  
  
}
