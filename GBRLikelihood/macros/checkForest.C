//#include "RooRealVar.h"
//#include "RooAbsPdf.h"
//#include "RooExponential.h"
//#include "RooGaussian.h"
//#include "RooPlot.h"
//#include "TCanvas.h"
//#include "RooConstVar.h"
//#include "RooDataSet.h"
//#include "HiggsAnalysis/GBRLikelihood/interface/RooHybridBDTAutoPdf.h"
//#include "RooFormulaVar.h"
//#include "RooProdPdf.h"
//#include "RooUniform.h"
//#include "TRandom.h"
//#include "TGraph.h"
//#include "RooAddPdf.h"
//#include "RooNDKeysPdf.h"
//#include "RooExtendPdf.h"
//#include "RooMinimizer.h"
#include "TFile.h"
#include "TNtuple.h"
//#include "HiggsAnalysis/GBRLikelihood/interface/HybridGBRForestFlex.h"
//#include "RooProduct.h"
//#include "RooChebychev.h"
//#include "RooBernstein.h"
//#include "RooPolynomial.h"
//#include "RooGenericPdf.h"
////#include "HZZ2L2QRooPdfs.h"
//#include "HiggsAnalysis/GBRLikelihood/interface/RooDoubleCBFast.h"
//#include "RooArgSet.h"
//#include "RooArgList.h"
//#include "RooCBShape.h"
//#include "RooWorkspace.h"
#include "TH1D.h"
#include "TChain.h"
#include "TCut.h"
#include "TLine.h"
#include "TLegend.h"
//#include "RooRandom.h"
//#include "RooAddition.h"
#include "TSystem.h"
#include "TTreeFormula.h"
//#include "RooLinearVar.h"
#include "GBRMath.h"

//#include "Cintex/Cintex.h"

//#include "HybridGBRForestFlex.h"

#include <iostream>
#include <vector>
#include <string>

#include <math.h>


#ifndef __CINT__
  #include "CondFormats/EgammaObjects/interface/GBRForestD.h"
#endif

//using namespace RooFit;

void checkForest() 
{
    //read workspace from training
    TFile *f = TFile::Open("wereg_test.root"); 
    GBRForestD* forest = (GBRForestD*)f->Get("regmean");

    TFile* finput = TFile::Open("/data_CMS/cms/sauvan/Regression/Ntuples/DoubleElectron_FlatPt-1To300/regression_ntuple_RunIISpring15DR74-AsymptFlat0to50bx25RawReco_MCRUN2_74_V9A-v2_2015_07_15/regression_ntuple_746_1.root");
    TTree* intree = (TTree*)finput->Get("gedGsfElectronTree/RegressionTree");

    std::vector<std::string> vars;
    vars.push_back("nVtx");
    vars.push_back("scRawEnergy");
    vars.push_back("scEta");
    vars.push_back("scPhi");
    vars.push_back("scEtaWidth");
    vars.push_back("scPhiWidth");
    vars.push_back("scSeedR9");
    vars.push_back("scSeedRawEnergy/scRawEnergy");
    vars.push_back("scSeedEmax/scRawEnergy");
    vars.push_back("scSeedE2nd/scRawEnergy");
    vars.push_back("scSeedLeftRightAsym");
    vars.push_back("scSeedTopBottomAsym");
    vars.push_back("scSeedSigmaIetaIeta");
    vars.push_back("scSeedSigmaIetaIphi");
    vars.push_back("scSeedSigmaIphiIphi");
    vars.push_back("N_ECALClusters");
    //vars.push_back("clusterMaxDR");
    //vars.push_back("clusterMaxDRDPhi");
    //vars.push_back("clusterMaxDRDEta");
    //vars.push_back("clusterMaxDRRawEnergy/scRawEnergy");
    //vars.push_back("clusterRawEnergy[0]/scRawEnergy");
    //vars.push_back("clusterRawEnergy[1]/scRawEnergy");
    //vars.push_back("clusterRawEnergy[2]/scRawEnergy");
    //vars.push_back("clusterDPhiToSeed[0]");
    //vars.push_back("clusterDPhiToSeed[1]");
    //vars.push_back("clusterDPhiToSeed[2]");
    //vars.push_back("clusterDEtaToSeed[0]");
    //vars.push_back("clusterDEtaToSeed[1]");
    //vars.push_back("clusterDEtaToSeed[2]");
    //vars.push_back("scSeedCryEta");
    //vars.push_back("scSeedCryPhi");
    //vars.push_back("scSeedCryIeta");
    //vars.push_back("scSeedCryIphi");


    int nvars = vars.size();

    //initialize TTreeFormulas to read variables from TTree
    std::vector<TTreeFormula*> inputforms;
    for (std::vector<std::string>::const_iterator it = vars.begin(); it != vars.end(); ++it) 
    {
        inputforms.push_back(new TTreeFormula(it->c_str(),it->c_str(),intree));
    }

    Float_t target = 0.;
    Float_t *vals = new Float_t[nvars];

    double low = 0.2;
    double high = 2.;
    double scale = 0.5*(high-low);
    double offset = low + 0.5*(high-low);


    for (Long64_t iev=0; iev<intree->GetEntries(); ++iev) {
        //if (iev%100000==0) printf("%i\n",int(iev));
        intree->LoadTree(iev);

        for (int i=0; i<nvars; ++i) {
            vals[i] = inputforms[i]->EvalInstance();
        }

        target = forest->GetResponse(vals);
        //double actualResponse = offset + scale*sin(target);
        double actualResponse = offset + scale*vdt::fast_sin(target);
        std::cout<<"Initial response = "<<offset + scale*sin(forest->InitialResponse())<<", Response = "<<actualResponse<<"\n";

    }


    //clear TTreeFormulas
    for (std::vector<TTreeFormula*>::const_iterator it = inputforms.begin(); 
            it != inputforms.end(); ++it) {
        delete *it;
    }

    delete[] vals;

    f->Close();
    finput->Close();


}
