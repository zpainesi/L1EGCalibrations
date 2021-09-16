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
#include "HybridGBRForest.h"
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
#include "RooCBExp.h"
#include "RooCBFast.h"
#include "RooGaussianFast.h"

 
using namespace RooFit;
  

void eregtrainingTest() 
{

    //build vectors with list of input variables
    std::vector<std::string> *varslist = new std::vector<std::string>;
    varslist->push_back("nVtx");
    varslist->push_back("scRawEnergy");
    varslist->push_back("scEta");
    varslist->push_back("scPhi");
    varslist->push_back("scEtaWidth");
    varslist->push_back("scPhiWidth");
    varslist->push_back("scSeedR9");
    varslist->push_back("scSeedRawEnergy/scRawEnergy");
    varslist->push_back("scSeedEmax/scRawEnergy");
    varslist->push_back("scSeedE2nd/scRawEnergy");
    varslist->push_back("scSeedLeftRightAsym");
    varslist->push_back("scSeedTopBottomAsym");
    varslist->push_back("scSeedSigmaIetaIeta");
    varslist->push_back("scSeedSigmaIetaIphi");
    varslist->push_back("scSeedSigmaIphiIphi");
    varslist->push_back("N_ECALClusters");
    //varslist->push_back("clusterMaxDR");
    //varslist->push_back("clusterMaxDRDPhi");
    //varslist->push_back("clusterMaxDRDEta");
    //varslist->push_back("clusterMaxDRRawEnergy/scRawEnergy");
    //varslist->push_back("clusterRawEnergy[0]/scRawEnergy");
    //varslist->push_back("clusterRawEnergy[1]/scRawEnergy");
    //varslist->push_back("clusterRawEnergy[2]/scRawEnergy");
    //varslist->push_back("clusterDPhiToSeed[0]");
    //varslist->push_back("clusterDPhiToSeed[1]");
    //varslist->push_back("clusterDPhiToSeed[2]");
    //varslist->push_back("clusterDEtaToSeed[0]");
    //varslist->push_back("clusterDEtaToSeed[1]");
    //varslist->push_back("clusterDEtaToSeed[2]");
    //varslist->push_back("scSeedCryEta");
    //varslist->push_back("scSeedCryPhi");
    //varslist->push_back("scSeedCryIeta");
    //varslist->push_back("scSeedCryIphi");

    //create RooRealVars for each input variable
    RooArgList vars;
    for (unsigned int ivar=0; ivar<varslist->size(); ++ivar) 
    {
        RooRealVar *var = new RooRealVar(TString::Format("var_%i",ivar),varslist->at(ivar).c_str(),0.);
        vars.addOwned(*var);
    }

    //make list of input variable RooRealVars
    RooArgList condvars(vars);

    //create RooRealVar for target
    RooRealVar *tgtvar = new RooRealVar("tgtvar","genEnergy/scRawEnergy",1.);


    //add target to full list
    vars.addOwned(*tgtvar);

    //RooRealVar for event weight 
    RooRealVar weightvar("weightvar","",1.);

    TChain *tree = new TChain("gedGsfElectronTree/RegressionTree");
    tree->Add("/data_CMS/cms/sauvan/Regression/Ntuples/DoubleElectron_FlatPt-1To300/regression_ntuple_RunIISpring15DR74-AsymptFlat0to50bx25RawReco_MCRUN2_74_V9A-v2_2015_07_15/regression_ntuple_746_1.root");

    //training selection cut
    TCut selcut("(isMatched==1)&&(scIsEB==1)");

    //weightvar title used for per-event weights and selection cuts
    weightvar.SetTitle(selcut);
    //create RooDataSet from TChain
    RooDataSet *hdata = RooTreeConvert::CreateDataSet("hdata",tree,vars,weightvar);   


    //RooRealVars corresponding to regressed parameters (sigma, mean, left tail parameter, right tail parameter)
    RooRealVar sigwidthtvar("sigwidthtvar","",0.01);
    sigwidthtvar.setConstant(false);

    RooRealVar sigmeantvar("sigmeantvar","",1.);
    sigmeantvar.setConstant(false); 

    RooRealVar signvar("signvar","",3.);
    signvar.setConstant(false);       

    RooRealVar sign2var("sign2var","",3.);
    sign2var.setConstant(false);     

    //define non-parametric functions for each regressed parameter
    RooGBRFunctionFlex *sigwidthtfunc = new RooGBRFunctionFlex("sigwidthtfunc","");
    RooGBRFunctionFlex *sigmeantfunc = new RooGBRFunctionFlex("sigmeantfunc","");
    RooGBRFunctionFlex *signfunc = new RooGBRFunctionFlex("signfunc","");
    RooGBRFunctionFlex *sign2func = new RooGBRFunctionFlex("sign2func","");

    //define mapping of input variables to non-parametric functions (in this case trivial since all 4 functions depend on the same inputs, but this is not a requirement)
    RooGBRTargetFlex *sigwidtht = new RooGBRTargetFlex("sigwidtht","",*sigwidthtfunc,sigwidthtvar,condvars);  
    RooGBRTargetFlex *sigmeant = new RooGBRTargetFlex("sigmeant","",*sigmeantfunc,sigmeantvar,condvars);  
    RooGBRTargetFlex *signt = new RooGBRTargetFlex("signt","",*signfunc,signvar,condvars);  
    RooGBRTargetFlex *sign2t = new RooGBRTargetFlex("sign2t","",*sign2func,sign2var,condvars);  

    //define list of mapped functions to regress
    RooArgList tgts;
    tgts.add(*sigwidtht);
    tgts.add(*sigmeant);
    tgts.add(*signt);
    tgts.add(*sign2t);  

    //define transformations corresponding to parameter bounds for non-parametric outputs  
    RooRealConstraint sigwidthlim("sigwidthlim","",*sigwidtht,0.0002,0.5);
    RooRealConstraint sigmeanlim("sigmeanlim","",*sigmeant,0.2,2.0);
    RooRealConstraint signlim("signlim","",*signt,1.01,5000.); 
    RooRealConstraint sign2lim("sign2lim","",*sign2t,1.01,5000.); 

    //define pdf, which depends on transformed outputs (and is intended to be treated as a conditional pdf over the
    //regression inputs in this case)
    //The actual pdf below is a double crystal ball, with crossover points alpha_1 and alpha_2 set constant, but all other
    //parameters regressed
    RooDoubleCBFast sigpdf("sigpdf","",*tgtvar,sigmeanlim,sigwidthlim,RooConst(2.),signlim,RooConst(1.),sign2lim);

    //dummy variable
    RooConstVar etermconst("etermconst","",0.);  

    //dummy variable
    RooRealVar r("r","",1.);
    r.setConstant();

    //define list of pdfs
    std::vector<RooAbsReal*> vpdf;
    vpdf.push_back(&sigpdf);  

    //define list of training datasets
    std::vector<RooAbsData*> vdata;
    vdata.push_back(hdata);     

    //define minimum event weight per tree node
    double minweight = 200;
    std::vector<double> minweights;
    minweights.push_back(minweight);

    //temp output file
    TFile *fres = new TFile("fres.root","RECREATE");

    //run training
    if (1) {
        RooHybridBDTAutoPdf bdtpdfdiff("bdtpdfdiff","",tgts,etermconst,r,vdata,vpdf);
        bdtpdfdiff.SetMinCutSignificance(5.);
        //bdtpdfdiff.SetPrescaleInit(100);
        bdtpdfdiff.SetShrinkage(0.1);
        bdtpdfdiff.SetMinWeights(minweights);
        bdtpdfdiff.SetMaxNodes(750);
        bdtpdfdiff.TrainForest(1e6);   
    }

    //create workspace and output to file
    RooWorkspace *wereg = new RooWorkspace("wereg");
    wereg->import(sigpdf);

    wereg->writeToFile("wereg_ele_eb.root");    

    return;


}
