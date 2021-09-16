#include <TTree.h>
#include <TString.h>
#include <TFile.h>
#include <TChain.h>
#include <TH1F.h>

//#include "Helpers.C"

using namespace std;

#ifdef __MAKECINT__
#pragma link C++ class vector<float>+;
#endif


void makeRegressionTree(){ //vector<TString> input_names, TString output_name){

/*
  TChain* intree = new TChain("produceNtuple/eIDSimpleTree");
  for(unsigned int i=0; i<input_names.size();i++)
    intree->Add(input_names[i]);
*/

  TFile *f=TFile::Open("/grid_mnt/t3storage3/athachay/l1egamma/emulationstuff/CMSSW_7_6_0/src/MergedFile.root");

  TTree *intree=(TTree*)f->Get("Ntuplizer/TagAndProbe");


  TString output_file("./test_mean.root");
  TFile* f_new = TFile::Open(output_file,"RECREATE");

  TTree* tree = new TTree("eIDSimpleTree","eIDSimpleTree");

  //Branches new tree

  ULong64_t _event;
  Int_t _lumi;
  Int_t _run;  
  int _ieta;
  int _E;
  int _shape;
  float _target;
  int _HoverERatio;
  int _nTT;
 // int _IsoSum;
  int _Run2IdLevel;
 
  tree->Branch("Run",&_run,"_run/I");
  tree->Branch("Event",&_event,"_event/l");
  tree->Branch("Lumi",&_lumi,"_lumi/I");
  tree->Branch("ieta",&_ieta,"_ieta/I");
  tree->Branch("E",&_E,"_E/I");
  tree->Branch("shape",&_shape,"_shape/I");
  tree->Branch("target",&_target,"_target/F");
  tree->Branch("HoverERatio",&_HoverERatio,"_HoverERatio/I");
  tree->Branch("nTT",&_nTT,"_nTT/I");
 // tree->Branch("IsoSum",&_IsoSum,"_IsoSum/I");
  tree->Branch("Run2IdLevel",&_Run2IdLevel,"_Run2IdLevel/I");

 
  //Branches in intree
  ULong64_t _nEvent;
  Int_t _nRun;
  Int_t _nLumi;
  float _ele_sclEt;//[10];
  int _ele_L1Stage2_emul_ieta;//[10];
  //int _ele_L1Stage2_emul_hwPt[10];
  int _ele_L1Stage2_emul_rawEt;//[10];
  int _ele_L1Stage2_emul_shape;//[10];
  int _ele_L1Stage2_emul_hOverERatio;//[10];
  int _ele_L1Stage2_emul_nTT;//[10];
  int _ele_L1Stage2_emul_hwIsoSum;//[10];
 // bool _ele_VetoIdDecisions[10];
  int _ele_LooseIdDecisions;//[10];
 // bool _ele_MediumIdDecisions[10];
 // bool _ele_TightIdDecisions[10];
 
  intree->SetBranchAddress("EventNumber",&_nEvent);
  intree->SetBranchAddress("RunNumber",&_nRun);
  intree->SetBranchAddress("lumi",&_nLumi);
  intree->SetBranchAddress("eleProbeSclEt",&_ele_sclEt);
  intree->SetBranchAddress("l1tEmuTowerIEta",&_ele_L1Stage2_emul_ieta);
  //intree->SetBranchAddress("ele_L1Stage2_emul_hwPt",&_ele_L1Stage2_emul_hwPt);
  intree->SetBranchAddress("l1tEmuRawEt",&_ele_L1Stage2_emul_rawEt);
  intree->SetBranchAddress("shape",&_ele_L1Stage2_emul_shape);
  intree->SetBranchAddress("TowerHoE",&_ele_L1Stage2_emul_hOverERatio);
  intree->SetBranchAddress("l1tEmuNTT",&_ele_L1Stage2_emul_nTT);
  //intree->SetBranchAddress("ele_L1Stage2_emul_hwIsoSum",&_ele_L1Stage2_emul_hwIsoSum);
  //intree->SetBranchAddress("ele_VetoIdDecisions",&_ele_VetoIdDecisions);
  intree->SetBranchAddress("isProbeLoose",&_ele_LooseIdDecisions);
  //intree->SetBranchAddress("ele_MediumIdDecisions",&_ele_MediumIdDecisions);
  //intree->SetBranchAddress("ele_TightIdDecisions",&_ele_TightIdDecisions);



  Long64_t nentries = intree->GetEntries();

  for(int i=0;i<nentries;i++){

    if(i%10000==0){
       std::cout<<"i="<<i<<endl;
    }

    intree->GetEntry(i);
    
  //  for(unsigned int i_eg=0;i_eg<10;i_eg++){
      
      _run = 0;
      _event = 0;
      _lumi = 0;
      _ieta = 0;
      _E = 0;
      _shape = 0;
      _target = 0;      
      _HoverERatio = -1;      
      _nTT = -1;      
   //   _IsoSum = -1; 
      _Run2IdLevel = -1;      
     
      if( _ele_L1Stage2_emul_rawEt>0 ){
	_run = _nRun;
	_event = _nEvent;
	_lumi = _nLumi;

	_ieta = _ele_L1Stage2_emul_ieta;
	_E = _ele_L1Stage2_emul_rawEt;
	_shape = _ele_L1Stage2_emul_shape;
	_target = _ele_sclEt/(0.5*_E);
	_HoverERatio = _ele_L1Stage2_emul_hOverERatio;
	_nTT = _ele_L1Stage2_emul_nTT;
//	_IsoSum = _ele_L1Stage2_emul_hwIsoSum[i_eg];

        if(_ele_LooseIdDecisions){
          _Run2IdLevel = 1;
        }
        //else if(_ele_MediumIdDecisions[i_eg])
        //  _Run2IdLevel = 2;
//      else if(_ele_LooseIdDecisions[i_eg])
//        _Run2IdLevel = 1;
//      else if(_ele_VetoIdDecisions[i_eg])
//        _Run2IdLevel = 0;


	tree->Fill();

      }
      
   // }      
	
  }

  tree->Write();
  f_new->Close();


  return;
}

/*

void calibrtree_maker_etastrip_excluded(vector<TString> input_names, TString output_name){

  TChain* intree = new TChain("produceNtuple/eIDSimpleTree");
  for(unsigned int i=0; i<input_names.size();i++)
    intree->Add(input_names[i]);


  TString output_file(output_name);
  TFile* f_new = TFile::Open(output_file,"RECREATE");

  TTree* tree = new TTree("eIDSimpleTree","eIDSimpleTree");

  //Branches new tree

  int _event;
  int _lumi;
  int _run;  
  int _ieta;
  int _E;
  int _shape;
  double _target;
  int _HoverERatio;
  int _nTT;
  int _IsoSum;
  int _Run2IdLevel;
 
  tree->Branch("Run",&_run,"_run/I");
  tree->Branch("Event",&_event,"_event/I");
  tree->Branch("Lumi",&_lumi,"_lumi/I");
  tree->Branch("ieta",&_ieta,"_ieta/I");
  tree->Branch("E",&_E,"_E/I");
  tree->Branch("shape",&_shape,"_shape/I");
  tree->Branch("target",&_target,"_target/D");
  tree->Branch("HoverERatio",&_HoverERatio,"_HoverERatio/I");
  tree->Branch("nTT",&_nTT,"_nTT/I");
  tree->Branch("IsoSum",&_IsoSum,"_IsoSum/I");
  tree->Branch("Run2IdLevel",&_Run2IdLevel,"_Run2IdLevel/I");

 
  //Branches in intree
  int _nEvent;
  int _nRun;
  int _nLumi;
  double _ele_sclEt;//[10];
  int _ele_L1Stage2_emul_ieta;//[10];
  int _ele_L1Stage2_emul_hwPt;//[10];
  int _ele_L1Stage2_emul_shape;//[10];
  int _ele_L1Stage2_emul_hOverERatio;//[10];
  int _ele_L1Stage2_emul_nTT;//[10];
  int _ele_L1Stage2_emul_hwIsoSum;//[10];
  bool _ele_VetoIdDecisions;//[10];
  bool _ele_LooseIdDecisions;//[10];
  bool _ele_MediumIdDecisions;//[10];
  bool _ele_TightIdDecisions;//[10];
 
  intree->SetBranchAddress("nEvent",&_nEvent);
  intree->SetBranchAddress("nRun",&_nRun);
  intree->SetBranchAddress("nLumi",&_nLumi);
  intree->SetBranchAddress("ele_sclEt",_ele_sclEt);
  intree->SetBranchAddress("ele_L1Stage2_emul_ieta",_ele_L1Stage2_emul_ieta);
  intree->SetBranchAddress("ele_L1Stage2_emul_hwPt",&_ele_L1Stage2_emul_hwPt);
  intree->SetBranchAddress("ele_L1Stage2_emul_shape",&_ele_L1Stage2_emul_shape);
  intree->SetBranchAddress("ele_L1Stage2_emul_hOverERatio",&_ele_L1Stage2_emul_hOverERatio);
  intree->SetBranchAddress("ele_L1Stage2_emul_nTT",&_ele_L1Stage2_emul_nTT);
  intree->SetBranchAddress("ele_L1Stage2_emul_hwIsoSum",&_ele_L1Stage2_emul_hwIsoSum);
  intree->SetBranchAddress("ele_VetoIdDecisions",&_ele_VetoIdDecisions);
  intree->SetBranchAddress("ele_LooseIdDecisions",&_ele_LooseIdDecisions);
  intree->SetBranchAddress("ele_MediumIdDecisions",&_ele_MediumIdDecisions);
  intree->SetBranchAddress("ele_TightIdDecisions",&_ele_TightIdDecisions);


  Long64_t nentries = intree->GetEntries();

  for(int i=0;i<nentries;i++){
    
    intree->GetEntry(i);
    
    //for(unsigned int i_eg=0;i_eg<10;i_eg++){
      
      _run = 0;
      _event = 0;
      _lumi = 0;
      _ieta = 0;
      _E = 0;
      _shape = 0;
      _target = 0;      
      _HoverERatio = -1;      
      _nTT = -1;      
      _IsoSum = -1; 
      _Run2IdLevel = -1;    

      if( _ele_L1Stage2_emul_ieta>=13 && _ele_L1Stage2_emul_ieta<=17)
	continue;
     
      if( _ele_L1Stage2_emul_hwPt>0 ){
	_run = _nRun;
	_event = _nEvent;
	_lumi = _nLumi;

	_ieta = _ele_L1Stage2_emul_ieta;
	_E = _ele_L1Stage2_emul_hwPt;
	_shape = _ele_L1Stage2_emul_shape;
	_target = _ele_sclEt/(0.5*_E);
	_HoverERatio = _ele_L1Stage2_emul_hOverERatio;
	_nTT = _ele_L1Stage2_emul_nTT;
	_IsoSum = _ele_L1Stage2_emul_hwIsoSum;

	if(_ele_TightIdDecisions[i_eg])
	  _Run2IdLevel = 3;
	else if(_ele_MediumIdDecisions[i_eg])
	  _Run2IdLevel = 2;
	else if(_ele_LooseIdDecisions[i_eg])
	  _Run2IdLevel = 1;
	else if(_ele_VetoIdDecisions[i_eg])
	  _Run2IdLevel = 0;

	tree->Fill();

      }
      
   // }      
	
  }

  tree->Write();
  f_new->Close();


  return;
}



void histoshape_maker(vector<TString> input_names, TString output_name){

  TChain* intree = new TChain("produceNtuple/eIDSimpleTree");
  for(unsigned int i=0; i<input_names.size();i++)
    intree->Add(input_names[i]);


  TString output_file(output_name);
  TFile* f_new = TFile::Open(output_file,"RECREATE");

  vector<TH1F*> histos_shape;
  for(unsigned int i=0; i<27;i++){
    TH1F* h=new TH1F(Form("h_ieta%i",i),Form("h_ieta%i",i),128,-0.5,127.5);
    histos_shape.push_back(h);
  }


  //Branches in intree
  int _nEvent;
  int _nRun;
  int _nLumi;
  double _ele_sclEt[10];
  int _ele_L1Stage2_emul_ieta[10];
  int _ele_L1Stage2_emul_hwPt[10];
  int _ele_L1Stage2_emul_shape[10];
  int _ele_L1Stage2_emul_hOverERatio[10];
  int _ele_L1Stage2_emul_nTT[10];
  int _ele_L1Stage2_emul_hwIsoSum[10];
  bool _ele_VetoIdDecisions[10];
  bool _ele_LooseIdDecisions[10];
  bool _ele_MediumIdDecisions[10];
  bool _ele_TightIdDecisions[10];
 
  intree->SetBranchAddress("nEvent",&_nEvent);
  intree->SetBranchAddress("nRun",&_nRun);
  intree->SetBranchAddress("nLumi",&_nLumi);
  intree->SetBranchAddress("ele_sclEt",_ele_sclEt);
  intree->SetBranchAddress("ele_L1Stage2_emul_ieta",_ele_L1Stage2_emul_ieta);
  intree->SetBranchAddress("ele_L1Stage2_emul_hwPt",&_ele_L1Stage2_emul_hwPt);
  intree->SetBranchAddress("ele_L1Stage2_emul_shape",&_ele_L1Stage2_emul_shape);
  intree->SetBranchAddress("ele_L1Stage2_emul_hOverERatio",&_ele_L1Stage2_emul_hOverERatio);
  intree->SetBranchAddress("ele_L1Stage2_emul_nTT",&_ele_L1Stage2_emul_nTT);
  intree->SetBranchAddress("ele_L1Stage2_emul_hwIsoSum",&_ele_L1Stage2_emul_hwIsoSum);
  intree->SetBranchAddress("ele_VetoIdDecisions",&_ele_VetoIdDecisions);
  intree->SetBranchAddress("ele_LooseIdDecisions",&_ele_LooseIdDecisions);
  intree->SetBranchAddress("ele_MediumIdDecisions",&_ele_MediumIdDecisions);
  intree->SetBranchAddress("ele_TightIdDecisions",&_ele_TightIdDecisions);


  Long64_t nentries = intree->GetEntries();

  for(int i=0;i<nentries;i++){
    
    intree->GetEntry(i);
    
    for(unsigned int i_eg=0;i_eg<10;i_eg++){      
     
      if( _ele_L1Stage2_emul_hwPt[i_eg]>0 && _ele_TightIdDecisions[i_eg]){
	
	int index=abs(_ele_L1Stage2_emul_ieta[i_eg])-1;
	if(index>26)
	  index=26;
	histos_shape[index]->Fill(_ele_L1Stage2_emul_shape[i_eg]);

      }
      
    }      
	
  }

  for(unsigned int i=0; i<histos_shape.size();i++)
    histos_shape[i]->Write();


  f_new->Close();


  return;
}









void Run260627(){

  vector<TString> input_names;
  input_names.push_back("/data_CMS/cms/strebler/L1_eg/Run260627/prod_151125/tree_Stage2_Run260627_1.root");
  input_names.push_back("/data_CMS/cms/strebler/L1_eg/Run260627/prod_151125/tree_Stage2_Run260627_2.root");
  input_names.push_back("/data_CMS/cms/strebler/L1_eg/Run260627/prod_151125/tree_Stage2_Run260627_3.root");
  TString output_name="/data_CMS/cms/strebler/L1_eg/Run260627/prod_151125/tree_Stage2_forcalibr.root";

  calibrtree_maker(input_names,output_name);

}





void Run260627_shapeID(){

  vector<TString> input_names;
  input_names.push_back("/data_CMS/cms/strebler/L1_eg/Run260627/prod_151127/tree_Stage2_Run260627_1.root");
  input_names.push_back("/data_CMS/cms/strebler/L1_eg/Run260627/prod_151127/tree_Stage2_Run260627_2.root");
  input_names.push_back("/data_CMS/cms/strebler/L1_eg/Run260627/prod_151127/tree_Stage2_Run260627_3.root");
  TString output_name="/data_CMS/cms/strebler/L1_eg/Run260627/prod_151127/tree_Stage2_forshapeID.root";

  calibrtree_maker(input_names,output_name);

}





void Run260627_HoverE(){

  vector<TString> input_names;
  input_names.push_back("/data_CMS/cms/strebler/L1_eg/Run260627/prod_151222/tree_Stage2_Run260627.root");
  TString output_name="/data_CMS/cms/strebler/L1_eg/Run260627/prod_151222/tree_Stage2_forHoverE.root";

  calibrtree_maker(input_names,output_name);

}






void Run260627_Iso(){

  vector<TString> input_names;
  input_names.push_back("/data_CMS/cms/strebler/L1_eg/Run260627/prod_160112/tree_Stage2.root");
  TString output_name="/data_CMS/cms/strebler/L1_eg/Run260627/prod_160112/eletree_Stage2_forIso.root";

  calibrtree_maker(input_names,output_name);

}







void Run260627_Shape_trimming(){

  vector<TString> input_names;
  input_names.push_back("/data_CMS/cms/strebler/L1_eg/Run260627/prod_160112/tree_Stage2.root");
  TString output_name="/data_CMS/cms/strebler/L1_eg/Run260627/prod_160112/histos_shape_forTrimming.root";

  histoshape_maker(input_names,output_name);

}






void Run260627_Iso_newTrim10(){

  vector<TString> input_names;
  input_names.push_back("/data_CMS/cms/strebler/L1_eg/Run260627/DoubleEG_newTrim10_19.01.2016_260627/tree_Stage2.root");
  TString output_name="/data_CMS/cms/strebler/L1_eg/Run260627/DoubleEG_newTrim10_19.01.2016_260627/eletree_Stage2_forIso.root";

  calibrtree_maker(input_names,output_name);

}






void Run260627_newLayer1(){

  vector<TString> input_names;
  input_names.push_back("/data_CMS/cms/strebler/L1_eg/Run260627/DoubleEG_l1t_integr_v7_Layer1_caloff_26.02.2016_260627/L1Ntuple.root");
  TString output_name="/data_CMS/cms/strebler/L1_eg/Run260627/DoubleEG_l1t_integr_v7_Layer1_caloff_26.02.2016_260627/eletree_Stage2_newLayer1.root";

  calibrtree_maker(input_names,output_name);

}







void Run260627_oldLayer1(){

  vector<TString> input_names;
  input_names.push_back("/data_CMS/cms/strebler/L1_eg/Run260627/DoubleEG_l1t_integr_v7_dummyLayer1_caloff_29.02.2016_260627/L1Ntuple.root");
  TString output_name="/data_CMS/cms/strebler/L1_eg/Run260627/DoubleEG_l1t_integr_v7_dummyLayer1_caloff_29.02.2016_260627/eletree_Stage2_newLayer1.root";

  calibrtree_maker(input_names,output_name);

}




void Run260627_Layer1_160310(){

  vector<TString> input_names;
  input_names.push_back("/data_CMS/cms/strebler/L1_eg/Run260627/DoubleEG_dasu-updates-for-layer1_10.03.2016_260627/L1Ntuple.root");
  TString output_name="/data_CMS/cms/strebler/L1_eg/Run260627/DoubleEG_dasu-updates-for-layer1_10.03.2016_260627/eletree_Stage2_Layer1_160310.root";

  calibrtree_maker(input_names,output_name);

}





void Run260627_Layer1_160313(){

  vector<TString> input_names;
  input_names.push_back("/data_CMS/cms/strebler/L1_eg/Run260627/DoubleEG_l1t_integration_v13_Layer1_caloff_13.03.2016_260627/L1Ntuple.root");
  TString output_name="/data_CMS/cms/strebler/L1_eg/Run260627/DoubleEG_l1t_integration_v13_Layer1_caloff_13.03.2016_260627/eletree_Stage2_Layer1_160313.root";

  calibrtree_maker(input_names,output_name);

}




void Run260627_Layer1_160313_etastrip_excl(){

  vector<TString> input_names;
  input_names.push_back("/data_CMS/cms/strebler/L1_eg/Run260627/DoubleEG_l1t_integration_v13_Layer1_caloff_13.03.2016_260627/L1Ntuple.root");
  TString output_name="/data_CMS/cms/strebler/L1_eg/Run260627/DoubleEG_l1t_integration_v13_Layer1_caloff_13.03.2016_260627/eletree_Stage2_Layer1_160313_etastrip_excl.root";

  calibrtree_maker_etastrip_excluded(input_names,output_name);

}





void Run260627_Layer1_160314(){

  vector<TString> input_names;
  input_names.push_back("/data_CMS/cms/strebler/L1_eg/Run260627/DoubleEG_dasu-asymmetry-fix_caloff_14.03.2016_260627/L1Ntuple.root");
  TString output_name="/data_CMS/cms/strebler/L1_eg/Run260627/DoubleEG_dasu-asymmetry-fix_caloff_14.03.2016_260627/eletree_Stage2_Layer1_160314.root";

  calibrtree_maker(input_names,output_name);

}





void Run260627_Layer1_caloff_160314(){

  vector<TString> input_names;
  input_names.push_back("/data_CMS/cms/strebler/L1_eg/Run260627/DoubleEG_dasu-asymmetry-fix_Layer1_caloff_14.03.2016_260627/L1Ntuple.root");
  TString output_name="/data_CMS/cms/strebler/L1_eg/Run260627/DoubleEG_dasu-asymmetry-fix_Layer1_caloff_14.03.2016_260627/eletree_Stage2_Layer1_caloff_160314.root";

  calibrtree_maker(input_names,output_name);

}






void DY_MC_020414(){

  vector<TString> input_names;
  input_names.push_back("/data_CMS/cms/strebler/L1_eg/DY_MC/DY_l1t_tsg_v3_firm_match_02.04.2016/L1Ntuple.root");
  TString output_name="/data_CMS/cms/strebler/L1_eg/DY_MC/DY_l1t_tsg_v3_firm_match_02.04.2016/eletree_Stage2_MC_020414.root";

  calibrtree_maker(input_names,output_name);

}





void DY_MC_layer1_050414(){

  vector<TString> input_names;
  input_names.push_back("/data_CMS/cms/strebler/L1_eg/DY_MC/DY_l1t_integr_v20_layer1_05.04.2016/L1Ntuple.root");
  TString output_name="/data_CMS/cms/strebler/L1_eg/DY_MC/DY_l1t_integr_v20_layer1_05.04.2016/eletree_Stage2_MC_layer1_050414.root";

  calibrtree_maker(input_names,output_name);

}




void Run260627_Layer1_trim15_caloff_160411(){

  vector<TString> input_names;
  input_names.push_back("/data_CMS/cms/strebler/L1_eg/Run260627/DoubleEG_l1t_integr_v20_layer1_trim15_nocal_11.04.2016_260627/L1Ntuple.root");
  TString output_name="/data_CMS/cms/strebler/L1_eg/Run260627/DoubleEG_l1t_integr_v20_layer1_trim15_nocal_11.04.2016_260627/eletree_Stage2_Run260627_layer1_trim15_110416.root";

  calibrtree_maker(input_names,output_name);

}





void Run260627_Layer1_trim20_caloff_160411(){

  vector<TString> input_names;
  input_names.push_back("/data_CMS/cms/strebler/L1_eg/Run260627/DoubleEG_l1t_integr_v20_layer1_trim20_nocal_12.04.2016_260627/L1Ntuple.root");
  TString output_name="/data_CMS/cms/strebler/L1_eg/Run260627/DoubleEG_l1t_integr_v20_layer1_trim20_nocal_12.04.2016_260627/eletree_Stage2_Run260627_layer1_trim20_110416.root";

  calibrtree_maker(input_names,output_name);

}




void Run260627_Layer1_trim25_caloff_160413(){

  vector<TString> input_names;
  input_names.push_back("/data_CMS/cms/strebler/L1_eg/Run260627/DoubleEG_l1t_integr_v20_layer1_trim25_shapeID99_nocal_13.04.2016_260627/L1Ntuple.root");
  TString output_name="/data_CMS/cms/strebler/L1_eg/Run260627/DoubleEG_l1t_integr_v20_layer1_trim25_shapeID99_nocal_13.04.2016_260627/eletree_Stage2_Run260627_layer1_trim25_130416.root";

  calibrtree_maker(input_names,output_name);

}




void Run260627_Layer1_trim20_50_caloff_160414(){

  vector<TString> input_names;
  input_names.push_back("/data_CMS/cms/strebler/L1_eg/Run260627/DoubleEG_l1t_integr_v34_trim_20_50_nocal_14.04.2016_260627/L1Ntuple.root");
  TString output_name="/data_CMS/cms/strebler/L1_eg/Run260627/DoubleEG_l1t_integr_v34_trim_20_50_nocal_14.04.2016_260627/eletree_Stage2_Run260627_layer1_trim20_50_140416.root";

  calibrtree_maker(input_names,output_name);

}






void DY_MC_layer1_150414(){

  vector<TString> input_names;
  input_names.push_back("/data_CMS/cms/strebler/L1_eg/DY_MC/DY_l1t_integr_v34_trim20_50_noshapeEE_cal_calmax_1.25_noshapeID_v16.04.15/L1Ntuple.root");
  TString output_name="/data_CMS/cms/strebler/L1_eg/DY_MC/DY_l1t_integr_v34_trim20_50_noshapeEE_cal_calmax_1.25_noshapeID_v16.04.15/eletree_Stage2_MC_layer1_150414.root";

  calibrtree_maker(input_names,output_name);

}






void DY_MC_layer1_160414(){

  vector<TString> input_names;
  input_names.push_back("/data_CMS/cms/strebler/L1_eg/DY_MC/DY_l1t_integr_v37.1_shapeEBrelax_16.04.2016/L1Ntuple.root");
  TString output_name="/data_CMS/cms/strebler/L1_eg/DY_MC/DY_l1t_integr_v37.1_shapeEBrelax_16.04.2016/eletree_Stage2_MC_layer1_160414.root";

  calibrtree_maker(input_names,output_name);

}







void DoubleEG_calibr_022017(){

  vector<TString> input_names;
  input_names.push_back("/data_CMS/cms/beschi/L1Ntuple_283478.root");
  TString output_name="/data_CMS/cms/strebler/L1_eg/Calibr_022017/eletree_Run283478.root";

  calibrtree_maker(input_names,output_name);

}








void SingleEle_calibr_062017_v1_9(){

  vector<TString> input_names;
  input_names.push_back("/data_CMS/cms/strebler/L1_eg/Calibr_062017/L1Ntuple_SingleEle_283548_v1_9.root");
  TString output_name="/data_CMS/cms/strebler/L1_eg/Calibr_062017/eletree_Run283548_v1_9.root";

  calibrtree_maker(input_names,output_name);

}






void SingleEle_calibr_062017_v1_9_mean(){

  vector<TString> input_names;
  input_names.push_back("/data_CMS/cms/strebler/L1_eg/Calibr_062017/L1Ntuple_SingleEle_283548_v1_9_mean.root");
  TString output_name="/data_CMS/cms/strebler/L1_eg/Calibr_062017/eletree_Run283548_v1_9_mean.root";

  calibrtree_maker(input_names,output_name);

}

*/
