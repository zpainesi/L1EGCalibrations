#include "Util.h"

using namespace std;

//Double_t getIntegral(TGraphAsymmErrors* aGraph, Double_t midX,Double_t lX,Double_t rX,TString);
//void tokenize(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters) ;

void processOptionFile(TString fileName,TString descHName,TString baselineFileName,TString baselineHistName,Double_t baselineEt=24,TString prefix="")
{
    std::cout<<"fileName : "<<fileName<<"\n";
    std::cout<<"descHName : "<<descHName<<"\n";
    std::cout<<"baselineFileName : "<<baselineFileName<<"\n";
    std::cout<<"baselineHistName : "<<baselineHistName<<"\n";
    std::cout<<"prefix : "<<prefix<<"\n";
    TGraphAsymmErrors * baselineTurnOn=nullptr;
    TFile * baselineFile=TFile::Open(baselineFileName,"READ");
    if(baselineFile)
    {
        baselineTurnOn=(TGraphAsymmErrors *)baselineFile->Get(baselineHistName);
        if( not baselineHistName)
        {
            std::cout<<"Baseline histo not found !! ["<<baselineHistName<<"] from "<<baselineFileName<<"\n";
        }
        else
        {
            std::cout<<" Obtained baseline as  : "<<baselineTurnOn->GetName()<<"\n";
        }
    }
    else
    {
        std::cout<<"Baseline File Not Found : "<<baselineFileName<<"\n";
    }

    TFile * aFile=TFile::Open(fileName,"READ");
    TGraphAsymmErrors * graphToIntegrate(nullptr);
    TString folderName;
    auto descriptionHistogram = (TH1F*) aFile->Get(descHName);
    Double_t eT_threshold,left_DX,right_DX,area;
    string cmd;
    Bool_t isBetterThanLoose(false);
    Bool_t isBetterThanTight(false);
    eT_threshold=30;
    left_DX=7.0;
    right_DX=7.0;

    // Getting the fixed Rate Metrics Printed
    auto xAxis = descriptionHistogram->GetXaxis();
    auto nBins = descriptionHistogram->GetNbinsX();
    std::vector<std::string> tockens;
    if(baselineTurnOn)
    {
        eT_threshold=baselineEt;
        area=getIntegral(baselineTurnOn, eT_threshold,eT_threshold-left_DX,eT_threshold+right_DX,folderName,nullptr);
        std::cout<<"Baseline Area = "<<area<<"\n";
    }
    std::cout<<"idx , option , eTMin, eff , eTMax, area, isBetterThanLoose, isBetterThanTight\n";
    int thresoldPassEvents(0),betterThanPassEvents(0);
    for(Int_t i=0;i<nBins;i++)
    {
        //if(i>5) break;
        TString histname(xAxis->GetBinLabel(i));
        std::string hname(xAxis->GetBinLabel(i));
        
        if(hname.size()==0) continue;
        graphToIntegrate=(TGraphAsymmErrors *) aFile->Get("turn_on_progression/"+histname);
        if( graphToIntegrate){
        true;
        }
        else 
        {
            std::cout<<"\n\tturn_on_progression/"+histname<<" not available "<<"\n";
            continue;
        }

        tockens.clear();
        tokenize(hname,tockens,"_");
        std::replace( tockens[4].begin(), tockens[4].end(), 'p', '.');
        std::replace( tockens[5].begin(), tockens[5].end(), 'p', '.');
        std::replace( tockens[6].begin(), tockens[6].end(), 'p', '.');
        eT_threshold=descriptionHistogram->GetBinContent(i);
        if(baselineTurnOn)
        {
            isBetterThanTight=isGoodTurnON(baselineTurnOn,graphToIntegrate,29,54);
        }   
 
       // std::cout<<histname<<"  auto x = (TGraphAsymmErrors *) aFile->Get(\"turn_on_progression/"<<histname<<"\")";
       // std::cout<<"\n\t"<<i<<" isBetterThanTight : "<<isBetterThanTight<<"\n";
        if( eT_threshold > baselineEt+2.0 ) continue;
        thresoldPassEvents++;
        if( not isBetterThanTight) continue;
        betterThanPassEvents++;
        
        folderName=prefix+"/"+tockens[3]+"/";
        cmd="mkdir -p "+folderName;
        system(cmd.c_str());
        
        area=getIntegral(graphToIntegrate, eT_threshold,eT_threshold-left_DX,eT_threshold+right_DX,folderName,baselineTurnOn);
        std::cout<<i<<" , "<<tockens[3]
                 <<" ,  "
                 <<tockens[4]<<" , "
                 <<tockens[5]<<" , "
                 <<tockens[6]<<"  "
                 <<" , "<<eT_threshold 
                 <<" , "<<area
                 <<" , "<<isBetterThanLoose
                 <<" , "<<isBetterThanTight
                 <<"\n";
    }
    
    aFile->Close();
    std::cout<<" Number of options passing the \" eT thresold \" : "<<thresoldPassEvents<<" / "<<nBins<<"\n";
    std::cout<<" Number of options passing the \" better than the baseline \" : "<<betterThanPassEvents<<" / "<<nBins<<"\n";
}

int  main(int argc,char *argv[])
{

    if(argc<3)
    {
        std::cout<<" Usage : \n"
                 <<"       ./optEval.exe <fname> <desc_hist>  prefix_to_save\n"
                 <<"eg. :  ./integrate.exe  HistgramFile_step3step4_eg_12X_recent_gs.root turn_on_progression/divide_pt_pass_option_1_Et_25_fr_by_pT_all 25 5 5 plots/ \n"
                 <<"\n";
                 
        exit(1);
    }

    processOptionFile(argv[1],argv[2],argv[3],argv[4],24.0,argv[5]);
    return 0;

}


