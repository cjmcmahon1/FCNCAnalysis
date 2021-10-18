#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <fstream>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/DataLoader.h"

using namespace std;

map<string, string> get_params(string path_to_csv) {
    ifstream fin;
    fin.open(path_to_csv);
    std::string line;
    char delimiter = ',';
    std::vector<string> parameters;
    string record;
    fin >> line;
    stringstream head(line);
    while (getline(head, record, delimiter)){
        parameters.push_back(record);
    }
    map<string, string> hyperparam_map;
    while (!fin.eof()){
        fin >> line;
        if (line.length() > 0) {
            stringstream s(line);
            for(string i : parameters){
                std::getline(s, record, delimiter);//get next entry in csv
                hyperparam_map[i] = record;
            }
        }
        line.clear();
    }
    return hyperparam_map;
}

//map<string, float> test = get_params("HCT_params.csv");
string get_BDT_options(map<string, string> hyperparam_map, bool debug=false){ 
    //TString option = "!H:V:NTrees=500:BoostType=Grad:Shrinkage=0.10:!UseBaggedGrad:nCuts=2000:MinNodeSize=0.1%:PruneStrength=5:PruneMethod=CostComplexity:MaxDepth=5:CreateMVAPdfs";
    string option = "!H:V:NTrees=" + hyperparam_map["n_estimators"]; //number of trees
    option += ":BoostType=Grad:Shrinkage=" + hyperparam_map["learning_rate"]; //learning rate
    option += "!UseBaggedGrad";
    option += ":nCuts=200"; //can't find an XGB eqivalent
    option += ":MinNodeSize=5%";
    option += ":MaxDepth=" + hyperparam_map["max_depth"];
    option += ":CreateMVAPdfs";
    if (debug) cout << "option = " << option << endl;
    return option;
}


void learn(int nTrain)
{
    // Initialize TMVA
    TMVA::Tools::Instance();
    TFile* outputFile = TFile::Open("BDT.root", "RECREATE");
    TMVA::Factory *factory = new TMVA::Factory("TMVA", outputFile, "V:DrawProgressBar=True:Transformations=I;D;P;G:AnalysisType=Classification");         
    TMVA::DataLoader *dataloader = new TMVA::DataLoader("datasetBkg");
    //TString path = "/hadoop/cms/store/user/smay/DeepIsolation/TTbar_DeepIso_v0.0.8/merged_ntuple_*.root"; // rounded kinematics
    TChain* chain_sig = new TChain("T");
    TChain* chain_bkg = new TChain("T");
    string base_dir = "/home/users/cmcmahon/fcnc/ana/analysis/helpers/BDT/babies/data_driven_v2/";
    for (int y : {2016, 2017, 2018}) {
        for (string s : {"tch", "tuh"}) {
            chain_sig->Add(Form("%s%d/MC/signal_%s.root", base_dir.c_str(), y, s.c_str()));
        }
        chain_bkg->Add(Form("%s%d/MC/fakes_mc.root", base_dir.c_str(), y));
        chain_bkg->Add(Form("%s%d/MC/flips_mc.root", base_dir.c_str(), y));
        chain_bkg->Add(Form("%s%d/MC/rares.root", base_dir.c_str(), y));
    }
    //TString signal_path = "/home/users/cmcmahon/fcnc/ana/analysis/helpers/BDT/babies/data_driven_v2/2016/MC/signal_tch.root"; // rounded kinematics
    //TString fakes_path = "/home/users/cmcmahon/fcnc/ana/analysis/helpers/BDT/babies/data_driven_v2/20*/MC/fakes_mc.root"; // rounded kinematics
    //TString flips_path = "/home/users/cmcmahon/fcnc/ana/analysis/helpers/BDT/babies/data_driven_v2/20*/MC/flips_mc.root"; // rounded kinematics
    //TString rares_path = "/home/users/cmcmahon/fcnc/ana/analysis/helpers/BDT/babies/data_driven_v2/20*/MC/rares.root"; // rounded kinematics
    //chain_sig->Add(signal_path);
    //chain_bkg->Add(fakes_path);
    //chain_bkg->Add(flips_path);
    //chain_bkg->Add(rares_path);
    TObjArray *listOfFiles_sig = chain_sig->GetListOfFiles();
    TObjArray *listOfFiles_bkg = chain_bkg->GetListOfFiles();
    cout << "num signal entries     = " << chain_sig->GetEntries() << endl;
    cout << "num background entries = " << chain_bkg->GetEntries() << endl;
    TIter fileIter_sig(listOfFiles_sig);
    TIter fileIter_bkg(listOfFiles_bkg);
    TFile *currentFile_sig = 0;
    TFile *currentFile_bkg = 0;
    
    vector<TString> vFiles_sig;
    //cout << currentFile_sig->GetTitle() << endl;
    cout << (TFile*)fileIter_sig.Next() << endl;

    while ( (currentFile_sig = (TFile*)fileIter_sig.Next()) ) {
        cout << currentFile_sig->GetTitle() << endl;
        vFiles_sig.push_back(currentFile_sig->GetTitle());
    }

    vector<TString> vFiles_bkg;
    //cout << currentFile_bkg->GetTitle();
    cout << (TFile*)fileIter_sig.Next() << endl;

    while ( (currentFile_bkg = (TFile*)fileIter_bkg.Next()) ) {
        cout << currentFile_bkg->GetTitle() << endl;
        vFiles_bkg.push_back(currentFile_bkg->GetTitle());
    }

    vector<TFile*> vFileSig;
    vector<TFile*> vFileBkg;
    vector<TTree*> vTreeSig;
    vector<TTree*> vTreeBkg;
                        
    for (unsigned int i = 0; i < vFiles_sig.size(); i++) {
        vFileSig.push_back(TFile::Open(vFiles_sig[i]));
        vTreeSig.push_back((TTree*)vFileSig[i]->Get("T"));
        dataloader->AddSignalTree(vTreeSig[i]);
    }
    for (unsigned int i = 0; i < vFiles_bkg.size(); i++) {
        vFileBkg.push_back(TFile::Open(vFiles_bkg[i]));
        vTreeBkg.push_back((TTree*)vFileBkg[i]->Get("T"));
        dataloader->AddBackgroundTree(vTreeBkg[i]);
    }

    std::vector<std::string> BDT_features = {
        "Most_Forward_pt",
        "HT",
        "LeadLep_eta",
        "LeadLep_pt",
        "LeadLep_dxy",
        "LeadLep_dz",
        "SubLeadLep_pt",
        "SubLeadLep_eta",
        "SubLeadLep_dxy",
        "SubLeadLep_dz",
        "nJets",
        "nBtag",
        "LeadJet_pt",
        "SubLeadJet_pt",
        "SubSubLeadJet_pt",
        "LeadJet_BtagScore",
        "SubLeadJet_BtagScore",
        "SubSubLeadJet_BtagScore",
        "nElectron",
        "MET_pt",
        "LeadBtag_pt",
        "MT_LeadLep_MET",
        "MT_SubLeadLep_MET",
        "LeadLep_SubLeadLep_Mass",
        "SubSubLeadLep_pt",
        "SubSubLeadLep_eta",
        "SubSubLeadLep_dxy",
        "SubSubLeadLep_dz",
        "MT_SubSubLeadLep_MET",
        "LeadBtag_score"
    };
    for (std::string feat : BDT_features) {
        dataloader->AddVariable(feat, 'F');
    }
    Double_t signalWeight     = 1.0;
    Double_t backgroundWeight = 1.0;
    //dataloader->AddVariable("lepton_eta", 'F');
    //dataloader->AddVariable("lepton_phi", 'F');
    //dataloader->AddVariable("lepton_pt", 'F');

    float nTrainF = nTrain;
    float nTrainSigF = nTrainF*0.5;
    float nTrainBkgF = nTrainF*0.5;

    int nTrainSig = (int) nTrainSigF;
    int nTrainBkg = (int) nTrainBkgF;

    int nTestSig = 100000;
    int nTestBkg = 100000;

    TString prepare_events = "nTrain_Signal=" + to_string(nTrainSig) + ":nTrain_Background=" + to_string(nTrainBkg) + ":nTest_Signal=" + to_string(nTestSig) + ":nTest_Background=" + to_string(nTestBkg) + ":SplitMode=Random:NormMode=NumEvents:!V";   

    //factory->PrepareTrainingAndTestTree("lepton_isFromW==1&&lepton_flavor==0", "lepton_isFromW==0&&lepton_flavor==0", prepare_events); // electrons

    dataloader->SetSignalWeightExpression("weight/2.5");
    dataloader->SetBackgroundWeightExpression("weight");
    
    TString option = "!H:V:NTrees=2:BoostType=Grad:Shrinkage=0.10:!UseBaggedGrad:nCuts=2:MinNodeSize=0.1%:MaxDepth=1";
    
    factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDT", option);
    factory->TrainAllMethods();
    factory->TestAllMethods();
    factory->EvaluateAllMethods();
    
    outputFile->Close();
    std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
    std::cout << "==> TMVAClassification is done!" << std::endl;
    delete factory;
}



