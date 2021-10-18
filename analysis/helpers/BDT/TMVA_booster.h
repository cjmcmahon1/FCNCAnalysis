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

using namespace std;

map<string, float> get_params(string path_to_csv) {
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
    map<string, float> hyperparam_map;
    while (!fin.eof()){
        fin >> line;
        if (line.length() > 0) {
            stringstream s(line);
            for(string i : parameters){
                std::getline(s, record, delimiter);//get next entry in csv
                hyperparam_map[i] = stof(record);
            }
        }
        line.clear();
    }
    return hyperparam_map;
}

map<string, float> test = get_params("HCT_params.csv");

void learn(int nTrain)
{
    // Initialize TMVA
    TMVA::Tools::Instance();
    TFile* outputFile = TFile::Open("BDT.root", "RECREATE");
    TMVA::Factory *factory = new TMVA::Factory("TMVA", outputFile, "V:DrawProgressBar=True:Transformations=I;D;P;G:AnalysisType=Classification");         
    TString path = "/hadoop/cms/store/user/smay/DeepIsolation/TTbar_DeepIso_v0.0.8/merged_ntuple_*.root"; // rounded kinematics
    TChain* chain = new TChain("t");
    chain->Add(path);
    TObjArray *listOfFiles = chain->GetListOfFiles();
    TIter fileIter(listOfFiles);
    TFile *currentFile = 0;

    vector<TString> vFiles;

    while ( (currentFile = (TFile*)fileIter.Next()) ) 
        vFiles.push_back(currentFile->GetTitle());

    vector<TFile*> vFileSig;
    vector<TFile*> vFileBkg;
    vector<TTree*> vTreeSig;
    vector<TTree*> vTreeBkg;
                        
    for (int i = 0; i < vFiles.size(); i++) {
        vFileSig.push_back(TFile::Open(vFiles[i]));
        vFileBkg.push_back(TFile::Open(vFiles[i]));
        vTreeSig.push_back((TTree*)vFileSig[i]->Get("t"));
        vTreeBkg.push_back((TTree*)vFileBkg[i]->Get("t"));
        factory->AddSignalTree(vTreeSig[i]);
        factory->AddBackgroundTree(vTreeBkg[i]);
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
        factory->AddVariable(feat, 'F');
    }
    Double_t signalWeight     = 1.0;
    Double_t backgroundWeight = 1.0;
    factory->AddVariable("lepton_eta", 'F');
    factory->AddVariable("lepton_phi", 'F');
    factory->AddVariable("lepton_pt", 'F');

    float nTrainF = nTrain;
    float nTrainSigF = nTrainF*0.5;
    float nTrainBkgF = nTrainF*0.5;

    int nTrainSig = (int) nTrainSigF;
    int nTrainBkg = (int) nTrainBkgF;

    int nTestSig = 100000;
    int nTestBkg = 100000;

    TString prepare_events = "nTrain_Signal=" + to_string(nTrainSig) + ":nTrain_Background=" + to_string(nTrainBkg) + ":nTest_Signal=" + to_string(nTestSig) + ":nTest_Background=" + to_string(nTestBkg) + ":SplitMode=Random:NormMode=NumEvents:!V";   

    //factory->PrepareTrainingAndTestTree("lepton_isFromW==1&&lepton_flavor==0", "lepton_isFromW==0&&lepton_flavor==0", prepare_events); // electrons

    factory->SetSignalWeightExpression("1");
    factory->SetBackgroundWeightExpression("1");
    
    TString option = "!H:V:NTrees=500:BoostType=Grad:Shrinkage=0.10:!UseBaggedGrad:nCuts=2000:MinNodeSize=0.1%:PruneStrength=5:PruneMethod=CostComplexity:MaxDepth=5:CreateMVAPdfs";
    
    factory->BookMethod(TMVA::Types::kBDT, "BDT", option);
    factory->TrainAllMethods();
    factory->TestAllMethods();
    factory->EvaluateAllMethods();
    
    outputFile->Close();
    std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
    std::cout << "==> TMVAClassification is done!" << std::endl;
    delete factory;
}



