#include <iostream>
#include "TMVA/Reader.h"
#include "../../../../NanoTools/NanoCORE/SSSelections.h"
#include "../../../../NanoTools/NanoCORE/Nano.h"

using namespace std;
using namespace std::chrono;

double convert_tmva_to_prob(double score) {
    // Undo TMVA transformation
    double raw_score = -0.5 * log( (2 / (score + 1)) - 1);
    // Apply logistic (sigmoid) transformation
    double prob = 1 / (1 + exp(-raw_score));
    return prob;
}

class BDT {
    // Boosted decision tree class. Allows evaluation of BDT using ROOT variables
    // by reading from a customized xml file. Uses the TMVA library:
    // https://root.cern.ch/download/doc/tmva/TMVAUsersGuide.pdf
    unique_ptr<TMVA::Reader> booster;
    std::map<std::string, Float_t> parameter_map;
    public:
        BDT(std::string);
        void set_features(std::map<std::string, Float_t>, bool);
        Float_t get_score();
        std::map<std::string, Float_t> calculate_features(Jets, Jets, float, Leptons);
};

BDT::BDT(std::string path_to_xml) {
    booster.reset( new TMVA::Reader( "!Color:Silent" ) );
    // Booster must be initialized with an xml file and the feature addresses.
    // The feature addresses cannot be changed, but the values can.
    // NOTE: you must add the booster features in the same order
    // that they are organized in the xml file
    booster->AddVariable("Most_Forward_pt", &(parameter_map["Most_Forward_pt"]));
    booster->AddVariable("HT", &(parameter_map["HT"]));
    booster->AddVariable("LeadLep_eta", &(parameter_map["LeadLep_eta"]));
    booster->AddVariable("LeadLep_pt", &(parameter_map["LeadLep_pt"]));
    booster->AddVariable("LeadLep_dxy", &(parameter_map["LeadLep_dxy"]));
    booster->AddVariable("LeadLep_dz", &(parameter_map["LeadLep_dz"]));
    booster->AddVariable("SubLeadLep_pt", &(parameter_map["SubLeadLep_pt"]));
    booster->AddVariable("SubLeadLep_eta", &(parameter_map["SubLeadLep_eta"]));
    booster->AddVariable("SubLeadLep_dxy", &(parameter_map["SubLeadLep_dxy"]));
    booster->AddVariable("SubLeadLep_dz", &(parameter_map["SubLeadLep_dz"]));
    booster->AddVariable("nJet", &(parameter_map["nJets"]));
    booster->AddVariable("nbtag", &(parameter_map["nBtag"]));
    booster->AddVariable("LeadJet_pt", &(parameter_map["LeadJet_pt"]));
    booster->AddVariable("SubLeadJet_pt", &(parameter_map["SubLeadJet_pt"]));
    booster->AddVariable("SubSubLeadJet_pt", &(parameter_map["SubSubLeadJet_pt"]));
    booster->AddVariable("nElectron", &(parameter_map["nElectron"]));
    booster->AddVariable("MET_pt", &(parameter_map["MET_pt"]));
    booster->AddVariable("LeadBtag_pt", &(parameter_map["LeadBtag_pt"]));
    booster->AddVariable("MT_LeadLep_MET", &(parameter_map["MT_LeadLep_MET"]));
    booster->AddVariable("MT_SubLeadLep_MET", &(parameter_map["MT_SubLeadLep_MET"]));
    booster->AddVariable("LeadLep_SubLeadLep_Mass", &(parameter_map["LeadLep_SubLeadLep_Mass"]));
    //booster->AddVariable("SubSubLeadLep_pt", &(parameter_map["SubSubLeadLep_pt"]));
    //booster->AddVariable("SubSubLeadLep_eta", &(parameter_map["SubSubLeadLep_eta"]));
    //booster->AddVariable("SubSubLeadLep_dxy", &(parameter_map["SubSubLeadLep_dxy"]));
    //booster->AddVariable("SubSubLeadLep_dz", &(parameter_map["SubSubLeadLep_dz"]));
    //booster->AddVariable("MT_SubSubLeadLep_MET", &(parameter_map["MT_SubSubLeadLep_MET"]));
    //booster->AddVariable("LeadBtag_score", &(parameter_map["LeadBtag_score"]));
    booster->BookMVA("BDT", path_to_xml);

}

std::map<std::string, Float_t> BDT::calculate_features(Jets good_jets, Jets good_bjets, float ht, Leptons ordered_leptons) {
    Float_t MET_pt = nt.MET_pt();
    Lepton LeadLep = ordered_leptons[0];
    Lepton SubLeadLep = ordered_leptons[1];
    Float_t LeadLep_pt = LeadLep.pt();
    Float_t SubLeadLep_pt = SubLeadLep.pt();
    Float_t LeadLep_eta = abs(LeadLep.eta());
    Float_t SubLeadLep_eta = abs(SubLeadLep.eta());
    Float_t LeadLep_dxy = abs(LeadLep.dxy());
    Float_t SubLeadLep_dxy = abs(SubLeadLep.dxy());
    Float_t LeadLep_dz = abs(LeadLep.dz());
    Float_t SubLeadLep_dz = abs(SubLeadLep.dz());
    Float_t LeadLep_SubLeadLep_Mass = (LeadLep.p4() + SubLeadLep.p4()).M();
    Float_t njets = good_jets.size();
    Float_t nbjets = good_bjets.size();
    Float_t LeadJet_pt=0., SubLeadJet_pt=0., SubSubLeadJet_pt=0., LeadBtag_pt=0., LeadBtag_score=0.;
    Float_t nElectron = 0;//ordered_leptons.size();
    //thiird lepton properties
    Float_t SubSubLeadLep_pt = 0., SubSubLeadLep_eta = 0., SubSubLeadLep_dxy = 0., SubSubLeadLep_dz = 0.;
    Float_t MT_SubSubLeadLep_MET = 0.;
    if (ordered_leptons.size() > 2) {
        Lepton SubSubLeadLep = ordered_leptons[2];
        SubSubLeadLep_pt = SubSubLeadLep.pt();
        SubSubLeadLep_eta = abs(SubSubLeadLep.eta());
        SubSubLeadLep_dxy = abs(SubSubLeadLep.dxy());
        SubSubLeadLep_dz = abs(SubSubLeadLep.dz());
        MT_SubSubLeadLep_MET = TMath::Sqrt(2*ordered_leptons[2].pt()*nt.MET_pt() * (1 - TMath::Cos(ordered_leptons[2].phi()-nt.MET_phi())));
        }

    for(auto lep: ordered_leptons){//count nElectrons
        if (abs(lep.id())==11) {
            nElectron++;
            }
        }
    if (njets > 2) {
        LeadJet_pt = good_jets[0].pt();
        SubLeadJet_pt = good_jets[1].pt();
        SubSubLeadJet_pt = good_jets[2].pt();
    }
    else if (njets == 2) {
        LeadJet_pt = good_jets[0].pt();
        SubLeadJet_pt = good_jets[1].pt();
    }
    if (nbjets > 0) { 
        LeadBtag_pt = good_bjets[0].pt();
        LeadBtag_score = good_bjets[0].bdisc();
    }
    Float_t Most_Forward_pt = 0., highest_abs_eta=0.;
    for(int i=0; i < njets; i++){
        if (abs(good_jets[i].eta()) >= highest_abs_eta) {
            highest_abs_eta = abs(good_jets[i].eta());
            Most_Forward_pt = good_jets[i].pt();
        }
    }
    Float_t BDT_HT = ht;
    Float_t MT_LeadLep_MET, MT_SubLeadLep_MET;
    MT_LeadLep_MET = TMath::Sqrt(2*ordered_leptons[0].pt()*nt.MET_pt() * (1 - TMath::Cos(ordered_leptons[0].phi()-nt.MET_phi())));
    MT_SubLeadLep_MET = TMath::Sqrt(2*ordered_leptons[1].pt()*nt.MET_pt() * (1 - TMath::Cos(ordered_leptons[1].phi()-nt.MET_phi())));
    std::map<std::string, Float_t> params = {
        {"nJets", njets},
        {"nBtag", nbjets},
        {"LeadJet_pt", LeadJet_pt},
        {"SubLeadJet_pt", SubLeadJet_pt},
        {"SubSubLeadJet_pt", SubSubLeadJet_pt},
        {"LeadBtag_pt", LeadBtag_pt},
        {"MT_LeadLep_MET", MT_LeadLep_MET},
        {"MT_SubLeadLep_MET", MT_SubLeadLep_MET},
        {"HT", BDT_HT},
        {"Most_Forward_pt", Most_Forward_pt},
        {"MET_pt", MET_pt},
        {"LeadLep_pt", LeadLep_pt},
        {"SubLeadLep_pt", SubLeadLep_pt},
        {"LeadLep_eta", LeadLep_eta},
        {"SubLeadLep_eta", SubLeadLep_eta},
        {"LeadLep_dxy", LeadLep_dxy},
        {"SubLeadLep_dxy", SubLeadLep_dxy},
        {"LeadLep_dz", LeadLep_dz},
        {"SubLeadLep_dz", SubLeadLep_dz},
        {"nElectron", nElectron},
        {"LeadLep_SubLeadLep_Mass", LeadLep_SubLeadLep_Mass},
        {"SubSubLeadLep_pt", SubSubLeadLep_pt},
        {"SubSubLeadLep_eta", SubSubLeadLep_eta},
        {"SubSubLeadLep_dxy", SubSubLeadLep_dxy},
        {"SubSubLeadLep_dz", SubSubLeadLep_dz},
        {"MT_SubSubLeadLep_MET", MT_SubSubLeadLep_MET},
        {"LeadBtag_score", LeadBtag_score}
    };
    //cout << "Setting nBtag = " << nbjets << endl;
    return params;
}

     
void BDT::set_features(std::map<std::string, Float_t> BDT_params, bool debug=false){
    parameter_map["Most_Forward_pt"] = BDT_params["Most_Forward_pt"];
    parameter_map["HT"] = BDT_params["HT"];
    parameter_map["LeadLep_eta"] = BDT_params["LeadLep_eta"];
    parameter_map["LeadLep_pt"] = BDT_params["LeadLep_pt"];
    parameter_map["LeadLep_dxy"] = BDT_params["LeadLep_dxy"];
    parameter_map["LeadLep_dz"] = BDT_params["LeadLep_dz"];
    parameter_map["SubLeadLep_pt"] = BDT_params["SubLeadLep_pt"];
    parameter_map["SubLeadLep_eta"] = BDT_params["SubLeadLep_eta"];
    parameter_map["SubLeadLep_dxy"] = BDT_params["SubLeadLep_dxy"];
    parameter_map["SubLeadLep_dz"] = BDT_params["SubLeadLep_dz"];
    parameter_map["nJets"] = BDT_params["nJets"];
    parameter_map["nBtag"] = BDT_params["nBtag"];
    parameter_map["LeadJet_pt"] = BDT_params["LeadJet_pt"];
    parameter_map["SubLeadJet_pt"] = BDT_params["SubLeadJet_pt"];
    parameter_map["SubSubLeadJet_pt"] = BDT_params["SubSubLeadJet_pt"];
    parameter_map["nElectron"] = BDT_params["nElectron"];
    parameter_map["MET_pt"] = BDT_params["MET_pt"];
    parameter_map["LeadBtag_pt"] = BDT_params["LeadBtag_pt"];
    parameter_map["MT_LeadLep_MET"] = BDT_params["MT_LeadLep_MET"];
    parameter_map["MT_SubLeadLep_MET"] = BDT_params["MT_SubLeadLep_MET"];
    parameter_map["LeadLep_SubLeadLep_Mass"] = BDT_params["LeadLep_SubLeadLep_Mass"];
    parameter_map["SubSubLeadLep_pt"] = BDT_params["SubSubLeadLep_pt"];
    parameter_map["SubSubLeadLep_eta"] = BDT_params["SubSubLeadLep_eta"];
    parameter_map["SubSubLeadLep_dxy"] = BDT_params["SubSubLeadLep_dxy"];
    parameter_map["SubSubLeadLep_dz"] = BDT_params["SubSubLeadLep_dz"];
    parameter_map["MT_SubSubLeadLep_MET"] = BDT_params["MT_SubSubLeadLep_MET"];
    parameter_map["LeadBtag_score"] = BDT_params["LeadBtag_score"];
    if (debug) {
        cout << "event: " << nt.event() << endl;
        cout << "Most_Forward_pt: " << BDT_params["Most_Forward_pt"] << endl;
        cout << "HT: " << BDT_params["HT"] << endl;
        cout << "LeadLep_eta: " << BDT_params["LeadLep_eta"] << endl;
        cout << "LeadLep_pt: " << BDT_params["LeadLep_pt"] << endl;
        cout << "LeadLep_dxy: " << BDT_params["LeadLep_dxy"] << endl;
        cout << "LeadLep_dz: " << BDT_params["LeadLep_dz"] << endl;
        cout << "SubLeadLep_pt: " << BDT_params["SubLeadLep_pt"] << endl;
        cout << "SubLeadLep_eta: " << BDT_params["SubLeadLep_eta"] << endl;
        cout << "SubLeadLep_dxy: " << BDT_params["SubLeadLep_dxy"] << endl;
        cout << "SubLeadLep_dz: " << BDT_params["SubLeadLep_dz"] << endl;
        cout << "nJets: " << parameter_map["nJets"] << endl;
        cout << "nBtag: " << parameter_map["nBtag"] << endl;
        cout << "LeadJet_pt: " << parameter_map["LeadJet_pt"] << endl;
        cout << "SubLeadJet_pt: " << parameter_map["SubLeadJet_pt"] << endl;
        cout << "SubSubLeadJet_pt: " << parameter_map["SubSubLeadJet_pt"] << endl;
        cout << "nElectron: " << BDT_params["nElectron"] << endl;
        cout << "MET_pt: " << BDT_params["MET_pt"] << endl;
        cout << "LeadBtag_pt: " << parameter_map["LeadBtag_pt"] << endl;
        cout << "MT_LeadLep_MET: " << parameter_map["MT_LeadLep_MET"] << endl;
        cout << "MT_SubLeadLep_MET: " << parameter_map["MT_SubLeadLep_MET"] << endl;
        cout << "LeadLep_SubLeadLep_Mass: " << BDT_params["LeadLep_SubLeadLep_Mass"] << endl;
        cout << "SubSubLeadLep_pt: " << BDT_params["SubSubLeadLep_pt"] << endl;
        cout << "SubSubLeadLep_eta: " << BDT_params["SubSubLeadLep_eta"] << endl;
        cout << "SubSubLeadLep_dxy: " << BDT_params["SubSubLeadLep_dxy"] << endl;
        cout << "SubSubLeadLep_dz: " << BDT_params["SubSubLeadLep_dz"] << endl;
        cout << "MT_SubSubLeadLep_MET: " << parameter_map["MT_SubSubLeadLep_MET"] << endl;
   }
}

Float_t BDT::get_score() {
    //convert the TMVA output to a bdt score [0,1] using a sigmoid
    Float_t booster_score = convert_tmva_to_prob(booster->EvaluateMVA("BDT"));
    return booster_score;
}

//old function (~20ms per event because it makes a new BDT for every event)
Float_t get_BDT_score(Leptons ordered_leptons, std::map<std::string, Float_t> BDT_params, bool debug=false){
    auto start_time = high_resolution_clock::now();
    unique_ptr<TMVA::Reader> FCNC_booster;
    Float_t event = 1.0;
    Float_t MET_pt = nt.MET_pt();
    Lepton LeadLep = ordered_leptons[0];
    Lepton SubLeadLep = ordered_leptons[1];
    Float_t LeadLep_pt = LeadLep.pt();
    Float_t SubLeadLep_pt = SubLeadLep.pt();
    Float_t LeadLep_eta = abs(LeadLep.eta());
    Float_t SubLeadLep_eta = abs(SubLeadLep.eta());
    Float_t LeadLep_dxy = abs(LeadLep.dxy());
    Float_t SubLeadLep_dxy = abs(SubLeadLep.dxy());
    Float_t LeadLep_dz = abs(LeadLep.dz());
    Float_t SubLeadLep_dz = abs(SubLeadLep.dz());
    Float_t nElectron = nt.nElectron();
    Float_t LeadLep_SubLeadLep_Mass = (LeadLep.p4() + SubLeadLep.p4()).M();
    FCNC_booster.reset( new TMVA::Reader( "!Color:Silent" ) );
    // NOTE: you must add the booster variables in the same order
    // that they are organized in the xml file
    FCNC_booster->AddVariable("Most_Forward_pt", &(BDT_params["Most_Forward_pt"]));
    FCNC_booster->AddVariable("HT", &(BDT_params["HT"]));
    FCNC_booster->AddVariable("LeadLep_eta", &LeadLep_eta);
    FCNC_booster->AddVariable("LeadLep_pt", &LeadLep_pt);
    FCNC_booster->AddVariable("LeadLep_dxy", &LeadLep_dxy);
    FCNC_booster->AddVariable("LeadLep_dz", &LeadLep_dz);
    FCNC_booster->AddVariable("SubLeadLep_pt", &SubLeadLep_pt);
    FCNC_booster->AddVariable("SubLeadLep_eta", &SubLeadLep_eta);
    FCNC_booster->AddVariable("SubLeadLep_dxy", &SubLeadLep_dxy);
    FCNC_booster->AddVariable("SubLeadLep_dz", &SubLeadLep_dz);
    FCNC_booster->AddVariable("nJet", &(BDT_params["nJets"]));
    FCNC_booster->AddVariable("nbtag", &(BDT_params["nBtag"]));
    FCNC_booster->AddVariable("LeadJet_pt", &(BDT_params["LeadJet_pt"]));
    FCNC_booster->AddVariable("SubLeadJet_pt", &(BDT_params["SubLeadJet_pt"]));
    FCNC_booster->AddVariable("SubSubLeadJet_pt", &(BDT_params["SubSubLeadJet_pt"]));
    FCNC_booster->AddVariable("nElectron", &nElectron);
    FCNC_booster->AddVariable("MET_pt", &MET_pt);
    FCNC_booster->AddVariable("LeadBtag_pt", &(BDT_params["LeadBtag_pt"]));
    FCNC_booster->AddVariable("MT_LeadLep_MET", &(BDT_params["MT_LeadLep_MET"]));
    FCNC_booster->AddVariable("MT_SubLeadLep_MET", &(BDT_params["MT_SubLeadLep_MET"]));
    FCNC_booster->AddVariable("LeadLep_SubLeadLep_Mass", &LeadLep_SubLeadLep_Mass);
    FCNC_booster->BookMVA("BDT", "./helpers/BDT/BDT.xml");

    if (debug) {
        cout << "event: " << nt.event() << endl;
        cout << "Most_Forward_pt: " << BDT_params["Most_Forward_pt"] << endl;
        cout << "HT: " << BDT_params["HT"] << endl;
        cout << "LeadLep_eta: " << LeadLep_eta << endl;
        cout << "LeadLep_pt: " << LeadLep_pt << endl;
        cout << "LeadLep_dxy: " << LeadLep_dxy << endl;
        cout << "LeadLep_dz: " << LeadLep_dz << endl;
        cout << "SubLeadLep_pt: " << SubLeadLep_pt << endl;
        cout << "SubLeadLep_eta: " << SubLeadLep_eta << endl;
        cout << "SubLeadLep_dxy: " << SubLeadLep_dxy << endl;
        cout << "SubLeadLep_dz: " << SubLeadLep_dz << endl;
        cout << "nJet: " << BDT_params["nJets"] << endl;
        cout << "nBtag: " << BDT_params["nBtag"] << endl;
        cout << "LeadJet_pt: " << BDT_params["LeadJet_pt"] << endl;
        cout << "SubLeadJet_pt: " << BDT_params["SubLeadJet_pt"] << endl;
        cout << "SubSubLeadJet_pt: " << BDT_params["SubSubLeadJet_pt"] << endl;
        cout << "nElectron: " << nElectron << endl;
        cout << "MET_pt: " << MET_pt << endl;
        cout << "LeadBtag_pt: " << BDT_params["LeadBtag_pt"] << endl;
        cout << "MT_LeadLep_MET: " << BDT_params["MT_LeadLep_MET"] << endl;
        cout << "MT_SubLeadLep_MET: " << BDT_params["MT_SubLeadLep_MET"] << endl;
        cout << "LeadLep_SubLeadLep_Mass: " << LeadLep_SubLeadLep_Mass << endl;
    }
    Float_t booster_score = convert_tmva_to_prob(FCNC_booster->EvaluateMVA("BDT"));
    cout << "BDT eval time: " << duration_cast<milliseconds>(high_resolution_clock::now() - start_time).count() << endl;
    return booster_score;
}
