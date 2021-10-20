#include <iostream>
#include <string>
#include "TMVA/Reader.h"
#include "../booster.h"

using namespace std;

bool in_vector(std::vector<string> vec, string query) {
    for (string s : vec) {
        if (s == query) {
            return true;
        }
    }
    return false;
}
vector<string> BDT_features = {
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

void print_mean_std(vector<float> res) {
    float sum = 0;
    float count = 0;
    for (float r : res) {
        count++;
        sum += r;
    }
    float mean = sum/count;
    float var = 0;
    for(float r : res) {
        var += pow(r - mean, 2.0);
        //std::cout << var << std::endl;
    }
    var = var / count;
    std::cout << "avg difference = " << mean << " +/- " << pow(var, 0.5) << std::endl;
}

void compare_BDT_score(string path_to_xml, string path_to_csv, string test_csv, bool debug=false, int early_stop_round=-1) {
    BDT booster(path_to_xml, path_to_csv);
    ifstream fin;
    fin.open(test_csv);
    string line;
    char delimiter = ',';
    vector<string> features;
    vector<float> results;
    map<string, Float_t> param_map; 
    string record;
    fin >> line;
    stringstream head(line);
    while (getline(head, record, delimiter)){
        if (in_vector(BDT_features, record)) {
            features.push_back(record);
        }
    }
    int counter = 0;
    bool early_stop = (early_stop_round > 0);
    while (!fin.eof() && (!early_stop || counter < early_stop_round)){
        fin >> line;
        counter += 1;
        int feature_idx = 0;
        if (line.length() > 0) {
            stringstream s(line);
            for(string f : features){
                getline(s, record, delimiter);//get next entry in csv
                if (debug) {
                    cout << "feature: " << f << "\t value: " << record << endl;
                }
                param_map[f] = Float_t(stof(record));
            }
            getline(s, record, delimiter);
            float python_score = stof(record);
            booster.set_features(param_map);
            float TMVA_score = float(booster.get_score());
            float diff = python_score - TMVA_score;
            results.push_back(abs(diff));
            if (debug) {
                cout << "python: " << python_score << "\tTMVA: " << TMVA_score << "\tdifference: " << diff << endl;
            }
        }
        line.clear();
    }
    print_mean_std(results);
}
