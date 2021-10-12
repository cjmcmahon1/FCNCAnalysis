#include "booster.h"
#include "test_booster.h"
#include <string>

using namespace std;

void compare_score() {
    string hct_bdt = "BDT_HCT.xml";
    string hut_bdt = "BDT_HUT.xml";
    string hct_bins = "BDT_HCT_bins.csv";
    string hut_bins = "BDT_HUT_bins.csv";
    string python_results = "python_test_results.csv";
    cout << "HCT abs(TMVA Score - XGB Score):\n"; 
    compare_BDT_score(hct_bdt, hct_bins, python_results, false, -1);
    cout << "HUT abs(TMVA Score - XGB Score):\n"; 
    compare_BDT_score(hut_bdt, hut_bins, python_results, false, -1);
}
