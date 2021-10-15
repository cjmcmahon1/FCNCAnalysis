//Basic looper to categorize FCNC events
//Compile and run with root -l event_looper.C+
//March 2021
//Kaitlin Salyer

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <chrono>
#include <dirent.h>
#include <TChain.h>
#include <TChainElement.h>
#include <TLeaf.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TPad.h>
#include <THStack.h>
#include <TLegend.h>
#include <TFile.h>
#include <signal.h>
#include "./helpers/weights.h"
#include "./helpers/histogrammingClass.h"
//#include "./helpers/flip_weights.h"
#include "./helpers/BDT/booster.h"
#include "./helpers/BDT/make_baby.h"
#include "./helpers/tqdm.h"
#include "../../NanoTools/NanoCORE/SSSelections.h"
#include "../../NanoTools/NanoCORE/IsolationTools.h"
#include "../../NanoTools/NanoCORE/Nano.h"
#include "../../NanoTools/NanoCORE/Config.h"
#include "../../NanoTools/NanoCORE/Tools/goodrun.h"
#include "../../NanoTools/NanoCORE/Tools/dorky.h"
#include "../../NanoTools/NanoCORE/Tools/btagsf/BTagCalibrationStandalone.h"
#include "../../NanoTools/NanoCORE/Tools/btagsf/BTagCalibrationStandalone.cc"
#include "../misc/common_utils.h"
//#include "MVAUtils/BDT.h"

using namespace std;
using namespace std::chrono;

void print_debug (ofstream& outfile, int hyp_type, Leptons leptons, Jets jets, Jets bjets, int category) {
/*
    outfile << hyp_type << ",";
    for (auto lepton : leptons) {
        outfile << lepton.id() << "," << lepton.pt() << "," << lepton.eta() << "," << lepton.is_loose() << "," << lepton.is_tight() << ",";
    }
    outfile << jets.size() << "," << bjets.size() << ",";
*/
    //int category=3;
    //if (is_fake) category=1;
    //if (is_flip) category=2;
    outfile << category << ","; 

    int nloose = 0;
    int ntight = 0;
    for (auto lepton : leptons) {if (lepton.is_loose()) nloose++; if (lepton.is_tight()) ntight++;}
    outfile << nloose << ","  << ntight << "," << hyp_type << "," << jets.size() << "," << bjets.size() << std::endl;
    //for (auto lepton: leptons) {outfile << "\tL " << lepton.id() << "," << lepton.pt() << ", " << lepton.eta() << "," << lepton.is_loose() << ", " << lepton.is_tight() << ", " << lepton.genPartFlav() << ", " << lepton.isFake() << ", " << lepton.isFlip() << ", " << lepton.mcid() << ", " << lepton.id() << std::endl;}
    for (auto jet : jets) {outfile << "\tJ " << jet.pt() << ", " << jet.eta() << "," << jet.isBtag() << std::endl;}
    for (auto bjet : bjets) {if (bjet.pt()>40.) continue; outfile << "\tJ " << bjet.pt() << ", " << bjet.eta() << "," << bjet.isBtag() << std::endl;}
}

void print_sample_yields (ofstream& outfile, TString sampleName, float nFlip, float nFake, float nRare) {
    outfile << sampleName << " " << nFake << " " << nFlip << " " << nRare << endl;
}

void event_looper(TObjArray* list, TString title, TString options="", int nevts=-1, TString outputdir="outputs/"){
    bool debugPrints = 0;
    //*************************************************************************//
    //*************************** begin set options ***************************//
    bool quiet = options.Contains("quiet");
    bool evaluateBDT = options.Contains("evaluateBDT") && not options.Contains("noBDT");
    bool write_tree = options.Contains("writeTree");
    bool minPtFake18 = options.Contains("minPtFake18");
    bool new2016FRBins = options.Contains("new2016FRBins");
    bool doStitch = options.Contains("doStitch");
    bool doTruthFake = options.Contains("doTruthFake");
    bool doTruthFlip = options.Contains("doTruthFlip");
    bool isData = 0;
    bool doFlips = 0;
    bool doFakes = 0;
    bool make_BDT_fakes_babies = options.Contains("doBDTfakesbabies");
    bool make_BDT_flips_babies = options.Contains("doBDTflipsbabies");
    bool make_BDT_MC_babies  = options.Contains("doBDTMCbabies");
    int exclude = 0;
    bool isGamma = 0;
    bool truthfake = 0;
    bool skipmatching = 0;
    bool ignoreFakeFactor = 0;
    float scaleLumi=1.;
    int nbdtbins = 17;
    int nsrdisc = nbdtbins+1; // this is supposed to be 1 more than nbdtbins (we add in CRZ as a "bin")
    bool print_debug_file = options.Contains("printDebugFile");
    bool print_sample_file = options.Contains("printSampleFile");
    bool iterativeBTag = options.Contains("iterativefit");
    ofstream debug_file;
    ofstream sample_file;

/*
    bool STOP_REQUESTED = false;

    signal(SIGINT, [](int){
            cout << "SIGINT Caught, stopping after current event" << endl;
            STOP_REQUESTED=true;
            });
    STOP_REQUESTED=false;
*/

    int year = -1;
    TString extra("");
    if (options.Contains("Data2016")) {
        year = 2016;
        if (options.Contains("94x")) extra = "_94x";
    } else if (options.Contains("Data2017")) {
        year = 2017;
    } else if (options.Contains("Data2018")) {
        year = 2018;
    } else {
        cout << "Need to specify year!\n";
        assert(year > 0);
    }

    int year_for_output = year;
    float lumi = getLumi(year);
    //cout << lumi << endl;

    float min_pt_fake = minPtFake18 ? 18. : -1;

    if (options.Contains("FakeLumi2017")) {
        if (!quiet) std::cout << "Faking 2017 luminosity" << std::endl;
        lumi = getLumi(2017);
    } else if (options.Contains("FakeLumi2018")) {
        if (!quiet) std::cout << "Faking 2018 luminosity" << std::endl;
        lumi = getLumi(2018);
    }

    uint exactLeps=0;
    uint exactJets=0;
    uint exactBs=0;
    if(options.Contains("nleps2")){exactLeps=2;}
    if(options.Contains("nleps3")){exactLeps=3;}
    
    if(options.Contains("njets1")){exactJets=1;}
    if(options.Contains("njets2")){exactJets=2;}
    if(options.Contains("njets3")){exactJets=3;}
    if(options.Contains("njets4")){exactJets=4;}
    
    if(options.Contains("nbtags0")){exactBs=0;}
    if(options.Contains("nbtags1")){exactBs=1;}
    if(options.Contains("nbtags2")){exactBs=2;}

    bool doVarArg = 0;
    if(options.Contains("variations")){doVarArg = 1;}

    bool partialUnblind = false;
    /*
    if (options.Contains("partialUnblind")) {
        set_goodrun_file("unblindjson/allsnt.txt");
        if (year == 2017 or options.Contains("FakeLumi2017")) {
            lumi = 10.;
        }
        if (year == 2018 or options.Contains("FakeLumi2018")) {
            lumi = 15.415;
        }
        partialUnblind = true;
    }
    */

    // Note the blocks and else/ifs
    if (options.Contains("doData")) {
        if (!quiet) std::cout << "Doing data" << std::endl;
        isData = 1;
    }

    // --------
    if (options.Contains("doFakesMCUnw")) {
        if (!quiet) std::cout << "Doing fakesMC unweighted" << std::endl;
        doFakes = 1;
        ignoreFakeFactor = 1;
    }
    else if (options.Contains("doFakesUnw")) {
        if (!quiet) std::cout << "Doing fakes unweighted" << std::endl;
        isData = 1;
        doFakes = 1;
        ignoreFakeFactor = 1;
    }
    else if (options.Contains("doFakesMC")) {
        if (!quiet) std::cout << "Doing fakesMC" << std::endl;
        doFakes = 1;
    }
    else if (options.Contains("doFakes")) {
        if (!quiet) std::cout << "Doing fakes" << std::endl;
        isData = 1;
        doFakes = 1;
    }

    // --------
    if (options.Contains("doFlipsMC")) {
        if (!quiet) std::cout << "Doing flipsMC" << std::endl;
        doFlips = 1;
    }

    // --------
    if (options.Contains("doFlips")) {
        if (!quiet) std::cout << "Doing flips" << std::endl;
        isData = 1;
        doFlips = 1;
    }

    // --------
    if (options.Contains("doTruthFake")) {
        if (!quiet) std::cout << "Doing truthfake" << std::endl;
        truthfake = 1;
    }

    // --------
    if (options.Contains("doSkipMatching")) {
        if (!quiet) std::cout << "Skipping gen matching" << std::endl;
        skipmatching = 1;
    }

    // --------
    if (options.Contains("doXgamma")) {
        if (!quiet) std::cout << "Doing x+gamma" << std::endl;
        isGamma = 1;
    }

    if (!quiet) cout << "Running with lumi=" << lumi*scaleLumi << endl;

    TString chainTitle = title;
    const char* chainTitleCh = chainTitle.Data();
    if (print_debug_file) {
        debug_file.open(Form("debug/%s.log",chainTitleCh));
        debug_file << "run,lumi,evt"/*,cat,nloose,ntight,hyp_type,njets,nbjets"*/ << std::endl;
    }
    if (print_sample_file) {
        sample_file.open(Form("sampleYields_%s.txt",chainTitleCh));
        sample_file << "sample,nFakes,nFlips,nRares" << std::endl;
    }
    if (!quiet) std::cout << "Working on " << chainTitle << std::endl;

    bool isFakes = (chainTitle=="fakes") || (chainTitle=="fakes_mc");
    bool isFlips = (chainTitle=="flips") || (chainTitle=="flips_mc");
    bool isRares = (chainTitle=="rares");
    bool isSignal = (chainTitle=="signal_tch") || (chainTitle=="signal_tuh");
    bool istt = (chainTitle=="ttjets") || (chainTitle=="fakes_mc") || (chainTitle=="fakes_mc_unw");
    bool isttH = (chainTitle=="tth");
    bool istttt = (chainTitle=="tttt");
    bool isttVV = (chainTitle=="ttvv");
    bool isttW = (chainTitle=="ttw");
    bool isWW = (chainTitle=="ww");
    bool isttZ = (chainTitle=="ttz_m1-10") || (chainTitle=="ttz_m10");
    bool isWZ = (chainTitle=="wz");
    bool isXgamma = (chainTitle=="xg");

    //**************************** end set options ****************************//
    //*************************************************************************//

    //*************************************************************************//
    //************************** setup good run list **************************//
    std::string goodrun_path = "/home/users/ksalyer/FranksFCNC/ana/samples/goodRunList/";
    std::map<int, std::string> goodrun_file = {
        {2016, "goldenJson_2016rereco_36p46ifb.txt"},
        {2017, "Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1_snt.txt"},
        {2018, "goldenJson_2018_final_59p76ifb_snt.txt"} };
    set_goodrun_file( (goodrun_path+goodrun_file[year]).c_str() );

    //open trig eff SF file here so it's only opened once
    TFile* f_trigEff = new TFile("/home/users/ksalyer/FranksFCNC/ana/misc/year_run2/triggeffcymapsRA5_Run2_ALL.root");


    // //open b SF files here so they are opened only once
    auto startBOpening = high_resolution_clock::now();
    // if(!iterativeBTag){
    //efficiency files
    TFile* f_btag_eff_2016 = new TFile("./outputs/2016_bEff.root");
    TFile* f_btag_eff_2017 = new TFile("./outputs/2017_bEff.root");
    TFile* f_btag_eff_2018 = new TFile("./outputs/2018_bEff.root");

    TH1D* h_btag_eff_b_temp_2016    = (TH1D*) f_btag_eff_2016->Get("h2_BTaggingEff_med_Eff_b");
    TH1D* h_btag_eff_c_temp_2016    = (TH1D*) f_btag_eff_2016->Get("h2_BTaggingEff_med_Eff_c");
    TH1D* h_btag_eff_udsg_temp_2016 = (TH1D*) f_btag_eff_2016->Get("h2_BTaggingEff_med_Eff_udsg");

    TH1D* h_btag_eff_b_temp_2017    = (TH1D*) f_btag_eff_2017->Get("h2_BTaggingEff_med_Eff_b");
    TH1D* h_btag_eff_c_temp_2017    = (TH1D*) f_btag_eff_2017->Get("h2_BTaggingEff_med_Eff_c");
    TH1D* h_btag_eff_udsg_temp_2017 = (TH1D*) f_btag_eff_2017->Get("h2_BTaggingEff_med_Eff_udsg");

    TH1D* h_btag_eff_b_temp_2018    = (TH1D*) f_btag_eff_2018->Get("h2_BTaggingEff_med_Eff_b");
    TH1D* h_btag_eff_c_temp_2018    = (TH1D*) f_btag_eff_2018->Get("h2_BTaggingEff_med_Eff_c");
    TH1D* h_btag_eff_udsg_temp_2018 = (TH1D*) f_btag_eff_2018->Get("h2_BTaggingEff_med_Eff_udsg");

    TH1D* h_btag_eff_b_2016    = (TH1D*) h_btag_eff_b_temp_2016->Clone("h_btag_eff_b");
    TH1D* h_btag_eff_c_2016    = (TH1D*) h_btag_eff_c_temp_2016->Clone("h_btag_eff_c");
    TH1D* h_btag_eff_udsg_2016 = (TH1D*) h_btag_eff_udsg_temp_2016->Clone("h_btag_eff_udsg");

    TH1D* h_btag_eff_b_2017    = (TH1D*) h_btag_eff_b_temp_2017->Clone("h_btag_eff_b");
    TH1D* h_btag_eff_c_2017    = (TH1D*) h_btag_eff_c_temp_2017->Clone("h_btag_eff_c");
    TH1D* h_btag_eff_udsg_2017 = (TH1D*) h_btag_eff_udsg_temp_2017->Clone("h_btag_eff_udsg");

    TH1D* h_btag_eff_b_2018    = (TH1D*) h_btag_eff_b_temp_2018->Clone("h_btag_eff_b");
    TH1D* h_btag_eff_c_2018    = (TH1D*) h_btag_eff_c_temp_2018->Clone("h_btag_eff_c");
    TH1D* h_btag_eff_udsg_2018 = (TH1D*) h_btag_eff_udsg_temp_2018->Clone("h_btag_eff_udsg");

    std::map<std::string,TH1D*> eff2016 {
        {"btag", h_btag_eff_b_2016},
        {"ctag", h_btag_eff_c_2016},
        {"udsgtag", h_btag_eff_udsg_2016},
    };

    std::map<std::string,TH1D*> eff2017 {
        {"btag", h_btag_eff_b_2017},
        {"ctag", h_btag_eff_c_2017},
        {"udsgtag", h_btag_eff_udsg_2017},
    };

    std::map<std::string,TH1D*> eff2018 {
        {"btag", h_btag_eff_b_2018},
        {"ctag", h_btag_eff_c_2018},
        {"udsgtag", h_btag_eff_udsg_2018},
    };

    h_btag_eff_b_temp_2016->Delete();
    h_btag_eff_c_temp_2016->Delete();
    h_btag_eff_udsg_temp_2016->Delete();
    h_btag_eff_b_temp_2017->Delete();
    h_btag_eff_c_temp_2017->Delete();
    h_btag_eff_udsg_temp_2017->Delete();
    h_btag_eff_b_temp_2018->Delete();
    h_btag_eff_c_temp_2018->Delete();
    h_btag_eff_udsg_temp_2018->Delete();

    //csv paths
    string csv_path_2016 = "../../NanoTools/NanoCORE/Tools/btagsf/csv/DeepJet_2016LegacySF_WP_V1.csv";
    string csv_path_2017 = "../../NanoTools/NanoCORE/Tools/btagsf/csv/DeepFlavour_94XSF_WP_V3_B_F.csv";
    string csv_path_2018 = "../../NanoTools/NanoCORE/Tools/btagsf/csv/DeepJet_102XSF_WP_V1.csv";
    // std::cout << "set csv path, elapsed time since b SF start: " << duration_cast<seconds>(high_resolution_clock::now() - startBOpening).count() << endl;
    // CSV objects
    BTagCalibration deepjet_csv_2016 = BTagCalibration("csvv1", csv_path_2016);
    BTagCalibration deepjet_csv_2017 = BTagCalibration("csvv1", csv_path_2017);
    BTagCalibration deepjet_csv_2018 = BTagCalibration("csvv1", csv_path_2018);
    // std::cout << "opened csv file, elapsed time since b SF start: " << duration_cast<seconds>(high_resolution_clock::now() - startBOpening).count() << endl;
    // Medium reader
    BTagCalibrationReader deepjet_medium_reader_2016 = BTagCalibrationReader(BTagEntry::OP_MEDIUM, "central");
    BTagCalibrationReader deepjet_medium_reader_2017 = BTagCalibrationReader(BTagEntry::OP_MEDIUM, "central");
    BTagCalibrationReader deepjet_medium_reader_2018 = BTagCalibrationReader(BTagEntry::OP_MEDIUM, "central");
    // std::cout << "made calibration reader, elapsed time since b SF start: " << duration_cast<seconds>(high_resolution_clock::now() - startBOpening).count() << endl;
    //up
    BTagCalibrationReader up_reader_2016 = BTagCalibrationReader(BTagEntry::OP_MEDIUM, "up");
    BTagCalibrationReader up_reader_2017 = BTagCalibrationReader(BTagEntry::OP_MEDIUM, "up");
    BTagCalibrationReader up_reader_2018 = BTagCalibrationReader(BTagEntry::OP_MEDIUM, "up");
    //down
    BTagCalibrationReader down_reader_2016 = BTagCalibrationReader(BTagEntry::OP_MEDIUM, "down");
    BTagCalibrationReader down_reader_2017 = BTagCalibrationReader(BTagEntry::OP_MEDIUM, "down");
    BTagCalibrationReader down_reader_2018 = BTagCalibrationReader(BTagEntry::OP_MEDIUM, "down");

    //central
    deepjet_medium_reader_2016.load(deepjet_csv_2016, BTagEntry::FLAV_B, "comb");
    deepjet_medium_reader_2016.load(deepjet_csv_2016, BTagEntry::FLAV_C, "comb");
    deepjet_medium_reader_2016.load(deepjet_csv_2016, BTagEntry::FLAV_UDSG, "incl");
    deepjet_medium_reader_2017.load(deepjet_csv_2017, BTagEntry::FLAV_B, "comb");
    deepjet_medium_reader_2017.load(deepjet_csv_2017, BTagEntry::FLAV_C, "comb");
    deepjet_medium_reader_2017.load(deepjet_csv_2017, BTagEntry::FLAV_UDSG, "incl");
    deepjet_medium_reader_2018.load(deepjet_csv_2018, BTagEntry::FLAV_B, "comb");
    deepjet_medium_reader_2018.load(deepjet_csv_2018, BTagEntry::FLAV_C, "comb");
    deepjet_medium_reader_2018.load(deepjet_csv_2018, BTagEntry::FLAV_UDSG, "incl");
    //up
    up_reader_2016.load(deepjet_csv_2016, BTagEntry::FLAV_B, "comb");
    up_reader_2016.load(deepjet_csv_2016, BTagEntry::FLAV_C, "comb");
    up_reader_2016.load(deepjet_csv_2016, BTagEntry::FLAV_UDSG, "incl");
    up_reader_2017.load(deepjet_csv_2017, BTagEntry::FLAV_B, "comb");
    up_reader_2017.load(deepjet_csv_2017, BTagEntry::FLAV_C, "comb");
    up_reader_2017.load(deepjet_csv_2017, BTagEntry::FLAV_UDSG, "incl");
    up_reader_2018.load(deepjet_csv_2018, BTagEntry::FLAV_B, "comb");
    up_reader_2018.load(deepjet_csv_2018, BTagEntry::FLAV_C, "comb");
    up_reader_2018.load(deepjet_csv_2018, BTagEntry::FLAV_UDSG, "incl");
    //down
    down_reader_2016.load(deepjet_csv_2016, BTagEntry::FLAV_B, "comb");
    down_reader_2016.load(deepjet_csv_2016, BTagEntry::FLAV_C, "comb");
    down_reader_2016.load(deepjet_csv_2016, BTagEntry::FLAV_UDSG, "incl");
    down_reader_2017.load(deepjet_csv_2017, BTagEntry::FLAV_B, "comb");
    down_reader_2017.load(deepjet_csv_2017, BTagEntry::FLAV_C, "comb");
    down_reader_2017.load(deepjet_csv_2017, BTagEntry::FLAV_UDSG, "incl");
    down_reader_2018.load(deepjet_csv_2018, BTagEntry::FLAV_B, "comb");
    down_reader_2018.load(deepjet_csv_2018, BTagEntry::FLAV_C, "comb");
    down_reader_2018.load(deepjet_csv_2018, BTagEntry::FLAV_UDSG, "incl");
    // }else{
    if(iterativeBTag){
        //b discriminator reshaping
        string csv_path_2016 = "../../NanoTools/NanoCORE/Tools/btagsf/csv/DeepJet_2016LegacySF_V1.csv";
        string csv_path_2017 = "../../NanoTools/NanoCORE/Tools/btagsf/csv/DeepJet_94XSF_V4_B_F.csv";
        string csv_path_2018 = "../../NanoTools/NanoCORE/Tools/btagsf/csv/DeepJet_102XSF_V3.csv";
        cout << "got paths" << endl;

        BTagCalibration deepjet_csv_2016 = BTagCalibration("csvv1", csv_path_2016);
        BTagCalibration deepjet_csv_2017 = BTagCalibration("csvv1", csv_path_2017);
        BTagCalibration deepjet_csv_2018 = BTagCalibration("csvv1", csv_path_2018);
        cout << "loaded calibration" << endl;

        BTagCalibrationReader deepjet_medium_reader_2016 = BTagCalibrationReader(BTagEntry::OP_RESHAPING,
            "central",
            {"up_lf","down_lf",
            "up_hf","down_hf",
            "up_hfstats1","down_hfstats1",
            "up_hfstats2","down_hfstats2",
            "up_lfstats1","down_lfstats1",
            "up_lfstats2","down_lfstats2",
            "up_cferr1","down_cferr1",
            "up_cferr2","down_cferr2"}
            );

        BTagCalibrationReader deepjet_medium_reader_2017 = BTagCalibrationReader(BTagEntry::OP_RESHAPING,
            "central",
            {"up_lf","down_lf",
            "up_hf","down_hf",
            "up_hfstats1","down_hfstats1",
            "up_hfstats2","down_hfstats2",
            "up_lfstats1","down_lfstats1",
            "up_lfstats2","down_lfstats2",
            "up_cferr1","down_cferr1",
            "up_cferr2","down_cferr2"}
            );

        BTagCalibrationReader deepjet_medium_reader_2018 = BTagCalibrationReader(BTagEntry::OP_RESHAPING,
            "central",
            {"up_lf","down_lf",
            "up_hf","down_hf",
            "up_hfstats1","down_hfstats1",
            "up_hfstats2","down_hfstats2",
            "up_lfstats1","down_lfstats1",
            "up_lfstats2","down_lfstats2",
            "up_cferr1","down_cferr1",
            "up_cferr2","down_cferr2"}
            );
        cout << "loaded OP and syst unc" << endl;

        deepjet_medium_reader_2016.load(deepjet_csv_2016, BTagEntry::FLAV_B, "iterativefit");
        deepjet_medium_reader_2016.load(deepjet_csv_2016, BTagEntry::FLAV_C, "iterativefit");
        deepjet_medium_reader_2016.load(deepjet_csv_2016, BTagEntry::FLAV_UDSG, "iterativefit");
        deepjet_medium_reader_2017.load(deepjet_csv_2017, BTagEntry::FLAV_B, "iterativefit");
        deepjet_medium_reader_2017.load(deepjet_csv_2017, BTagEntry::FLAV_C, "iterativefit");
        deepjet_medium_reader_2017.load(deepjet_csv_2017, BTagEntry::FLAV_UDSG, "iterativefit");
        deepjet_medium_reader_2018.load(deepjet_csv_2018, BTagEntry::FLAV_B, "iterativefit");
        deepjet_medium_reader_2018.load(deepjet_csv_2018, BTagEntry::FLAV_C, "iterativefit");
        deepjet_medium_reader_2018.load(deepjet_csv_2018, BTagEntry::FLAV_UDSG, "iterativefit");
        cout << "loaded flavors for iterativefit" << endl;
    }

    auto outFileName = outputdir+"/"+chainTitle+"_"+TString(std::to_string(year).c_str())+"_hists.root";
    std::cout << "Will write histograms to " << outFileName.Data() << std::endl;

    //BDT constructor
    BDT hct_booster("./helpers/BDT/BDT_HCT.xml", "./helpers/BDT/BDT_HCT_bins.csv");
    BDT hut_booster("./helpers/BDT/BDT_HUT.xml", "./helpers/BDT/BDT_HUT_bins.csv");
    std::string tmp_yr_str = std::to_string(year);
    BDTBabyMaker bdt_fakes_baby;
    BDTBabyMaker bdt_flips_baby;
    //BDTBabyMaker bdt_MC_baby_training;
    //BDTBabyMaker bdt_MC_baby_testing;
    BDTBabyMaker bdt_MC_baby;
    std::string BDT_base_dir = "./helpers/BDT/babies/split_train_eval/";
    if (make_BDT_fakes_babies){
        char* output_baby_name = Form("%s/%s/data_driven/%s_fakes.root", BDT_base_dir.c_str(), tmp_yr_str.c_str(), chainTitleCh);
        bdt_fakes_baby.Initialize(output_baby_name);
        if (debugPrints) {
            cout << output_baby_name << endl;
        }
    }
    if (make_BDT_flips_babies){
        char* output_baby_name = Form("%s/%s/data_driven/%s_flips.root", BDT_base_dir.c_str(), tmp_yr_str.c_str(), chainTitleCh);
        bdt_flips_baby.Initialize(output_baby_name);
        if (debugPrints) {
            cout << output_baby_name << endl;
        }
    }
    if (make_BDT_MC_babies){
        //char* output_training_baby_name = Form("%s/%s/MC/training/%s.root", BDT_base_dir.c_str(), tmp_yr_str.c_str(), chainTitleCh);
        //char* output_testing_baby_name = Form("%s/%s/MC/testing/%s.root", BDT_base_dir.c_str(), tmp_yr_str.c_str(), chainTitleCh);
        //bdt_MC_baby_training.Initialize(output_training_baby_name);
        //bdt_MC_baby_testing.Initialize(output_testing_baby_name);
        char* output_baby_name = Form("%s/%s/MC/%s.root", BDT_base_dir.c_str(), tmp_yr_str.c_str(), chainTitleCh);
        bdt_MC_baby.Initialize(output_baby_name);
    }

    HistContainer hists;
    hists.loadHists(chainTitleCh, hct_booster.get_BDT_bins(year), hut_booster.get_BDT_bins(year));

    //nEvents in chain
    unsigned int nEventsTotal = 0;
    //unsigned int nEventsChain = chain->GetEntries();
    //std::cout << "Events in Chain: " << nEventsChain << endl;

    //Set up iterator
    TObjArray *listOfFiles = list;
    TIter fileIter(listOfFiles);
    TFile *currentFile = 0;

    //delete chain;

    //Number of selected events
    //int nSelected = 0;

    tqdm bar;
    // bar.set_theme_braille();
    
    auto start = high_resolution_clock::now();

    //File Loop
    // int nLooseFlips=0;
    // int nTightFlips=0;

    // int nDilep_ss=0;
    // int nDilep_fake_ss=0;
    // int nDilep_flip_ss=0;
    // int n1lep_ss=0;

    // int nDilep_sf=0;
    // int nDilep_fake_sf=0;
    // int nDilep_flip_sf=0;
    // int n1lep_sf=0;
    while ( (currentFile = (TFile*)fileIter.Next()) ){
        TFile *file = new TFile(currentFile->GetName());
        TString filename = currentFile->GetName();
        TTree *tree = (TTree*)file->Get("Events");
        nt.Init(tree);
        if (debugPrints){std::cout << "working on file " << filename << endl;}
        if (debugPrints){std::cout << "elapsed time since start: " << duration_cast<seconds>(high_resolution_clock::now() - start).count() << endl;}
        if (debugPrints){std::cout << "elapsed time since b SF start: " << duration_cast<seconds>(high_resolution_clock::now() - startBOpening).count() << endl;}

        float nFakes = 0.;
        float nFlips = 0.;
        float nRares = 0.;

        for ( unsigned int counter = 0; counter < tree->GetEntries(); counter++ ){ 
        //for ( unsigned int counter = 0; counter < 1000; counter++ ){ //for testing only!!
            //if (counter%100000==0){std::cout << "processed " << counter << " events and elapsed time since start: " << duration_cast<seconds>(high_resolution_clock::now() - start).count() << endl;}

            nt.GetEntry(counter);
            if ( nevts>0 && nEventsTotal>=nevts ) break;
            ++nEventsTotal;

            //if (!quiet) bar.progress(nEventsTotal, nEventsChain);

            // filter lumi blocks not in the good run list
            if ( isData && !goodrun( nt.run(), nt.luminosityBlock() ) ) continue;

            // remove duplicates
            if ( isData ) {
                duplicate_removal::DorkyEventIdentifier id( nt.run(), nt.event(), nt.luminosityBlock() );
                if ( duplicate_removal::is_duplicate(id) ) continue;
            }


            // cutflow counter
            int cutflow_counter=1;
            hists.fill1d("cutflow","br",chainTitleCh,1,1);
            cutflow_counter++;

            //if(nt.event()!=313915097){continue;}

            //MET CUT
            //if (nt.MET_pt() < 50.){continue;}
            //jes central
            //jes up
            // if (nt.MET_T1_pt_jesTotalUp() < 50.){continue;}
            //jes down
            // if (nt.MET_T1_pt_jesTotalDown() < 50.){continue;}
            hists.fill1d("cutflow","br",chainTitleCh,2,1);
            if (debugPrints){std::cout << "passed MET cut for event " << nt.event() << ": " << nt.MET_pt() << endl;}
            if (debugPrints){std::cout << "elapsed time since start: " << duration_cast<seconds>(high_resolution_clock::now() - start).count() << endl;}
            if (debugPrints){std::cout << "elapsed time since b SF start: " << duration_cast<seconds>(high_resolution_clock::now() - startBOpening).count() << endl;}

            //get event weight based on sample!
            // cout << nt.event() << endl;
            double weight = 1.;
            if (!isData){
                weight = getEventWeight( file->GetName(), chainTitle.Data(), nt.year());
                //cout << "scale1fb: " << weight << endl;
                //double weight = 1.;
                double genWeight = nt.Generator_weight();
                //cout << "genWeight: " << genWeight << endl;
                weight = (weight*genWeight*lumi)/(abs(genWeight));
                //cout << "weight: " << weight << endl;
                //cout << "lumi: " << lumi << endl;
            }else {weight = 1.;}
            //cout << weight << endl;
            // double evtWeight = weight; 
            if (debugPrints){std::cout << "passed eventWeight for event " << nt.event() << endl;}
            if (debugPrints){std::cout << "elapsed time since start: " << duration_cast<seconds>(high_resolution_clock::now() - start).count() << endl;}
            if (debugPrints){std::cout << "elapsed time since b SF start: " << duration_cast<seconds>(high_resolution_clock::now() - startBOpening).count() << endl;}


            // load the leptons
            Leptons leptons = getLeptons();
            Leptons tight_leptons = getTightLeptons();
            Leptons loose_leptons = getLooseLeptons();
            unsigned int nleps_tight = tight_leptons.size();
            unsigned int nleps_loose = loose_leptons.size();

            // for (auto lep : leptons){cout << "lepton " << lep.id() << " " << lep.eta() << " " << lep.pt() << endl;}
            // for (auto lep : tight_leptons){cout << "tight lepton " << lep.id() << " " << lep.eta() << " " << lep.pt() << endl;}
            // for (auto lep : loose_leptons){cout << "loose lepton " << lep.id() << " " << lep.eta() << " " << lep.pt() << endl;}
            // cout << "*************" << endl;
            
            if (!quiet) {
                std::cout << "Event has " << nleps_loose << " loose leptons and "
                          << nleps_tight << " tight leptons." << std::endl;
                for (auto lep: leptons) {
                    std::cout << "pdgId: " << lep.id() << ", pt: " << lep.pt() << ", eta: " << lep.eta()
                              << ", loose: " << lep.is_loose() << ", tight: " << lep.is_tight() << std::endl;
                }
            }
            if (debugPrints){std::cout << "loaded leptons for event " << nt.event() << endl;}
            if (debugPrints){std::cout << "elapsed time since start: " << duration_cast<seconds>(high_resolution_clock::now() - start).count() << endl;}
            if (debugPrints){std::cout << "elapsed time since b SF start: " << duration_cast<seconds>(high_resolution_clock::now() - startBOpening).count() << endl;}


            // for (auto lepton: leptons) {cout << "\tL " << lepton.id() << "," << lepton.pt() << ", " << lepton.eta() << "," << lepton.is_loose() << ", " << lepton.is_tight() << ", " << lepton.genPartFlav() << ", " << lepton.isFake() << ", " << lepton.isFlip() << ", " << lepton.mcid() << ", " << lepton.id() << ", " << lepton.idlevel() << std::endl;}
            // for (auto lepton: loose_leptons) {cout << "\tL " << lepton.id() << "," << lepton.pt() << ", " << lepton.eta() << "," << lepton.is_loose() << ", " << lepton.is_tight() << ", " << lepton.genPartFlav() << ", " << lepton.isFake() << ", " << lepton.isFlip() << ", " << lepton.mcid() << ", " << lepton.id() << std::endl;}
            // for (auto lepton: tight_leptons) {cout << "\tL " << lepton.id() << "," << lepton.pt() << ", " << lepton.eta() << "," << lepton.is_loose() << ", " << lepton.is_tight() << ", " << lepton.genPartFlav() << ", " << lepton.isFake() << ", " << lepton.isFlip() << ", " << lepton.mcid() << ", " << lepton.id() << std::endl;}
            // if(nleps_tight!=2){continue;}
            // // cout << "event has TT: " << nt.event() << endl;
            // if(tight_leptons[0].charge()!=tight_leptons[1].charge()){continue;}
            // // cout << "event has SS: " << nt.event() << endl;

            // if (print_debug_file) {
            //     debug_file << nt.run() << "," << nt.luminosityBlock() << "," << nt.event() << endl;
            //     continue;
            // }

            // let's first figure out what kind of lepton hypothesis we have, by priority: 3L, TTL, TT, OS, TL, LL
            std::pair<int,Leptons> best_hyp_info = getBestHypFCNC(leptons,!quiet);
            int best_hyp_type = best_hyp_info.first;

            // if (best_hyp_type < 0 && print_debug_file) {
            //     debug_file << nt.run() << "," << nt.luminosityBlock() << "," << nt.event() << ",";
            //     debug_file << " ," << loose_leptons.size() << "," << tight_leptons.size();
            //     debug_file << std::endl;
            //     for (auto lepton: loose_leptons) {debug_file << "\tL " << lepton.id() << "," << lepton.pt() << ", " << lepton.eta() << "," << lepton.is_loose() << ", " << lepton.is_tight() << std::endl;}
            //     debug_file << std::endl;
            //     continue;
            // }
            
            if (best_hyp_type < 0) continue;
            else{hists.fill1d("cutflow","br",chainTitleCh,5,weight);}
            
            if(!quiet) {std::cout << "best_hyp_type: " << best_hyp_type << std::endl;}
            Leptons best_hyp = best_hyp_info.second;
            if (!quiet) std::cout << "best hyp type: " << best_hyp_type << std::endl;

            sort(best_hyp.begin(),best_hyp.end(),lepsort); // sort hyp leptons by pt (may already be done in getBestHyp)
            Lepton leading_lep = best_hyp[0];
            Lepton trailing_lep = best_hyp[1];
            Lepton third_lep;
            if (best_hyp.size()>2) third_lep = best_hyp[2];


            if (!quiet) {
                std::cout << "best hyp type: " << best_hyp_type << "; lepton ids: " << leading_lep.id()
                                                                << ", " << trailing_lep.id() << std::endl;
            }
            if (debugPrints){std::cout << "got best hyp for event " << nt.event() << endl;}
            if (debugPrints){std::cout << "elapsed time since start: " << duration_cast<seconds>(high_resolution_clock::now() - start).count() << endl;}
            if (debugPrints){std::cout << "elapsed time since b SF start: " << duration_cast<seconds>(high_resolution_clock::now() - startBOpening).count() << endl;}

            // now let's check if there are additional requirements we need to make
            bool is_fake = false;
            bool is_flip = false;
            bool is_os = false;
            bool is_rare = false;
            int nfakes=0;
            int nflips=0;
            if(!isData){
                for ( auto lep : best_hyp ) {if (lep.isFake()) nfakes++;}
                for ( auto lep : best_hyp ) {
                    if (lep.isFlip() && lep.absid()==11){
                        nflips++;
                        // if(lep.idlevel()==SS::IDLevel::IDtight){nTightFlips++;}
                        // else if(lep.idlevel()==SS::IDLevel::IDfakable){nLooseFlips++;}
                    }
                }
                if (best_hyp.size() == 2 && best_hyp[0].charge()*best_hyp[1].charge()<0) is_os=true;
                if (nfakes>0 && !is_os) is_fake=true;
                if (best_hyp.size() == 2 && nflips==1) is_flip=true;
                if (!is_fake && (best_hyp.size()>2 || (best_hyp.size()==2 && !is_flip && best_hyp[0].charge()*best_hyp[1].charge()>0))) is_rare = true;
            }
            int category=4;
            if (is_flip) category=2;
            if (is_rare) category=3;
            if (is_fake) category=1;
            // cout << category << endl;

            //count dilep vs 1lep events
            // int ngenleps=0;
            // for(int i=0; i<nt.GenPart_pdgId().size(); i++){
            //     if(abs(nt.GenPart_pdgId()[i])==11||abs(nt.GenPart_pdgId()[i])==13){
            //         if(abs(nt.GenPart_pdgId()[nt.GenPart_genPartIdxMother()[i]])==24){
            //             ngenleps++;
            //         }
            //     }
            // }
            // if(best_hyp_type==4){
            //     if(ngenleps==2){nDilep_ss++;}
            //     if(ngenleps==1){n1lep_ss++;}
            //     if(ngenleps==2&&category==1){nDilep_fake_ss++;}
            //     if(ngenleps==2&&category==2){nDilep_flip_ss++;}
            // }else if(best_hyp_type==6){
            //     if(ngenleps==2){
            //         nDilep_sf++;
            //         for(auto lep : best_hyp){
            //             if(lep.idlevel()==SS::IDLevel::IDfakable){
            //                 cout << nDilep_sf << " loose lep reco id: " << lep.id() << " and gen id: " << lep.mcid() << " it is category " << category << endl;
            //                 if(lep.id()==-1*lep.mcid()){cout << "     above is a flip"<< endl;}
            //             }
            //         }
            //     }
            //     if(ngenleps==1){n1lep_sf++;}
            //     if(ngenleps==2&&category==1){
            //         nDilep_fake_sf++;
            //         // cout << "sf dilep fake: " << nt.event() << endl;
            //     }
            //     if(ngenleps==2&&category==2){
            //         nDilep_flip_sf++;
            //         // cout << "sf dilep flip: " << nt.event() << endl;
            //     }
            // }
            if (debugPrints){std::cout << "got category for event " << nt.event() << endl;}
            if (debugPrints){std::cout << "elapsed time since start: " << duration_cast<seconds>(high_resolution_clock::now() - start).count() << endl;}
            if (debugPrints){std::cout << "elapsed time since b SF start: " << duration_cast<seconds>(high_resolution_clock::now() - startBOpening).count() << endl;}


            //for fake validation: gen matching for best_hyp leptons
            bool isVR_CR_fake = 0;
            bool isVR_SR_fake = 0;
            bool isEE = 0;
            bool isMM = 0;
            bool isEM = 0;
            bool isME = 0;
            bool isEFake = 0;
            bool isMFake = 0;
            //for flip validation: gen matching for best_hyp leptons
            bool isVR_CR_flip = 0;
            bool isVR_SR_flip = 0;
            bool isEE_flip = 0;
            bool isEM_flip = 0;
            if (!isData){
                vector<int> nPromptTight;
                vector<int> nNonPromptTight;
                vector<int> nPromptLoose;
                vector<int> nNonPromptLoose;
                for ( auto lep : loose_leptons){
                    if (lep.idlevel()==SS::IDLevel::IDtight && (lep.genPartFlav()==1 || lep.genPartFlav()==15)){nPromptTight.push_back(lep.absid());}
                    if (lep.idlevel()==SS::IDLevel::IDtight && !(lep.genPartFlav()==1 || lep.genPartFlav()==15)){nNonPromptTight.push_back(lep.absid());}
                    if (lep.idlevel()==SS::IDLevel::IDfakable && (lep.genPartFlav()==1 || lep.genPartFlav()==15)){nPromptLoose.push_back(lep.absid());}
                    if (lep.idlevel()==SS::IDLevel::IDfakable && !(lep.genPartFlav()==1 || lep.genPartFlav()==15)){nNonPromptLoose.push_back(lep.absid());}
                }
                vector<int> nTightNotFlip;
                vector<int> nTightFlip;
                //cout << "******************" << endl;
                for ( auto lep : tight_leptons ){
                    //cout << lep.id() << " " << lep.mcid() << " " << lep.pt() << " " << lep.eta() << " " << lep.isFlip() << endl;
                    if (lep.idlevel()==SS::IDLevel::IDtight && (lep.genPartFlav()==1||lep.genPartFlav()==15) && !lep.isFlip()){nTightNotFlip.push_back(lep.id());}
                    if (lep.idlevel()==SS::IDLevel::IDtight && (lep.genPartFlav()==1||lep.genPartFlav()==15) && lep.absid()==11 && lep.isFlip()){nTightFlip.push_back(lep.id());}
                }

                if(best_hyp.size()==2){
                    //fake categorization
                    if ((nPromptTight.size()==1&&nNonPromptTight.size()==1)||(nNonPromptTight.size()==2)){
                        isVR_SR_fake=1;
                        if (nPromptTight.size()==1){
                            if (nPromptTight[0]==11 && nNonPromptTight[0]==11){isEE=1;}
                            else if (nPromptTight[0]==11 && nNonPromptTight[0]==13){isEM=1;}
                            else if (nPromptTight[0]==13 && nNonPromptTight[0]==11){isME=1;}
                            else if (nPromptTight[0]==13 && nNonPromptTight[0]==13){isMM=1;}
                        }else{
                            if (nNonPromptTight[0]==11 && nNonPromptTight[1]==11){isEE=1;}
                            else if (nNonPromptTight[0]==11 && nNonPromptTight[1]==13){isEM=1;}
                            else if (nNonPromptTight[0]==13 && nNonPromptTight[1]==11){isME=1;}
                            else if (nNonPromptTight[0]==13 && nNonPromptTight[1]==13){isMM=1;}
                        }
                    }
                    if ((nPromptTight.size()==1&&nNonPromptLoose.size()==1)||(nNonPromptLoose.size()==2)){
                        isVR_CR_fake=1;
                        if (nPromptTight.size()==1){
                            if (nPromptTight[0]==11 && nNonPromptLoose[0]==11){isEE=1;}
                            else if (nPromptTight[0]==11 && nNonPromptLoose[0]==13){isEM=1;}
                            else if (nPromptTight[0]==13 && nNonPromptLoose[0]==11){isME=1;}
                            else if (nPromptTight[0]==13 && nNonPromptLoose[0]==13){isMM=1;}
                        }else{
                            if (nNonPromptLoose[0]==11 && nNonPromptLoose[1]==11){isEE=1;}
                            else if (nNonPromptLoose[0]==11 && nNonPromptLoose[1]==13){isEM=1;}
                            else if (nNonPromptLoose[0]==13 && nNonPromptLoose[1]==11){isME=1;}
                            else if (nNonPromptLoose[0]==13 && nNonPromptLoose[1]==13){isMM=1;}
                        }
                    }
                    //flip categorization
                    if(nTightNotFlip.size()==2 && (nTightNotFlip[0]/abs(nTightNotFlip[0]))!=(nTightNotFlip[1]/abs(nTightNotFlip[1]))){
                        //check if there's an electron in the CR event
                        if(abs(nTightNotFlip[0])==11){
                            if(abs(nTightNotFlip[1])==11){
                                isVR_CR_flip=1;
                                isEE_flip=1;
                            }else{
                                isVR_CR_flip=1;
                                isEM_flip=1;
                            }
                        }else if(abs(nTightNotFlip[1])==11){
                            isVR_CR_flip=1;
                            isEM_flip=1;
                        }
                    }
                    //can't enter VR SR without containing an electron
                    //because nTightFlip isn't filled for muons
                    if(nTightNotFlip.size()==1 && nTightFlip.size()==1 && (nTightNotFlip[0]/abs(nTightNotFlip[0]))==(nTightFlip[0]/abs(nTightFlip[0]))){
                        isVR_SR_flip=1;
                        //cout << "found a VR_SR event" << endl;
                        if(abs(nTightNotFlip[0])==11&&abs(nTightFlip[0])==11){isEE_flip=1;}
                        else{isEM_flip=1;}
                    }
                }
                if(best_hyp.size()==3){
                    isVR_CR_flip = 0;
                    isVR_SR_flip = 0;
                    isEE_flip = 0;
                    isEM_flip = 0;
                    
                    if ((nPromptTight.size()==2&&nNonPromptTight.size()==1)||(nPromptTight.size()==1&&nNonPromptTight.size()==2)){
                        isVR_SR_fake=1;
                        if (nNonPromptTight.size()==1){
                            if (nNonPromptTight[0]==11){isEFake=1;}
                            else if (nNonPromptTight[0]==13){isMFake=1;}
                        }else{
                            if (nNonPromptTight[0]==11 && nNonPromptTight[1]==11){isEE=1;}
                            else if (nNonPromptTight[0]==11 && nNonPromptTight[1]==13){isEM=1;}
                            else if (nNonPromptTight[0]==13 && nNonPromptTight[1]==11){isME=1;}
                            else if (nNonPromptTight[0]==13 && nNonPromptTight[1]==13){isMM=1;}
                        }
                    }
                    if ((nPromptTight.size()==2&&nNonPromptLoose.size()==1)||(nPromptTight.size()==1&&nNonPromptLoose.size()==2)){
                        isVR_CR_fake=1;
                        if (nNonPromptLoose.size()==1){
                            if (nNonPromptLoose[0]==11){isEFake=1;}
                            else if (nNonPromptLoose[0]==13){isMFake=1;}
                        }else{
                            if (nNonPromptLoose[0]==11 && nNonPromptLoose[1]==11){isEE=1;}
                            else if (nNonPromptLoose[0]==11 && nNonPromptLoose[1]==13){isEM=1;}
                            else if (nNonPromptLoose[0]==13 && nNonPromptLoose[1]==11){isME=1;}
                            else if (nNonPromptLoose[0]==13 && nNonPromptLoose[1]==13){isMM=1;}
                        }
                    }
                }
            }
            if (debugPrints){std::cout << "passed validation categorization for event " << nt.event() << endl;}
            if (debugPrints){std::cout << "elapsed time since start: " << duration_cast<seconds>(high_resolution_clock::now() - start).count() << endl;}
            if (debugPrints){std::cout << "elapsed time since b SF start: " << duration_cast<seconds>(high_resolution_clock::now() - startBOpening).count() << endl;}

            // get jets and bjets
            // cout << "before get jets" << endl;
            //getJets parameters: Leptons &leps,float min_jet_pt=40., float min_bjet_pt=25.
            //central
            std:pair<Jets, Jets> good_jets_and_bjets = getJets(best_hyp,30.,25.,0);
            // //jes down
            //std:pair<Jets, Jets> good_jets_and_bjets = getJets(best_hyp,30.,25.,-1);
            // //jes up
            //std:pair<Jets, Jets> good_jets_and_bjets = getJets(best_hyp,30.,25.,1);
            // std:pair<Jets, Jets> good_jets_and_bjets = getJets(best_hyp);
            Jets good_jets = good_jets_and_bjets.first;
            Jets good_bjets = good_jets_and_bjets.second;
            int njets = good_jets.size();
            int nbjets = good_bjets.size();
            sort(good_jets.begin(),good_jets.end(),jetptsort);
            sort(good_bjets.begin(),good_bjets.end(),jetptsort);
            float ht = 0.;
            for (Jet j : good_jets){
                ht += j.pt();
            }
            
            //BJET CUT
            // if(nbjets<1){continue;}

            if (debugPrints){std::cout << "loaded jets for event " << nt.event() << endl;}
            if (debugPrints){std::cout << "elapsed time since start: " << duration_cast<seconds>(high_resolution_clock::now() - start).count() << endl;}
            if (debugPrints){std::cout << "elapsed time since b SF start: " << duration_cast<seconds>(high_resolution_clock::now() - startBOpening).count() << endl;}




            //apply SFs to MC events
            // cout << nt.event() << endl;
            // cout << weight << endl;
            if (!isData){
                //apply lepton SFs to hypothesis leptons
                for ( auto lep : best_hyp ) {weight = weight * leptonScaleFactor(nt.year(), lep.id(), lep.pt(), lep.eta(), ht);}
                // cout << "after lepton SF: " << weight << endl;

                //apply PU weight
                weight = weight * nt.puWeight();
                // cout << "after PU weight: " << weight << endl;

                //trigger scale factors to dilepton events
                weight = weight * triggerScaleFactor(f_trigEff, nt.year(), best_hyp[0].id(), best_hyp[1].id(), best_hyp[0].pt(), best_hyp[1].pt(), best_hyp[0].eta(), best_hyp[1].eta(), ht, 0);
                // cout << triggerScaleFactor(f_trigEff, nt.year(), best_hyp[0].id(), best_hyp[1].id(), best_hyp[0].pt(), best_hyp[1].pt(), best_hyp[0].eta(), best_hyp[1].eta(), ht, 0) << endl;
                // cout << "after trigger SF: " << weight << endl;

                // //applying b-tag SFs
                if(!iterativeBTag){
                    if (nt.year() == 2016){weight = weight * getBSF(nt.year(),good_jets,good_bjets,eff2016,deepjet_medium_reader_2016);}
                    else if (nt.year() == 2017){weight = weight * getBSF(nt.year(),good_jets,good_bjets,eff2017,deepjet_medium_reader_2017);}
                    else if (nt.year() == 2018){weight = weight * getBSF(nt.year(),good_jets,good_bjets,eff2018,deepjet_medium_reader_2018);}
                }else{
                    if (nt.year() == 2016){weight = weight * getIterativeBSF(nt.year(),good_jets,good_bjets,deepjet_medium_reader_2016,"central");}
                    else if (nt.year() == 2017){weight = weight * getIterativeBSF(nt.year(),good_jets,good_bjets,deepjet_medium_reader_2017,"central");}
                    else if (nt.year() == 2018){weight = weight * getIterativeBSF(nt.year(),good_jets,good_bjets,deepjet_medium_reader_2018,"central");}
                }
                // cout << "after b tag SF: " << weight << endl;
            }
            if(isnan(weight)||isinf(weight)){
                cout << "found a wrong weight!!!!!" << endl;
                cout << nt.event() << endl;
                cout << "final weight: " << weight << endl;
                cout << "******************" << endl;
            }
            // if ((category == 1 && chainTitle=="fakes_mc")||
            //     (category == 2 && chainTitle=="flips_mc")||
            //     (category == 3 && chainTitle=="rares")||
            //     // (category == 4 && chainTitle!="fakes_mc" && chainTitle!="flips_mc" && chainTitle!="rares" && !(isSignal||isData))||
            //     (isSignal||isData)){
                
            //     hists.fill1d("cutflow","br",chainTitleCh,3,weight);
            //     if(best_hyp_type==4){hists.fill1d("cutflow","br",chainTitleCh,6,weight);}
            //     else{continue;}
            // }
            if (debugPrints){std::cout << "passed sfs for event " << nt.event() << endl;}
            if (debugPrints){std::cout << "elapsed time since start: " << duration_cast<seconds>(high_resolution_clock::now() - start).count() << endl;}
            if (debugPrints){std::cout << "elapsed time since b SF start: " << duration_cast<seconds>(high_resolution_clock::now() - startBOpening).count() << endl;}


            if (print_debug_file) {
                debug_file << nt.run() << "," << nt.luminosityBlock() << "," << nt.event() << ",";
                print_debug(debug_file,best_hyp_type,best_hyp,good_jets,good_bjets,category);
                //debug_file << std::endl;
            }
            
            // if there isn't a good lepton hypothesis
            if (best_hyp_type<0) continue;
            //else {hists.fill1d("cutflow","br",chainTitleCh,2,weight);}


            if (!quiet) {
                std::cout << "Event has " << good_jets.size() << " jets and "
                          << good_bjets.size() << " b-tagged jets." << std::endl;
            }

            //NJETS CUT
            if (best_hyp.size()==2&&njets < 2) continue;
            if (best_hyp.size()>2&&njets < 1) continue;
            // if (best_hyp.size()>2&&nbjets<2&&njets < 1) continue;
            // if (best_hyp.size()>2&&nbjets>=2&&njets < 2) continue;
            if ((category == 1 && chainTitle=="fakes_mc")||
                (category == 2 && chainTitle=="flips_mc")||
                (category == 3 && chainTitle=="rares")||
                // (category == 4 && chainTitle!="fakes_mc" && chainTitle!="flips_mc" && chainTitle!="rares" && !(isSignal||isData))||
                (isSignal||isData)){
                hists.fill1d("cutflow","br",chainTitleCh,7,weight);
            }

            if (!quiet) {
                std::cout << "leading jet pt: " << good_jets[0].pt() << std::endl;
            }

            if (doFakes && !isData && doTruthFake && !is_fake) continue;
            if (doFlips && !isData && doTruthFlip && !is_flip) continue;
            if (nt.year()==2016){
                if(!(
                     nt.HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL()||
                     nt.HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ()||
                     nt.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ()||
                     nt.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL()||

                     nt.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL()||
                     nt.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ()||
                     nt.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ()//||
                     // nt.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ()||
                     // nt.HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ()
                     )) {continue;}
            }else if (nt.year()==2017){
                if(!(nt.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL()||
                     nt.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8()||
                     nt.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL()||
                     nt.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ()||
                     nt.HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ())) {continue;}
            }else if (nt.year()==2018){
                if(!(nt.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8()||
                     nt.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL()||
                     nt.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ()||
                     nt.HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ())) {continue;}
            }
            if ((category == 1 && chainTitle=="fakes_mc")||
                (category == 2 && chainTitle=="flips_mc")||
                (category == 3 && chainTitle=="rares")||
                // (category == 4 && chainTitle!="fakes_mc" && chainTitle!="flips_mc" && chainTitle!="rares" && !(isSignal||isData))||
                (isSignal||isData)){
                hists.fill1d("cutflow","br",chainTitleCh,8,weight);
            }
            if (debugPrints){std::cout << "passed triggers for event " << nt.event() << endl;}
            if (debugPrints){std::cout << "elapsed time since start: " << duration_cast<seconds>(high_resolution_clock::now() - start).count() << endl;}
            if (debugPrints){std::cout << "elapsed time since b SF start: " << duration_cast<seconds>(high_resolution_clock::now() - startBOpening).count() << endl;}

            if (nt.year()==2016){
                if(!(nt.Flag_globalSuperTightHalo2016Filter()&&
                     nt.Flag_HBHENoiseFilter()&&
                     nt.Flag_HBHENoiseIsoFilter()&&
                     nt.Flag_EcalDeadCellTriggerPrimitiveFilter()&&
                     nt.Flag_BadPFMuonFilter()&&
                     nt.Flag_goodVertices())) {continue;}
            }else if (nt.year()==2017){
                if(!(nt.Flag_globalSuperTightHalo2016Filter()&&
                     nt.Flag_HBHENoiseFilter()&&
                     nt.Flag_HBHENoiseIsoFilter()&&
                     nt.Flag_EcalDeadCellTriggerPrimitiveFilter()&&
                     nt.Flag_BadPFMuonFilter()&&
                     nt.Flag_goodVertices()&&
                     nt.Flag_ecalBadCalibFilterV2())) {continue;}
            }else if (nt.year()==2018){
                if(!(nt.Flag_globalSuperTightHalo2016Filter()&&
                     nt.Flag_HBHENoiseFilter()&&
                     nt.Flag_HBHENoiseIsoFilter()&&
                     nt.Flag_EcalDeadCellTriggerPrimitiveFilter()&&
                     nt.Flag_BadPFMuonFilter()&&
                     nt.Flag_goodVertices()&&
                     nt.Flag_ecalBadCalibFilterV2())) {continue;}
            }

            if ((category == 1 && chainTitle=="fakes_mc")||
                (category == 2 && chainTitle=="flips_mc")||
                (category == 3 && chainTitle=="rares")||
                // (category == 4 && chainTitle!="fakes_mc" && chainTitle!="flips_mc" && chainTitle!="rares" && !(isSignal||isData))||
                (isSignal||isData)){
                hists.fill1d("cutflow","br",chainTitleCh,9,weight);
            }
            if (debugPrints){std::cout << "passed filters for event " << nt.event() << endl;}
            if (debugPrints){std::cout << "elapsed time since start: " << duration_cast<seconds>(high_resolution_clock::now() - start).count() << endl;}
            if (debugPrints){std::cout << "elapsed time since b SF start: " << duration_cast<seconds>(high_resolution_clock::now() - startBOpening).count() << endl;}







            //CUTS FOR FAKE SYSTEMATICS
            if(exactLeps!=0){
                if(best_hyp.size()!=exactLeps){continue;}
                if(exactJets!=4){
                    if(good_jets.size()!=exactJets){continue;}
                }else{
                    if(good_jets.size()<exactJets){continue;}
                }
                if(exactBs!=2){
                    if(good_bjets.size()!=exactBs){continue;}
                }else{
                    if(good_bjets.size()<exactBs){continue;}
                }
            }




            //calculate control region weights for background estimates
            //we calculate the weight
            float crWeight = 0;
            //start with fake rate weight
            if (best_hyp_type==3||best_hyp_type>=6){
                float fakeWeight = 1;
                for (auto lep : best_hyp){
                    if (lep.idlevel()==SS::IDLevel::IDtight) continue;
                    float fakeRateValue = fakeRate(year, lep.id(), lep.conecorrpt(), lep.eta(), isData);
                    fakeWeight = fakeWeight * (fakeRateValue/(1-fakeRateValue));
                }
                crWeight = weight * fakeWeight;
            }

            //now, we move to the flip weight
            if (best_hyp_type==5){
                //cout << "******************" << endl;
                float flipWeight = 0;
                for (auto lep : tight_leptons){
                    //cout << lep.id() << " " << lep.mcid() << " " << lep.pt() << " " << lep.eta() << " " << lep.isFlip() << endl;
                    if (lep.absid()==13){continue;}
                    float flipRateValue = flipRate(nt.year(),lep.pt(),lep.eta());
                    if(nt.year()==2016){flipRateValue *= 1.1128218361769364;}
                    if(nt.year()==2017){flipRateValue *= 1.4447739510004296;}
                    if(nt.year()==2018){flipRateValue *= 1.3061565063535605;}
                    //std::cout << "flip rate: " << flipRateValue << endl;
                    flipWeight = flipWeight + (flipRateValue/(1-flipRateValue));
                    //std::cout << "flip Weight: " << flipWeight << endl;
                }
                crWeight = weight * flipWeight;
                //std::cout << "event crWeight: " << crWeight << endl;
            }
            if (debugPrints){std::cout << "passed crWeight for event " << nt.event() << endl;}
            if (debugPrints){std::cout << "elapsed time since start: " << duration_cast<seconds>(high_resolution_clock::now() - start).count() << endl;}
            if (debugPrints){std::cout << "elapsed time since b SF start: " << duration_cast<seconds>(high_resolution_clock::now() - startBOpening).count() << endl;}
            
            bool doVariations = 0;
            std::map<std::string, float> variationalWeights;

            if(!isData && doVarArg){
                doVariations = 1;
                float weightLepSFup = weight;
                float weightLepSFdown = weight;
                float weightTrigup = weight;
                float weightTrigdown = weight;

                for ( auto lep : best_hyp ) {weightLepSFup =    (weightLepSFup/leptonScaleFactor(nt.year(), lep.id(), lep.pt(), lep.eta(), ht))*
                                                                (leptonScaleFactor(nt.year(), lep.id(), lep.pt(), lep.eta(), ht)+
                                                                leptonScaleFactorError(nt.year(), lep.id(), lep.pt(), lep.eta(), ht));}
                for ( auto lep : best_hyp ) {weightLepSFdown =  (weightLepSFdown/leptonScaleFactor(nt.year(), lep.id(), lep.pt(), lep.eta(), ht))*
                                                                (leptonScaleFactor(nt.year(), lep.id(), lep.pt(), lep.eta(), ht)-
                                                                leptonScaleFactorError(nt.year(), lep.id(), lep.pt(), lep.eta(), ht));}


                if(!iterativeBTag){
                    float weightBtagSFup = weight;
                    float weightBtagSFdown = weight;
                    if (nt.year() == 2016){ weightBtagSFup = (weightBtagSFup/getBSF(nt.year(),good_jets,good_bjets,eff2016,deepjet_medium_reader_2016))*
                                                            getBSF(nt.year(),good_jets,good_bjets,eff2016,up_reader_2016);
                                            weightBtagSFdown = (weightBtagSFdown/getBSF(nt.year(),good_jets,good_bjets,eff2016,deepjet_medium_reader_2016))*
                                                            getBSF(nt.year(),good_jets,good_bjets,eff2016,down_reader_2016);}
                    
                    else if (nt.year() == 2017){weightBtagSFup = (weightBtagSFup/getBSF(nt.year(),good_jets,good_bjets,eff2017,deepjet_medium_reader_2017))*
                                                            getBSF(nt.year(),good_jets,good_bjets,eff2017,up_reader_2017);
                                                weightBtagSFdown = (weightBtagSFdown/getBSF(nt.year(),good_jets,good_bjets,eff2017,deepjet_medium_reader_2017))*
                                                            getBSF(nt.year(),good_jets,good_bjets,eff2017,down_reader_2017);}
                    
                    else if (nt.year() == 2018){weightBtagSFup = (weightBtagSFup/getBSF(nt.year(),good_jets,good_bjets,eff2018,deepjet_medium_reader_2018))*
                                                            getBSF(nt.year(),good_jets,good_bjets,eff2018,up_reader_2018);
                                                weightBtagSFdown = (weightBtagSFdown/getBSF(nt.year(),good_jets,good_bjets,eff2018,deepjet_medium_reader_2018))*
                                                            getBSF(nt.year(),good_jets,good_bjets,eff2018,down_reader_2018);}
                    variationalWeights["bTag_up"]=      weightBtagSFup;
                    variationalWeights["bTag_down"]=    weightBtagSFdown;
                }else{

                    vector<string> bVariations = {"lf","hf","hfstats1","hfstats2","lfstats1","lfstats2","cferr1","cferr2"};

                    for (uint b = 0; b < bVariations.size(); b++){
                        float weightBtagSFup = weight;
                        float weightBtagSFdown = weight;
                        if (nt.year() == 2016){ weightBtagSFup = (weightBtagSFup/getIterativeBSF(nt.year(),good_jets,good_bjets,deepjet_medium_reader_2016,"central"))*
                                                                getIterativeBSF(nt.year(),good_jets,good_bjets,deepjet_medium_reader_2016,"up_"+bVariations[b]);
                                                weightBtagSFdown = (weightBtagSFdown/getIterativeBSF(nt.year(),good_jets,good_bjets,deepjet_medium_reader_2016,"central"))*
                                                                getIterativeBSF(nt.year(),good_jets,good_bjets,deepjet_medium_reader_2016,"down_"+bVariations[b]);}
                        
                        else if (nt.year() == 2017){weightBtagSFup = (weightBtagSFup/getIterativeBSF(nt.year(),good_jets,good_bjets,deepjet_medium_reader_2017,"central"))*
                                                                getIterativeBSF(nt.year(),good_jets,good_bjets,deepjet_medium_reader_2017,"up_"+bVariations[b]);
                                                    weightBtagSFdown = (weightBtagSFdown/getIterativeBSF(nt.year(),good_jets,good_bjets,deepjet_medium_reader_2017,"central"))*
                                                                getIterativeBSF(nt.year(),good_jets,good_bjets,deepjet_medium_reader_2017,"down_"+bVariations[b]);}
                        
                        else if (nt.year() == 2018){weightBtagSFup = (weightBtagSFup/getIterativeBSF(nt.year(),good_jets,good_bjets,deepjet_medium_reader_2018,"central"))*
                                                                getIterativeBSF(nt.year(),good_jets,good_bjets,deepjet_medium_reader_2018,"up_"+bVariations[b]);
                                                    weightBtagSFdown = (weightBtagSFdown/getIterativeBSF(nt.year(),good_jets,good_bjets,deepjet_medium_reader_2018,"central"))*
                                                                getIterativeBSF(nt.year(),good_jets,good_bjets,deepjet_medium_reader_2018,"down_"+bVariations[b]);}

                        variationalWeights[bVariations[b]+"_up"] = weightBtagSFup;
                        variationalWeights[bVariations[b]+"_down"] = weightBtagSFdown;
                    }
                }
                                    

                weightTrigup = weightTrigup/triggerScaleFactor(f_trigEff, nt.year(), best_hyp[0].id(), best_hyp[1].id(), best_hyp[0].pt(), best_hyp[1].pt(), best_hyp[0].eta(), best_hyp[1].eta(), ht, 0);
                weightTrigup = weightTrigup*triggerScaleFactor(f_trigEff, nt.year(), best_hyp[0].id(), best_hyp[1].id(), best_hyp[0].pt(), best_hyp[1].pt(), best_hyp[0].eta(), best_hyp[1].eta(), ht, 1);

                weightTrigdown = weightTrigdown/triggerScaleFactor(f_trigEff, nt.year(), best_hyp[0].id(), best_hyp[1].id(), best_hyp[0].pt(), best_hyp[1].pt(), best_hyp[0].eta(), best_hyp[1].eta(), ht, 0);
                weightTrigdown = weightTrigdown*triggerScaleFactor(f_trigEff, nt.year(), best_hyp[0].id(), best_hyp[1].id(), best_hyp[0].pt(), best_hyp[1].pt(), best_hyp[0].eta(), best_hyp[1].eta(), ht, -1);
                
                variationalWeights["PU_up"]=        (weight/nt.puWeight())*nt.puWeightUp();
                variationalWeights["PU_down"]=      (weight/nt.puWeight())*nt.puWeightDown();
                variationalWeights["LepSF_up"]=     weightLepSFup;
                variationalWeights["LepSF_down"]=   weightLepSFdown;
                variationalWeights["Trigger_up"]=   weightTrigup;
                variationalWeights["Trigger_down"]= weightTrigdown;
                variationalWeights["btag_central"]= weight;
                auto pdfWeights = nt.LHEPdfWeight();
                auto scaleWeights = nt.LHEScaleWeight();
                // cout << nt.event() << endl;
                // cout << pdfWeights.size() << "   " << scaleWeights.size() << endl;
                // if(pdfWeights.size()>0){
                bool doPdfs = 0;
                bool doScale = 0;
                if(nt.year()==2016){
                    doPdfs=1;
                    doScale=1;
                }
                // for(uint i=0; i<pdfWeights.size(); i++){
                //     // cout << i << "   " << pdfWeights[i] << endl;
                //     if(pdfWeights[i]>2){doPdfs=0;}
                // }
                if (doPdfs){
                    for(uint i=0; i<pdfWeights.size(); i++){
                        // cout << i << "   " << pdfWeights[i] << endl;
                        variationalWeights[std::to_string(i)+"_"+std::to_string(i)+"_pdf_scale"]= weight*pdfWeights[i];
                        // cout << "pdf_scale_"+std::to_string(i) << "  " << variationalWeights["pdf_scale_"+std::to_string(i)] << "  " << weight << endl;
                    }
                }
                // }
                //for 2016 need indices 0,1,3,5,7,8
                //for 2017/2018 need indices 0,5,15,24,34,39
                vector<int> scaleIdx(6);
                if(nt.year()==2016){scaleIdx = {4,0,1,3,5,7,8};}
                else{scaleIdx = {0,5,15,24,34,39};}
                // if(scaleWeights.size()>0){
                if (doScale){
                    for(uint j=0; j<scaleIdx.size(); j++){
                        // cout << j << "   " << scaleIdx[j] << "   " << scaleWeights[scaleIdx[j]] << endl;
                        variationalWeights[std::to_string(j)+"_renorm_scale"]= weight*scaleWeights[scaleIdx[j]];
                    }
                }
                // }
                // if(isnan(weightTrigup)){cout << weight<< " "<< triggerScaleFactor(f_trigEff, nt.year(), best_hyp[0].id(), best_hyp[1].id(), best_hyp[0].pt(), best_hyp[1].pt(), best_hyp[0].eta(), best_hyp[1].eta(), ht, 0) << " " << triggerScaleFactor(f_trigEff, nt.year(), best_hyp[0].id(), best_hyp[1].id(), best_hyp[0].pt(), best_hyp[1].pt(), best_hyp[0].eta(), best_hyp[1].eta(), ht, 1) << endl;}
            }
            bool fill_BDT_MC = (((category==1) && (chainTitle=="fakes_mc")) ||
                                ((category==2) && (chainTitle=="flips_mc")) ||
                                ((category==3) && (chainTitle=="rares"   )) ||
                                ((chainTitle=="signal_tch"))                ||
                                ((chainTitle=="signal_tuh"))                );
            bool is_MC_fakes_flips = (chainTitle=="fakes_mc") || (chainTitle=="flips_mc");
            bool is_sig_rares = (chainTitle=="rares"   ) || (chainTitle=="signal_tch") || (chainTitle=="signal_tuh"); 
            bool fill_BDT_data_driven = (abs(crWeight) > 0.0);
            if (fill_BDT_MC) {
                if (((best_hyp_type == 4) || (((best_hyp.size() > 2) && (best_hyp_type==2)))) && make_BDT_MC_babies) {
                    std::map<std::string, Float_t> BDT_params = hct_booster.calculate_features(good_jets, good_bjets, best_hyp);
                    bdt_MC_baby.set_features(BDT_params, weight, variationalWeights);
                    //if (is_MC_fakes_flips) {//split fakes/flips into 2 datasets for training/testing
                    //    if (nt.event() % 2 == 0) {
                    //        bdt_MC_baby_training.set_features(BDT_params, weight*2, variationalWeights);
                    //    }
                    //    else if (nt.event() % 2 == 1) {
                    //        bdt_MC_baby_testing.set_features(BDT_params, weight*2, variationalWeights);
                    //    }
                    //}
                    //if (is_sig_rares) { //split signal and rares into 3 sections, training, testing, and evaluation of limits
                    //    if (nt.event() % 3 == 0) {
                    //        bdt_MC_baby_training.set_features(BDT_params, weight*3, variationalWeights);
                    //    }
                    //    else if (nt.event() % 3 == 1) {
                    //        bdt_MC_baby_testing.set_features(BDT_params, weight*3, variationalWeights);
                    //    }
                    //}
                }
            }
            else if (fill_BDT_data_driven) {
                if ((best_hyp_type == 5) && (make_BDT_flips_babies)){
                    Float_t BDT_crWeight = crWeight;
                    std::map<std::string, Float_t> BDT_params = hct_booster.calculate_features(good_jets, good_bjets, best_hyp);
                    bdt_flips_baby.set_features(BDT_params, BDT_crWeight, variationalWeights);
            
                } else if (make_BDT_fakes_babies && ((best_hyp_type==3) || (best_hyp_type>=6))) { 
                    Float_t BDT_crWeight = crWeight;
                    std::map<std::string, Float_t> BDT_params = hct_booster.calculate_features(good_jets, good_bjets, best_hyp);
                    bdt_fakes_baby.set_features(BDT_params, BDT_crWeight, variationalWeights);
                }
            }
            
        std::map<std::string, Float_t> HCT_BDT_params = hct_booster.calculate_features(good_jets, good_bjets, best_hyp);
        std::map<std::string, Float_t> HUT_BDT_params = hut_booster.calculate_features(good_jets, good_bjets, best_hyp);
        hct_booster.set_features(HCT_BDT_params);
        hut_booster.set_features(HUT_BDT_params);
        float HCT_BDT_score = float(hct_booster.get_score());
        float HUT_BDT_score = float(hut_booster.get_score());
        if (debugPrints) {
            std::cout << "HCT BDT score: " << HCT_BDT_score << std::endl;
            std::cout << "HUT BDT score: " << HUT_BDT_score << std::endl;
        }
            // if we've reached here we've passed the baseline selection
            // cout << category << endl;
            // fill histograms
            /*if ((chainTitle == "top"||chainTitle =="dy") && (best_hyp_type==4) ){
                hists.fill(chainTitleCh,best_hyp_type,best_hyp,good_jets,good_bjets,nt.MET_pt(),isVR_SR_fake,isVR_CR_fake,isVR_SR_flip,isVR_CR_flip,isEE,isEM,isME,isMM,isEFake,isMFake,isEE_flip,isEM_flip,weight,crWeight,doVariations,variationalWeights);
            }*/

            if (category == 1 && !(isSignal||isData)){
                // hists.fill("fakes_mc",best_hyp_type,best_hyp,good_jets,good_bjets,nt.MET_pt(),isVR_SR_fake,isVR_CR_fake,isVR_SR_flip,isVR_CR_flip,isEE,isEM,isME,isMM,isEFake,isMFake,isEE_flip,isEM_flip,weight,crWeight,doVariations,variationalWeights);
                hists.fill(chainTitleCh,best_hyp_type,best_hyp,good_jets,good_bjets,nt.MET_pt(),isVR_SR_fake,isVR_CR_fake,isVR_SR_flip,isVR_CR_flip,isEE,isEM,isME,isMM,isEFake,isMFake,isEE_flip,isEM_flip,weight,crWeight,doVariations,variationalWeights, HCT_BDT_score, HUT_BDT_score);
                if(print_sample_file){nFakes+=weight;}
            }else if (category == 2 && !(isSignal||isData)){
                hists.fill("flips_mc",best_hyp_type,best_hyp,good_jets,good_bjets,nt.MET_pt(),isVR_SR_fake,isVR_CR_fake,isVR_SR_flip,isVR_CR_flip,isEE,isEM,isME,isMM,isEFake,isMFake,isEE_flip,isEM_flip,weight,crWeight,doVariations,variationalWeights, HCT_BDT_score, HUT_BDT_score);
                // hists.fill(chainTitleCh,best_hyp_type,best_hyp,good_jets,good_bjets,nt.MET_pt(),isVR_SR_fake,isVR_CR_fake,isVR_SR_flip,isVR_CR_flip,isEE,isEM,isME,isMM,isEFake,isMFake,isEE_flip,isEM_flip,weight,crWeight,doVariations,variationalWeights);
                if(print_sample_file){nFlips+=weight;}
            }else if (category == 3 && !(isSignal||isData)){
                hists.fill("rares",best_hyp_type,best_hyp,good_jets,good_bjets,nt.MET_pt(),isVR_SR_fake,isVR_CR_fake,isVR_SR_flip,isVR_CR_flip,isEE,isEM,isME,isMM,isEFake,isMFake,isEE_flip,isEM_flip,weight,crWeight,doVariations,variationalWeights, HCT_BDT_score, HUT_BDT_score);
                // hists.fill(chainTitleCh,best_hyp_type,best_hyp,good_jets,good_bjets,nt.MET_pt(),isVR_SR_fake,isVR_CR_fake,isVR_SR_flip,isVR_CR_flip,isEE,isEM,isME,isMM,isEFake,isMFake,isEE_flip,isEM_flip,weight,crWeight,doVariations,variationalWeights);
                if(print_sample_file){nRares+=weight;}
            }else if (category == 4 && chainTitle!="fakes_mc" && chainTitle!="flips_mc" && chainTitle!="rares" && !(isSignal||isData)){
                hists.fill(chainTitleCh,best_hyp_type,best_hyp,good_jets,good_bjets,nt.MET_pt(),isVR_SR_fake,isVR_CR_fake,isVR_SR_flip,isVR_CR_flip,isEE,isEM,isME,isMM,isEFake,isMFake,isEE_flip,isEM_flip,weight,crWeight,doVariations,variationalWeights, HCT_BDT_score, HUT_BDT_score);
            }else if (isSignal||isData){
                hists.fill(chainTitleCh,best_hyp_type,best_hyp,good_jets,good_bjets,nt.MET_pt(),isVR_SR_fake,isVR_CR_fake,isVR_SR_flip,isVR_CR_flip,isEE,isEM,isME,isMM,isEFake,isMFake,isEE_flip,isEM_flip,/*weight*/1.0,crWeight,doVariations,variationalWeights, HCT_BDT_score, HUT_BDT_score);
            }
            //cout << "**********" << endl;
        }//loop over events
        if(print_sample_file){print_sample_yields(sample_file, filename ,nFlips, nFakes, nRares);}
        std::cout << "closing file " << file->GetName() << endl;
        if (debugPrints){std::cout << "elapsed time since start: " << duration_cast<seconds>(high_resolution_clock::now() - start).count() << endl;}
        if (debugPrints){std::cout << "elapsed time since b SF start: " << duration_cast<seconds>(high_resolution_clock::now() - startBOpening).count() << endl;}
        file->Close();
        delete file;
    } // loop over files
    // cout << "nLooseFlips: " << nLooseFlips << endl;
    // cout << "nTightFlips: " << nTightFlips << endl;
    // cout << "SS: " << endl;
    // cout << "    nDilep: " << nDilep_ss << endl;
    // cout << "    n1lep: " << n1lep_ss << endl;
    // cout << "    nDilep fakes: " << nDilep_fake_ss << endl;
    // cout << "    nDilep flips: " << nDilep_flip_ss << endl;
    // cout << "TL: " << endl;
    // cout << "    nDilep: " << nDilep_sf << endl;
    // cout << "    n1lep: " << n1lep_sf << endl;
    // cout << "    nDilep fakes: " << nDilep_fake_sf << endl;
    // cout << "    nDilep flips: " << nDilep_flip_sf << endl;

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start);
    if (debug_file.is_open()) debug_file.close();

    if (!quiet) cout << "processed " << nEventsTotal << " events in " << duration.count() << " seconds!!" << endl;
    cout << "processed " << nEventsTotal << " events in " << duration.count() << " seconds!!" << endl;
    //write histograms
    auto outFile = new TFile(outFileName.Data(), "recreate");
    if (!quiet) std::cout << "Writing " << chainTitleCh << " histograms to " << outFile->GetName() << std::endl;
    outFile->cd();
    hists.write();
    outFile->Close();
    if (make_BDT_fakes_babies){
        bdt_fakes_baby.close();
    }
    if (make_BDT_flips_babies){
        bdt_flips_baby.close();
    }
    if (make_BDT_MC_babies){
        bdt_MC_baby.close();
        //bdt_MC_baby_training.close();
        //bdt_MC_baby_testing.close();
    }
}
