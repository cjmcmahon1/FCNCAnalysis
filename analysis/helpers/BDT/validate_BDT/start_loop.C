{
    gROOT->ProcessLine(".L ../../../../../NanoTools/NanoCORE/NANO_CORE.so");
    gROOT->ProcessLine(".L compare_score.C+g");
    compare_score();
    return;
}
