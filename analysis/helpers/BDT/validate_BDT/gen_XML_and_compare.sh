#!/bin/bash
cp ~/public_html/BDT/HCT/random_search_optimized_all/booster_HCT.model ~/public_html/BDT/HUT/random_search_optimized_all/booster_HUT.model .
printf "Copied over booster models from public folder."
conda activate coffeadev
pushd /home/users/cmcmahon/CMSSW_10_2_9/src/tW_scattering
source activate_conda.sh
popd
printf "Generating test CSV files from XGB BDTs.\n"
python3 BDT_test.py
printf "Done.\nConverting XGB model to XML.\n"
python3 make_XML.py
printf "Done.\nSwitching to ROOT environment.\n"
conda deactivate
pushd .
source ~/scripts/setup_looper.sh
popd
printf "Done.\nEvaluating TMVA BDTs in ROOT.\n"
root -b -l start_loop.C+
