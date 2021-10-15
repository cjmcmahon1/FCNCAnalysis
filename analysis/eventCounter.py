import glob
import tqdm
import os

#VERSION = "fcnc_v10_SRonly_11june2021"
VERSION = "fcnc_v14_SRonly_30august2021"

BABY_DIR = "/hadoop/cms/store/user/ksalyer/FCNC_NanoSkim/{}".format(VERSION)
#BABY_DIR = "/hadoop/cms/store/user/dspitzba/FCNC_NanoSkim/{}".format(VERSION)
N_EVENTS_DIR = "./n_events"

if not os.path.exists(N_EVENTS_DIR): os.makedirs(N_EVENTS_DIR)

sample_dirs = glob.glob(BABY_DIR+"/*")

for sample_dir in tqdm.tqdm(sample_dirs):
    print(sample_dir)
    n_events_files = glob.glob(sample_dir+"/output_*_nevents.txt")
    n_total_events = 0
    n_effective_events = 0
    for n_events_file in n_events_files:
        counted_total = False
        counted_pos = False
        counted_neg = False
        with open(n_events_file, "r") as f_in:
            lines = f_in.readlines()
            for line in lines:
                line = line.split("\n")[0]
                if unicode(line, "utf-8").isnumeric():
                    if not counted_total:
                        n_total_events += int(line)
                        counted_total = True
                    elif not counted_pos:
                        n_effective_events += int(line)
                        counted_pos = True
                    elif not counted_neg:
                        n_effective_events -= int(line)
                        counted_neg = True
                counted_effective = counted_pos and counted_neg
                if counted_total and counted_effective:
                    break
            if not counted_total or not counted_effective:
                print(lines)
    sample = sample_dir.split("/")[-1]
    # if sample == "ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1_RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1_NANOAODSIM_fcnc_v10_SRonly_11june2021":
    #     print("found sample")
    with open(N_EVENTS_DIR+"/"+sample+"_n_events.txt", "w") as f_out:
        f_out.write(str(n_total_events)+"\n")
        f_out.write(str(n_effective_events)+"\n")
