#code to write datacard from tableMaker.py output
#run with python writeCards.py
#Kaitlin Salyer, April 2021

#import some useful packages
import ROOT
import time
import numpy as np
import sys
import pandas as pd
import os

#hardcoded variables other users should customize
outdir = "/home/users/ksalyer/FranksFCNC/ana/analysis/datacards/"
indir = "test"
outdir = outdir+indir+"/"
if not os.path.exists(outdir): os.makedirs(outdir)
# years = [2016, 2017, 2018]
years = [2016, 2017, 2018]
includeSignalInObs = False


fakeCRStatYld_2016 = [  1507,1342,104,
                        757,1162,155,
                        439,847,222,
                        282,376,5,
                        167,210,15,
                        46,80,15,
                        22,48,8
                        ]
fakeCRStatYld_2017 = [  2304,2282,222,
                        1191,1908,278,
                        629,1302,416,
                        328,508,5,
                        171,300,30,
                        62,137,32,
                        38,50,28
                        ]
fakeCRStatYld_2018 = [  2919,3112,297,
                        1431,2452,425,
                        766,1774,612,
                        503,766,13,
                        211,467,54,
                        93,167,49,
                        35,89,35
                        ]


fakeCRStatYldDict_2016 = {  "2_2_0":1507,"2_3_0":757,"2_4_0":439,
                            "2_2_1":1342,"2_3_1":1162,"2_4_1":847,
                            "2_2_2":104,"2_3_2":155,"2_4_2":222,
                            "3_1_0":282,"3_2_0":167,"3_3_0":46,"3_4_0":22,
                            "3_1_1":376,"3_2_1":210,"3_3_1":80,"3_4_1":48,
                            "3_1_2":5,"3_2_2":15,"3_3_2":15,"3_4_2":8
                            }
fakeCRStatYldDict_2017 = {  "2_2_0":2304,"2_3_0":1191,"2_4_0":629,
                            "2_2_1":2282,"2_3_1":1908,"2_4_1":1302,
                            "2_2_2":222,"2_3_2":278,"2_4_2":416,
                            "3_1_0":328,"3_2_0":171,"3_3_0":62,"3_4_0":38,
                            "3_1_1":508,"3_2_1":300,"3_3_1":137,"3_4_1":50,
                            "3_1_2":5,"3_2_2":30,"3_3_2":32,"3_4_2":28
                            }
fakeCRStatYldDict_2018 = {  "2_2_0":2919,"2_3_0":1431,"2_4_0":766,
                            "2_2_1":3112,"2_3_1":2452,"2_4_1":1774,
                            "2_2_2":297,"2_3_2":425,"2_4_2":612,
                            "3_1_0":503,"3_2_0":211,"3_3_0":93,"3_4_0":35,
                            "3_1_1":766,"3_2_1":467,"3_3_1":167,"3_4_1":89,
                            "3_1_2":13,"3_2_2":54,"3_3_2":49,"3_4_2":35
                            }


fakeCRStatErr_2016 = {  "2_2_0":3.50365,"2_3_0":4.99353,"2_4_0":6.52687,
                        "2_2_1":3.60495,"2_3_1":3.88372,"2_4_1":4.56676,
                        "2_2_2":13.1561,"2_3_2":10.5987,"2_4_2":8.83165,
                        "3_1_0":8.31807,"3_2_0":10.6145,"3_3_0":17.0339,"3_4_0":29.254,
                        "3_1_1":7.11159,"3_2_1":9.03126,"3_3_1":16.6831,"3_4_1":19.3909,
                        "3_1_2":53.7054,"3_2_2":33.4056,"3_3_2":31.4495,"3_4_2":66.0103
                        }
fakeCRStatErr_2017 = {  "2_2_0":3.33519,"2_3_0":4.41385,"2_4_0":6.08993,
                        "2_2_1":3.11701,"2_3_1":3.30036,"2_4_1":4.05022,
                        "2_2_2":10.4076,"2_3_2":9.20855,"2_4_2":7.40525,
                        "3_1_0":8.19731,"3_2_0":11.5854,"3_3_0":19.1633,"3_4_0":23.6831,
                        "3_1_1":6.19351,"3_2_1":8.0818,"3_3_1":12.7294,"3_4_1":25.1242,
                        "3_1_2":53.354,"3_2_2":26.9341,"3_3_2":25.3507,"3_4_2":26.0582,
                        }
fakeCRStatErr_2018 = {  "2_2_0":3.22747,"2_3_0":4.51919,"2_4_0":6.10372,
                        "2_2_1":2.9054,"2_3_1":3.22509,"2_4_1":3.87638,
                        "2_2_2":9.39522,"2_3_2":8.25866,"2_4_2":6.87119,
                        "3_1_0":7.3366,"3_2_0":11.4367,"3_3_0":17.4301,"3_4_0":27.1185,
                        "3_1_1":5.62841,"3_2_1":7.37149,"3_3_1":13.5088,"3_4_1":19.717,
                        "3_1_2":59.9883,"3_2_2":23.2358,"3_3_2":23.367,"3_4_2":26.6004
                        }


fakeSystErr_2016 = [    4.46067,4.76551,5.03131,
                        5.62097,4.84993,5.98169,
                        6.07675,5.59812,6.29573,
                        6.79787,6.49155,8.72948,
                        6.70704,5.32922,5.47438,
                        9.27242,8.69589,7.82636,
                        8.82801,6.4532,19.5941
                        ]
print(sum(fakeSystErr_2016)/len(fakeSystErr_2016))
fakeSystErr_2017 = [    6.63146,8.08039,7.74266,
                        7.24731,7.56606,7.26889,
                        9.42543,9.85621,8.45383,
                        6.98604,8.04912,40.6501,
                        9.96351,9.18872,7.21745,
                        10.7108,7.79917,10.4827,
                        8.9602,13.3502,10.228
                        ]
print(sum(fakeSystErr_2017)/len(fakeSystErr_2017))
fakeSystErr_2018 = [    7.71618,8.71154,7.85155,
                        7.72384,9.50327,7.64877,
                        7.31554,10.0837,8.39966,
                        7.96532,9.96087,16.6508,
                        10.14,10.1111,10.7022,
                        14.7398,10.4591,13.3193,
                        21.8245,13.9245,11.5897
                        ]
print(sum(fakeSystErr_2018)/len(fakeSystErr_2018))


if not os.path.exists(outdir): os.makedirs(outdir)

for y in years:
    #first, we load the txt output from the tableMaker.py script into a dataframe
    #we will manipulate these data, save it into a different dataframe, and print to an output file
    df = pd.read_csv("/home/users/ksalyer/FranksFCNC/ana/analysis/outputs/tables/"+indir+"/tableMaker_"+str(y)+".txt")
    fakeEst_df = pd.read_csv("/home/users/ksalyer/FranksFCNC/ana/analysis/outputs/tables/"+indir+"/fakeEstyields_"+str(y)+".txt")
    flipEst_df = pd.read_csv("/home/users/ksalyer/FranksFCNC/ana/analysis/outputs/tables/"+indir+"/flipEstyields_"+str(y)+".txt")
    raresSyst_df = pd.read_csv("/home/users/ksalyer/FranksFCNC/ana/analysis/outputs/variations/rares_"+str(y)+".txt")
    tchSyst_df = pd.read_csv("/home/users/ksalyer/FranksFCNC/ana/analysis/outputs/variations/signal_tch_"+str(y)+".txt")
    tuhSyst_df = pd.read_csv("/home/users/ksalyer/FranksFCNC/ana/analysis/outputs/variations/signal_tuh_"+str(y)+".txt")
    # print (fakeEst_df)
    # print (flipEst_df)
    # df = pd.read_csv("/home/users/ksalyer/FranksFCNC/ana/analysis/outputs/tables/aug09_withST/tableMaker_"+str(y)+".txt")
    # fakeEst_df = pd.read_csv("/home/users/ksalyer/FranksFCNC/ana/analysis/outputs/tables/aug09_withST/fakeEstyields_"+str(y)+".txt")
    # flipEst_df = pd.read_csv("/home/users/ksalyer/FranksFCNC/ana/analysis/outputs/tables/aug09_withST/flipEstyields_"+str(y)+".txt")
    # # print df
    # print fakeEst_df
    # print flipEst_df
    fakeCRStatYld = []
    fakeCRStatYldDict = {}
    fakeCRStatErr = {}
    fakeSystErr = []
    if y == 2016: 
        fakeCRStatYld = fakeCRStatYld_2016
        fakeCRStatYldDict = fakeCRStatYldDict_2016
        fakeCRStatErr = fakeCRStatErr_2016
        fakeSystErr = fakeSystErr_2016
        yr = "16"
    if y == 2017: 
        fakeCRStatYld = fakeCRStatYld_2017
        fakeCRStatYldDict = fakeCRStatYldDict_2017
        fakeCRStatErr = fakeCRStatErr_2017
        fakeSystErr = fakeSystErr_2017
        yr = "17"
    if y == 2018: 
        fakeCRStatYld = fakeCRStatYld_2018
        fakeCRStatYldDict = fakeCRStatYldDict_2018
        fakeCRStatErr = fakeCRStatErr_2018
        fakeSystErr = fakeSystErr_2018
        yr = "18"
    print("got yields and stat errors")

    #now we have imported the data and manipulated it into the categories we want
    #we will do the rest in a loop over signals
    signals = ["signal_tuh","signal_tch"]
    for s in signals:
        if "tch" in s: signalSyst_df = tchSyst_df
        else: signalSyst_df = tuhSyst_df
        outfileName = "datacard_"
        outfileName+=s
        outfileName+="_"
        outfileName+=str(y)
        outfileName+=".txt"

        #signal region bins
        nLeps = [2, 3]
        nJets = [2,3,4]
        nMLJets = [1,2,3,4]
        nBtags = [0,1,2]
        # nBtags = [1,2]
        numBins = (len(nMLJets)+len(nJets))*len(nBtags)
        nProc = ["signal", "rares", "fakes_mc", "flips_mc"]
        numBackgrounds = len(nProc)-1
        # systematicSources = ["LepSF","PU","Trigger","bTag","jes"]
        systematicSources = [   "LepSF","PU","Trigger","jes",
                                "lf","lfstats1","lfstats2",
                                "hf","hfstats1","hfstats2",
                                "cferr1","cferr2"]

        #make some headers for my dataframe columns
        dcColumns = []
        binNames = []
        binNamesObs = []
        procNames = []
        procIndex = []
        for l in nLeps:
            if l == 2:
                lep = "dilep"
                numJets = list(nJets)
            if l == 3:
                lep = "trilep"
                numJets = list(nMLJets)
            for j in numJets:
                for b in nBtags:
                    counter = 0
                    srName_base = lep+"_"+str(j)+"_"+str(b)
                    binName = srName_base
                    while len(binName) < 20:
                        binName+=" "
                    binNamesObs.append(binName)
                    for p in nProc:
                        srName = srName_base+"_"+p

                        binNames.append(binName)

                        while len(srName) < 20:
                            srName+=" "
                        dcColumns.append(srName)

                        while len(p) < 20:
                            p+=" "
                        procNames.append(p)

                        counterString = str(counter)
                        while len(counterString) < 20:
                            counterString+=" "
                        procIndex.append(counterString)
                        counter+=1

        # ok, now I have headers, I can start making the titles for my rows
        rowTitles = []
        numParameters = 0
        for p in nProc:
            iterator = 0
            while iterator < numBins:
                if p == "fakes_mc":
                    title = "fkStat"+yr+"_"+str(iterator)
                    numParameters+=1
                    yld = str(fakeCRStatYld[iterator])
                    while len(title) <(16-len(yld)):
                        title += " "
                    title += "gmN "+yld

                else:
                    if p=="flips_mc": title = "flp_stat_"+yr+"_"+str(iterator)
                    elif p=="signal": title = "sig_stat_"+yr+"_"+str(iterator)
                    else: title = p+"_stat_"+yr+"_"+str(iterator)
                    # title = p+"_stat_"+yr+"_"+str(iterator)
                    numParameters+=1
                    while len(title) <17:
                        title+=" "
                    title+="lnN"
                rowTitles.append(title)
                iterator+=1
            # print(rowTitles)

        for p in nProc:
            # if "rares" in p: continue
            # if "signal" in p: continue
            title = p+"_syst"
            while len(title)<17:
                title+=" "
            title+="lnN"
            rowTitles.append(title)
            numParameters+=1

        title = "fakeRate_syst"
        while len(title)<17:
            title+=" "
        title+="lnN"
        rowTitles.append(title)
        numParameters+=1

        for source in systematicSources:
            if source == "lf" or source == "hf" or "cferr" in source:
                title = source
            else:
                title = source+"_"+yr
            while len(title)<17:
                title+=" "
            title+="lnN"
            rowTitles.append(title)
            numParameters+=1


        #dataframe I will print to datacard file
        dcard_df = pd.DataFrame(index = rowTitles, columns = dcColumns)
        # print(dcColumns)
        print("defined output dataframe")
        # print(dcard_df)


        #ok, now I want to know the stat uncertainty as a percentage of the yield
        for p in nProc:
            if p == "signal":
                proc = s
            else:
                proc = p
            statUnc = []
            for l in nLeps:
                if l == 2:
                    numJets = list(nJets)
                if l == 3:
                    numJets = list(nMLJets)
                for j in numJets:
                    for b in nBtags:
                        # print (p, l, j, b)
                        #calculate signal percentage for statistical unc.
                        if proc == "fakes_mc":
                            srbin = str(l)+"_"+str(j)+"_"+str(b)
                            yld = fakeCRStatYldDict[srbin]
                            row = fakeEst_df.loc[ (df["nLeptons"]==l) & (df["nJets"]==j) & (df["nBtags"]==b) ]
                            err = row["data estimate"].values[0]
                            # print(yld, err, err/yld)
                        elif proc == "flips_mc":
                            if l == 3:
                                yld = 0.0
                                err = 0
                            else:
                                row = flipEst_df.loc[ (df["nLeptons"]==l) & (df["nJets"]==j) & (df["nBtags"]==b) ]
                                yld = row['data estimate'].values[0]
                                err = row["data estimate error"].values[0]
                        else:
                            # print df
                            # print df.loc[ (df["nLeptons"]==l) ]
                            # print df.loc[ (df["nJets"]==j) ]
                            # print df.loc[ (df["nBtags"]==b) ]
                            row = df.loc[ (df["nLeptons"]==l) & (df["nJets"]==j) & (df["nBtags"]==b) ]
                            # print row
                            yld = row[proc].values[0]
                            err = row[proc+" error"].values[0]
                        
                        if yld > 0:
                            dcPercentage = round(err/yld,3)
                        # else:
                        #     dcPercentage = 1

                        if dcPercentage>=1:
                            dcPercentage = 1

                        if proc == "fakes_mc": statUnc.append([l,j,b,dcPercentage])
                        else: statUnc.append([l,j,b,1+dcPercentage])
            # print(statUnc)
            for i in range(len(statUnc)):
                lep = statUnc[i][0]
                if lep == 2: l = "dilep"
                if lep == 3: l = "trilep"
                jet = statUnc[i][1]
                btag = statUnc[i][2]
                unc = statUnc[i][3]

                cTitle = str(l)+"_"+str(jet)+"_"+str(btag)+"_"+p
                while len(cTitle) < 20:
                    cTitle+=" "

                if p =="fakes_mc":
                    rTitle = "fkStat"+yr+"_"+str(i)
                    yld = str(fakeCRStatYld[i])
                    while len(rTitle) <(16-len(yld)):
                        rTitle += " "
                    rTitle += "gmN "+yld
                else:
                    if p=="flips_mc": rTitle = "flp_stat_"+yr+"_"+str(i)
                    elif p=="signal": rTitle = "sig_stat_"+yr+"_"+str(i)
                    else: rTitle = p+"_stat_"+yr+"_"+str(i)
                    # rTitle = p+"_stat_"+yr+"_"+str(i)
                    while len(rTitle) < 17:
                        rTitle+=" "
                    rTitle+="lnN"
                #print("*************")
                # if p=="flips_mc":
                #     print(cTitle)
                #     print(rTitle)
                #     print(unc)


                for column in dcard_df:
                    # print column
                    if column==cTitle:
                        filler = str(unc)
                        while len(filler)<20:
                            filler+=" "
                        dcard_df.at[rTitle,column] = filler
                        # if p=="flips_mc": print filler
                        # if p=="flips_mc": print dcard_df.at[rTitle,column]
                    else:
                        filler = "-"
                        while len(filler) < 20:
                            filler += " "
                        dcard_df.at[rTitle,column] = filler
                        # print filler
                # if p == "flips_mc": print dcard_df
        print("filled stat uncertainties")
        #print(dcard_df)

        #get MC yields in correct order of bins/processes
        rates = []
        observation = []
        for l in nLeps:
            if l == 2:
                numJets = list(nJets)
            if l == 3:
                numJets = list(nMLJets)
            # print(numJets)
            for j in numJets:
                for b in nBtags:
                    # print fakeEst_df
                    fakerow = fakeEst_df.loc[ (df["nLeptons"]==l) & (df["nJets"]==j) & (df["nBtags"]==b) ]
                    fliprow = flipEst_df.loc[ (df["nLeptons"]==l) & (df["nJets"]==j) & (df["nBtags"]==b) ]
                    row     = df.loc[ (df["nLeptons"]==l) & (df["nJets"]==j) & (df["nBtags"]==b) ]
                    obsYld = row["rares"].values[0] + fakerow['data estimate'].values[0] + fliprow['data estimate'].values[0]
                    if includeSignalInObs: obsYld += row[s].values[0]
                    obsYld = round(obsYld,0)
                    obsString = str(obsYld)
                    while len(obsString)<20:
                        obsString+=" "
                    observation.append(obsString)

                    for p in nProc:
                        if p == "signal":
                            p = s
                            yld = row[p].values[0]
                        elif p == "fakes_mc":
                            fakerow = fakeEst_df.loc[ (df["nLeptons"]==l) & (df["nJets"]==j) & (df["nBtags"]==b) ]
                            yld = fakerow["data estimate"].values[0]
                        elif p == "flips_mc":
                            fliprow = flipEst_df.loc[ (df["nLeptons"]==l) & (df["nJets"]==j) & (df["nBtags"]==b) ]
                            yld = fliprow["data estimate"].values[0]
                        else:
                            yld = row[p].values[0]

                        # if yld<0:
                        # if yld<=0:
                        #     yld = 0.01

                        yldString = str(yld)
                        while len(yldString)<20:
                            yldString+=" "
                        rates.append(yldString)



        #filling dummy systematic uncertainties
        for p in nProc:
            if "signal" in p:
                unc = "0.8/1.2"
                # continue
            if p == "flips_mc":
                # if y == 2016: unc = '1.1'
                # if y == 2017: unc = '1.4'
                # if y == 2018: unc = '1.3'
                unc = "0.7/1.3"
            elif p == "rares":
                unc = "0.7/1.3"
                # continue
            elif p == "fakes_mc":
                unc = "0.6/1.4"

            while len(unc)<20:
                unc+=" "
            
            rTitle = p+"_syst"
            while len(rTitle)<17:
                rTitle+=" "
            rTitle+="lnN"

            for column in dcard_df:
                if ((p in column) and not (p=="flips_mc" and "trilep" in column)):
                    filler = unc
                    while len(filler)<20:
                        filler+=" "
                    dcard_df.at[rTitle,column] = filler
                else:
                    filler = "-"
                    while len(filler) < 20:
                        filler += " "
                    dcard_df.at[rTitle,column] = filler

        fakeRateIter = 0
        for column in dcard_df:
            rTitle = "fakeRate_syst"
            while len(rTitle)<17:
                rTitle+=" "
            rTitle+="lnN"
            if "fakes_mc" in column:
                filler = str(1+fakeSystErr[fakeRateIter]/100)
                while len(filler)<20:
                    filler+=" "
                dcard_df.at[rTitle,column] = filler
                fakeRateIter+=1
            else:
                filler = "-"
                while len(filler)<20:
                    filler+=" "
                dcard_df.at[rTitle,column]=filler

        for source in systematicSources:
            if source == "lf" or source == "hf" or "cferr" in source:
                rTitle = source
            else:
                rTitle = source+"_"+yr
            # rTitle = source+"_"+yr
            while len(rTitle)<17:
                rTitle+=" "
            rTitle+="lnN"
            for column in dcard_df:
                if ("rares" in column) or ("signal" in column):
                    if "dilep" in column: l=2
                    else: l=3
                    srbin = column.split("_")
                    j = int(srbin[1])
                    b = int(srbin[2])
                    if "rares" in column:
                        row = raresSyst_df.loc[(raresSyst_df["nLeptons"]==l) & (raresSyst_df["nJets"]==j) & (raresSyst_df["nBtags"]==b)]
                        colTitleUp = "rares_"+source+"_up"
                        colTitleDown = "rares_"+source+"_down"
                        filler = str(round(row[colTitleDown],3))+"/"+str(round(row[colTitleUp],3))
                        while len(filler)<20:
                            filler+=" "
                        dcard_df.at[rTitle,column]=filler
                    elif "signal" in column:
                        row = signalSyst_df.loc[(signalSyst_df["nLeptons"]==l) & (signalSyst_df["nJets"]==j) & (signalSyst_df["nBtags"]==b)]
                        colTitleUp = s+"_"+source+"_up"
                        colTitleDown = s+"_"+source+"_down"
                        filler = str(round(row[colTitleDown],3))+"/"+str(round(row[colTitleUp],3))
                        while len(filler)<20:
                            filler+=" "
                        dcard_df.at[rTitle,column]=filler
                else:
                    filler = "-"
                    while len(filler)<20:
                        filler+=" "
                    dcard_df.at[rTitle,column]=filler

        print("filled syst uncertainties")


        #define output file and write to output file
        outfile = open(outdir+outfileName,"w")
        binHeadersObs = "bin                 \t"
        for b in binNamesObs:
            binHeadersObs+=b
            binHeadersObs+="\t"
        binHeadersObs+="\n"
        binHeaders = "bin                 \t"
        for b in binNames:
            binHeaders+=b
            binHeaders+="\t"
        binHeaders+="\n"
        procHeaders = "process             \t"
        for p in procNames:
            procHeaders+=p
            procHeaders+="\t"
        procHeaders+="\n"
        pInxHeaders = "process             \t"
        for i in procIndex:
            pInxHeaders+=i
            pInxHeaders+="\t"
        pInxHeaders+="\n"
        rateHeaders = "rate                \t"
        for r in rates:
            rateHeaders+=r
            rateHeaders+="\t"
        rateHeaders+="\n"
        obsHeaders = "observation         \t"
        for o in observation:
            obsHeaders+=o
            obsHeaders+="\t"
        obsHeaders+="\n"
        imaxHeader = "imax "+str(numBins)+" number of channels\n"
        jmaxHeader = "jmax "+str(numBackgrounds)+" number of backgrounds\n"
        kmaxHeader = "kmax "+str(numParameters)+" number of nuisance parameters\n"
        outfile.write(imaxHeader)
        outfile.write(jmaxHeader)
        outfile.write(kmaxHeader)
        outfile.write("shapes * * FAKE\n\n")
        outfile.write(binHeadersObs)
        outfile.write(obsHeaders)
        outfile.write("\n")
        outfile.write(binHeaders)
        outfile.write(procHeaders)
        outfile.write(pInxHeaders)
        outfile.write(rateHeaders)
        outfile.write(dcard_df.to_csv(sep="\t", index=True, header=False))
        outfile.close()

