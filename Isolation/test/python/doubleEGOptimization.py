#!/usr/bin/env python

import numpy as np
import os, argparse
import ROOT as rt

import warnings

warnings.filterwarnings("ignore")


from doubleEGUtil import *


def main():
    # RateStudyNtuple='Data , EpiZB , 356X148'
    # EffStudyNtuple='TagAndProbe Data Era G'

    # DatsetInfoStr='Rates studied using \n\t:\t '+RateStudyNtuple+'\n'
    # DatsetInfoStr+='Efficiency studied using \n\t:\t '+EffStudyNtuple+'\n'

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f",
        "--inputFileList",
        help="Input File List",
        default="data/grid_eraG_default_nTTExtrap_newOpt.fls",
    )
    parser.add_argument(
        "-o", "--dest", help="destination To Use", default="results/DoublEG_default"
    )

    args = parser.parse_args()

    RateStudyNtuple = "Data , EpiZB , 356X148"
    EffStudyNtuple = "TagAndProbe Data Era G"
    DatsetInfoStr = "Rates studied using \n\t:\t " + RateStudyNtuple + "\n"
    DatsetInfoStr += "Efficiency studied using \n\t:\t " + EffStudyNtuple + "\n"

    # Run2RateFile='data/HistgramFile_Rate_Run3MC_Run3IsoLUT_PU44to52.root'
    # Run2RateFile='/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/sbaradia/LUTs_2023/Baselines/0p5/HistogramFile_Rate_Run3DataEraG_Run362616_calo_v6_ZS0p5.root'
    # Run2TurnOnFile='/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/sbaradia/LUTs_2023/Baselines/0p5/HistogramFile_Eff_Run3DataEraG_calo_v6_ZS0p5.root'
    # Run2RateFile='/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/sbaradia/LUTs_2023/Baselines/default/HistogramFile_Rate_Run3DataEraG_Run362616_calo_v6_default.root'
    # Run2RateFile='/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/sbaradia/LUTs_2023/Baselines/default/HistogramFile_Rate_Run3DataEraG_Run362616_calo_v6_default.root'
    Run2RateFile = "HistogramFile_Rate_Run3DataEraG_Run362616_calo_v6_default.root"
    Run2TurnOnFile = "HistogramFile_Eff_Run3DataEraG_calo_v6_default.root"
    OtherInfo = "\tBaseline Turnons : " + str(os.path.abspath(Run2TurnOnFile))
    OtherInfo += "\n\tBaseline Rates : " + str(os.path.abspath(Run2TurnOnFile))

    acceptanceTol = 0.85
    areaTol = 0.85

    effMin = 0.85
    effMax = 0.95

    # condorFileList='data/gridC_v2.fls'

    prefixBase = args.dest + "/DoubleEG/FixedRateScans"
    condorFileList = args.inputFileList

    DatsetInfoStr += "Input fileList :  " + condorFileList + "\n"

    fileRate = rt.TFile(Run2RateFile)
    Run2RateHist = fileRate.Get("DoubleEG_rate_LooseIso")

    fileRun2TurnOns = rt.TFile(Run2TurnOnFile)

    Run2TurnOnMap = {}
    Run2AcceptanceMap = {}
    PassOptions = {}
    FailOptions = {}

    for E in range(10, 40):
        Dir = "LooseIsoEt" + str(E)
        #    print(dire)
        PFDir = fileRun2TurnOns.Get(Dir)
        # fileRun2TurnOns.Get(dire)
        for ky in PFDir.GetListOfKeys():
            hName = ky.GetName()
            if "pT_pass" in hName:
                PassOptions[int(E)] = PFDir.Get(hName)
            elif "pT_fail" in hName:
                FailOptions[int(E)] = PFDir.Get(hName)
            else:
                continue

    GraphDir = fileRun2TurnOns.Get("TGraphs")
    for ky in GraphDir.GetListOfKeys():
        hName = ky.GetName()
        if "LooseIso" not in hName:
            continue
            # print(hName)
        eT = int(hName.split("_")[3].replace("LooseIsoEt", ""))
        # print(hName)
        Run2TurnOnMap[eT] = GraphDir.Get(hName)
        acceptance = -1.0
        if eT not in PassOptions or eT not in FailOptions:
            print("\t\t Problem !! eT not found in pass and fain maps !! eT = ", eT)
        else:
            acceptance = (
                PassOptions[eT].Integral()
                * 1.0
                / (PassOptions[eT].Integral() + FailOptions[eT].Integral())
            )
            Run2AcceptanceMap[eT] = acceptance
            #   print('For ',eT,' Adding ',hName, "Acceptance = ",acceptance)

    # #### Load the Files and data

    # In[19]:

    fList = open(condorFileList)
    txt = fList.readlines()
    fileNames = []
    for l in txt:
        fileNames.append(l[:-1])
    fList.close()

    # In[20]:

    data = {}
    TFILE_STORE = {}
    count = 0
    for fName in fileNames:
        count += 1
        print(count, " / ", len(fileNames), " Adding File ", fName)
        file = rt.TFile(fName)
        RateForTurnons = file.Get("RateDForTurnons")
        EtForTurnons = file.Get("EtForTurnons")
        AcceptanceForTurnons = file.Get("AcceptanceForTurnons")
        TurnOnProgressionDir = file.Get("turn_on_progression")
        RateProgressionDir = file.Get("rate_progression_double")
        SingleEGRateProgressionDir = file.Get("rate_progression")
        data = updateDataStore(
            data,
            RateForTurnons,
            AcceptanceForTurnons,
            TurnOnProgressionDir,
            RateProgressionDir,
            SingleEGRateProgressionDir,
        )
        TFILE_STORE[fName] = file
    #     if count>2:
    #         break
    print()
    print("Loaded ", count, " Files")
    print("Loaded ", len(data), " options")

    # ### Testing the root macro functions

    # In[21]:

    print(data[892]["rateProgression"].Integral())
    c = rt.TCanvas()
    data[892]["rateProgression"].Draw()
    Run2RateHist.Draw("same")
    c.Draw()

    # In[22]:

    opt = 101
    opt2 = 101
    cvs = rt.TCanvas("cvs")
    cvs.cd()
    data[opt]["FixedEtMetrics"][10.0]["turnOn"].Draw()
    data[opt2]["FixedEtMetrics"][10.0]["turnOn"].Draw("same")
    cvs.Draw()

    # isGoodTurnOn

    t1 = data[opt]["FixedEtMetrics"][10.0]["turnOn"]
    t2 = data[opt2]["FixedEtMetrics"][10.0]["turnOn"]

    y = rt.isGoodTurnON(t2, t1, 10, 19)
    print(y)

    ### Integral

    eT = 10
    print(data[opt]["FixedEtMetrics"][eT])
    turnOn = data[opt]["FixedEtMetrics"][eT]["turnOn"]

    area = rt.getIntegral(turnOn, eT, eT - 4.0, eT + 4.0, "pref", rt.nullptr, False)
    print("area = ", area)

    # ## Get Baseline and Best Option for Fixed Rates

    # In[23]:

    # ## Scanning the fixed rates

    # In[24]:

    fixedRates = [0.5, 1, 1.5, 2, 2.5, 3, 3.5, 3.9, 4.0, 4, 5, 5, 6]
    fixedRates = [3.9]

    # In[25]:

    dataForFixedRates = {}
    effMin = 0.85
    effMax = 0.95
    doFixed = False
    fEt = -1
    OtherInfo += (
        "Efficeiecy range in which isGoodTurnOn is ran : [ "
        + str(effMin)
        + ","
        + str(effMax)
        + " ]\n"
    )
    if doFixed:
        OtherInfo += "Doing at fixed Et threshold ! eT : " + str(fEt) + " \n"
    for fR in fixedRates:
        if not doFixed:
            print("Doing for Fixed Rate : ", fR)
        else:
            print("Doing for Fixed Et   : ", fEt)
        Baseline = getBaselineForFixedRate(
            fR,
            Run2RateHist,
            Run2TurnOnMap,
            Run2AcceptanceMap,
            effMin,
            effMax,
            doFixedEt=doFixed,
            fixedEt=fEt,
        )
        #     print(Baseline)
        #     continue
        dataForFixedRates[fR] = {}
        dataForFixedRates[fR]["Baseline"] = Baseline
        dataForFixedRates[fR]["scan"] = getAllOptionsDataForFixedRate(
            data, fR, Baseline, doFixedEt=doFixed, fixedEt=fEt
        )

    # ## Export the results

    # In[26]:

    # metricStr='area'
    metricStr = "acceptance"
    areaTol = 0.85
    accepTol = 0.85
    if not os.path.exists(prefixBase):
        print("making : ", prefixBase)
        os.system("mkdir -p " + prefixBase)
    else:
        print("removing : ", prefixBase)
        os.system("rm -r " + prefixBase)
        print("making : ", prefixBase)
        os.system("mkdir -p " + prefixBase)

    summaryFileName = prefixBase + "/" + "summaryForFixedRateScans.txt"
    summaryFile = open(summaryFileName, "w")
    summaryFile.write(DatsetInfoStr + "\n\n")
    summaryFile.write(OtherInfo + "\n\n")

    for fRate in dataForFixedRates:
        summaryFile.write("=" * 20 + "\n")
        summaryFile.write("Summary at Fixed Rate = " + str(fRate) + " : \n \n")

        prefix = prefixBase + "/" + str(fRate) + "kHz"
        print("For Fixed rate ", fRate, " Saving results to ", prefix)

        if not os.path.exists(prefix):
            print("making : ", prefix)
            os.system("mkdir -p " + prefix)

        #################   Baseline Details to the summary ###########

        summaryFile.write("  Baseline details the fixed rate \n")
        BaselineTurnOn = dataForFixedRates[fRate]["Baseline"][0]
        eTBaseLine = dataForFixedRates[fRate]["Baseline"][3]
        areaBaseline = rt.getIntegral(
            BaselineTurnOn,
            eTBaseLine,
            eTBaseLine - 7.0,
            eTBaseLine + 7.0,
            prefix + "/Baseline_",
            rt.nullptr,
            True,
        )
        acceptanceBaseline = dataForFixedRates[fRate]["Baseline"][4]
        baselineString = "\tBaseline eT   : " + str(eTBaseLine)
        baselineString += "\tBaseline Area : " + str(areaBaseline)
        baselineString += "\tBaseline Acceptance   : " + str(acceptanceBaseline)
        summaryFile.write(baselineString)
        print(baselineString)
        summaryFile.write("\n\n")

        #####################################################################

        fixedRateData = dataForFixedRates[fRate]["scan"]
        #        print(fixedRateData['option'])

        metric = np.asarray(fixedRateData[metricStr])
        sortedIdx = np.argsort(-1 * metric)

        isSaturatingEarly = np.asarray(fixedRateData["isGoodTurnOn"]) == True
        hasBetterAreaMask = np.asarray(fixedRateData["area"]) >= areaBaseline
        isFeasibleTurnONMask = (np.asarray(fixedRateData["eT"]) - eTBaseLine) < 2.2

        isGoodTurnONMask = np.logical_and(isSaturatingEarly, isFeasibleTurnONMask)
        isGoodTurnONMask = hasBetterAreaMask
        isGoodTurnONMask = np.logical_and(hasBetterAreaMask, isFeasibleTurnONMask)
        isGoodTurnONMask = isSaturatingEarly
        print(
            "\t\t isGoodTurnONMask  : ",
            sum(isGoodTurnONMask),
            " / ",
            len(isGoodTurnONMask),
        )

        allOptResults = open(prefix + "/allOptionsScan.txt", "w")
        allOptResults.write("Full options scans at fixed Rate " + str(fRate) + " kHz\n")
        allOptResults.write("In decending order of the metric " + metricStr + "\n")
        allOptResults.write(baselineString + "\n")
        allOptResults.write(OtherInfo + "\n")

        allOptResultsWithSelection = open(prefix + "/selectedOptions.txt", "w")
        allOptResultsWithSelection.write(
            "Full options scans at fixed Rate " + str(fRate) + " kHz\n"
        )
        allOptResultsWithSelection.write(
            "In decending order of the metric " + metricStr + "\n"
        )
        allOptResultsWithSelection.write(baselineString + "\n")
        sectedOptionsFile = rt.TFile(prefix + "/selectedOptions.root", "RECREATE")
        baselineDir = sectedOptionsFile.mkdir("baselineDir")
        baselineDir.cd()
        BaselineTurnOn.Write()
        Run2RateHist.Write()

        allOptionsFile = rt.TFile(prefix + "/allOptions.root", "RECREATE")
        baselineDir = allOptionsFile.mkdir("baselineDir")
        baselineDir.cd()
        BaselineTurnOn.Write()
        keysToPrint = [
            "option",
            "option_pars",
            "eTLowEtLeg",
            "eT",
            "rate",
            "isGoodTurnOn",
            "area",
            "acceptance",
        ]
        strToW = ""

        for j in range(len(keysToPrint)):
            ky = keysToPrint[j]
            strToW += ky
            strToW += "\t"
        strToW += "\n"
        allOptResults.write(strToW)
        allOptResultsWithSelection.write(strToW)

        # eTsInTheScan=np.unique(np.asarray(fixedRateData['eTLowEtLeg'])[isGoodTurnONMask])
        eTsInTheScan = np.unique(np.asarray(fixedRateData["eTLowEtLeg"]))
        topOptionIdxForEachEt = {i: [] for i in eTsInTheScan}
        count = 0
        for i in sortedIdx:
            strToW = ""
            for j in range(len(keysToPrint)):
                ky = keysToPrint[j]
                strToW += str(fixedRateData[ky][i])
                strToW += "\t"
            strToW += "\n"
            allOptResults.write(strToW)

            opt = fixedRateData["option"][i]
            cDir = allOptionsFile.mkdir(str(opt))
            cDir.cd()
            fixedRateData["turnOn"][i].Write()
            fixedRateData["rateProgression"][i].Write()
            fixedRateData["singleEGRateProgression"][i].Write()
            allTurnOns = cDir.mkdir("allTurnons")
            allTurnOns.cd()
            #         for eT in data[opt]['FixedEtMetrics']:
            #             if eT<22 or eT>28:
            #                 continue
            #             if data[opt]['FixedEtMetrics'][eT]['area']< 0.7*areaBaseline:
            #                 continue
            #             data[opt]['FixedEtMetrics'][eT]['turnOn'].Write()

            eT = fixedRateData["eT"][i]
            eTL = fixedRateData["eTLowEtLeg"][i]
            if (fixedRateData["option_pars"][i][0] > 10.0) & (
                fixedRateData["option_pars"][i][1] < 0.3
            ):
                topOptionIdxForEachEt[eTL].append(i)

            #         if not isGoodTurnONMask[i]:
            #             continue;

            #         if fixedRateData['eT'][i] > eTBaseLine+3:
            #             continue
            if fixedRateData["eTLowEtLeg"][i] > eTBaseLine + 3:
                continue
            if fixedRateData["area"][i] < areaTol * areaBaseline:
                continue
            if fixedRateData["acceptance"][i] < acceptanceTol * acceptanceBaseline:
                continue
            count += 1
            strToW = ""

            #         print(fixedRateData['eT'][i]," / ",eTBaseLine+3 ,
            #               fixedRateData['area'][i] ," / ",areaTol*areaBaseline,
            #               fixedRateData['acceptance'][i]," / ",acceptanceTol*acceptanceBaseline )

            for j in range(len(keysToPrint)):
                ky = keysToPrint[j]
                strToW += str(fixedRateData[ky][i])
                strToW += "\t"
            strToW += "\n"
            allOptResultsWithSelection.write(strToW)
        allOptResults.close()
        allOptResultsWithSelection.close()
        allOptionsFile.Close()
        print("  Number of options passing the section ", count)
        summaryFile.write(
            "\n  Number of options passing the selection : " + str(count) + "\n\n"
        )

        #################   TOP n Options   to Summary ################
        summaryFile.write("  Best options for the fixed rate \n")
        count = 0
        strToW = "\t"
        top5opt_dir = sectedOptionsFile.mkdir("top_five_opions")
        top5opt_dir.cd()
        for j in range(len(keysToPrint)):
            strToW += keysToPrint[j] + "\t"
        summaryFile.write(strToW + "\n")

        nMax = 10
        print(sortedIdx, len(sortedIdx))

        # maximum1 = float('-inf')
        # for i in sortedIdx:
        #   if(fixedRateData['acceptance'][i] > maximum1):
        #      maximum1=fixedRateData['acceptance'][i]

        for i in sortedIdx:
            #         if not isGoodTurnONMask[i]:
            #             continue;
            if fixedRateData["option_pars"][i][0] <= 10.0:
                continue
                # NEW
            # if(fixedRateData['acceptance'][i] < 0.98 * maximum1): continue                #NEW
            if fixedRateData["eTLowEtLeg"][i] > eTBaseLine + 3:
                continue
            if fixedRateData["area"][i] < areaTol * areaBaseline:
                continue
            if fixedRateData["acceptance"][i] < acceptanceTol * acceptanceBaseline:
                continue

            count += 1
            eT = fixedRateData["eT"][i]
            if doFixed:
                eT = fixedRateData["fEt"][i]

            print("Doing Et", eT)
            area = rt.getIntegral(
                fixedRateData["turnOn"][i],
                eT,
                eT - 7.0,
                eT + 7.0,
                prefix + "/top" + str(nMax) + "_",
                BaselineTurnOn,
                True,
            )
            if count > nMax:
                break
            strToW = "\t"
            for j in range(len(keysToPrint)):
                strToW += str(fixedRateData[keysToPrint[j]][i]) + "\t"
            summaryFile.write(strToW + "\n")

            cDir = top5opt_dir.mkdir(str(count))
            cDir.cd()
            fixedRateData["turnOn"][i].Write()
            fixedRateData["rateProgression"][i].Write()

            allTurnOns = cDir.mkdir("allTurnons")
            allTurnOns.cd()
            opt = fixedRateData["option"][i]
            for eT in data[opt]["FixedEtMetrics"]:
                data[opt]["FixedEtMetrics"][eT]["turnOn"].Write()

        summaryFile.write("\n")
        #     break
        #################   TOP n Options for unique eT  to Summary ################
        nMax = 5
        top3OptInPt_dir = sectedOptionsFile.mkdir("top_three_options_in_Et")
        for eT in eTsInTheScan:
            summaryFile.write(
                "  Best options for eT : [ " + str(eT) + "," + str(eT + 10) + " ]\n"
            )
            print("  Best options for eT : [ " + str(eT) + "," + str(eT + 10) + " ]\n")
            count = 0
            strToW = "\t"
            for j in range(len(keysToPrint)):
                strToW += keysToPrint[j] + "\t"
            summaryFile.write(strToW + "\n")
            forCurrentPt_dir = top3OptInPt_dir.mkdir("eT" + str(eT))
            forCurrentPt_dir.cd()

            maximum = float("-inf")  # NEW
            for i in topOptionIdxForEachEt[eT]:
                if fixedRateData["acceptance"][i] > maximum:
                    maximum = fixedRateData["acceptance"][i]
            print(maximum)

            for i in topOptionIdxForEachEt[eT]:
                if fixedRateData["acceptance"][i] < 0.99 * maximum:
                    continue  # NEW

                #             if not isGoodTurnONMask[i]:
                #                 continue;
                eTx = eT + 10
                if doFixed:
                    eTx = fixedRateData["fEt"][i]
                count += 1
                print(
                    "eTx ",
                    eTx,
                    " : ",
                    prefix + "/top" + str(nMax) + "inPtBin" + str(int(eT)) + "_",
                )
                area = rt.getIntegral(
                    fixedRateData["turnOn"][i],
                    eTx,
                    eTx - 7.0,
                    eTx + 7.0,
                    prefix + "/top" + str(nMax) + "inPtBin" + str(int(eT)) + "_",
                    BaselineTurnOn,
                    True,
                )
                if count > nMax:
                    break
                strToW = "\t"
                for j in range(len(keysToPrint)):
                    strToW += str(fixedRateData[keysToPrint[j]][i]) + "\t"
                print(
                    "Area : ",
                    fixedRateData["area"][i],
                    area,
                    " opt : ",
                    fixedRateData["option"][i],
                    " eTx : ",
                    eTx,
                    "",
                    fixedRateData["option"][i],
                )
                summaryFile.write(strToW + "\n")

                cDir = forCurrentPt_dir.mkdir(str(count))
                cDir.cd()
                fixedRateData["turnOn"][i].Write()
                fixedRateData["rateProgression"][i].Write()

                allTurnOns = cDir.mkdir("allTurnons")
                allTurnOns.cd()
                opt = fixedRateData["option"][i]
                for e in data[opt]["FixedEtMetrics"]:
                    data[opt]["FixedEtMetrics"][e]["turnOn"].Write()
        sectedOptionsFile.Close()

    summaryFile.close()
    print("\nsummary written into ", summaryFileName)


# In[ ]:

if __name__ == "__main__":
    main()
