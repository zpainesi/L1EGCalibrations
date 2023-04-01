import numpy as np
import os
import json, argparse
import ROOT as rt

from singleEGUtil import *


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f",
        "--inputFileList",
        help="Input File List",
        default="data/grid_eraG_default_nTTExtrap_newOpt.fls",
    )
    parser.add_argument(
        "-o",
        "--dest",
        help="destination To Use",
        default="results/SingleEG_nTTExtrap_jugad_check",
    )

    args = parser.parse_args()

    RateStudyNtuple = "Data , EpiZB , 3562616"
    EffStudyNtuple = "TagAndProbe Data Era G"
    DatsetInfoStr = "Rates studied using \n\t:\t " + RateStudyNtuple + "\n"
    DatsetInfoStr += "Efficiency studied using \n\t:\t " + EffStudyNtuple + "\n"

    Run2RateFile = "HistogramFile_Rate_Run3DataEraG_Run362616_calo_v6_default.root"
    Run2TurnOnFile = "HistogramFile_Eff_Run3DataEraG_calo_v6_default.root"
    Run2RateFile = "HistogramFile_Rate_Run3DataEraG_Run362616_calo_v6_default.root"
    Run2TurnOnFile = "HistogramFile_Eff_Run3DataEraG_calo_v6_default.root"
    OtherInfo = "\tBaseline Turnons : " + str(os.path.abspath(Run2TurnOnFile))
    OtherInfo += "\n\tBaseline Rates : " + str(os.path.abspath(Run2TurnOnFile))

    effMin = 0.85
    effMax = 0.95

    prefixBase = args.dest + "/SingleEG/FixedRateScans"
    condorFileList = args.inputFileList

    fileRate = rt.TFile(Run2RateFile)
    Run2RateHist = fileRate.Get("SingleEG_rate_TightIso")

    fileRun2TurnOns = rt.TFile(Run2TurnOnFile)

    Run2TurnOnMap = {}
    Run2AcceptanceMap = {}

    PassOptions = {}
    FailOptions = {}

    for E in range(10, 40):
        Dir = "TightIsoEt" + str(E)
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
        if "TightIso" not in hName:
            continue
            # print(hName)
        eT = int(hName.split("_")[3].replace("TightIsoEt", ""))
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

    fList = open(condorFileList)
    txt = fList.readlines()
    fileNames = []
    for l in txt:
        fileNames.append(l[:-1])
    fList.close()

    data = {}
    TFILE_STORE = {}
    count = 0
    for fName in fileNames:
        count += 1
        print(count, " / ", len(fileNames), " Adding File ", fName)
        file = rt.TFile(fName)
        RateForTurnons = file.Get("RateForTurnons")
        EtForTurnons = file.Get("EtForTurnons")
        AcceptanceForTurnons = file.Get("AcceptanceForTurnons")
        TurnOnProgressionDir = file.Get("turn_on_progression")
        RateProgressionDir = file.Get("rate_progression")
        data = updateDataStore(
            data,
            RateForTurnons,
            AcceptanceForTurnons,
            TurnOnProgressionDir,
            RateProgressionDir,
        )
        TFILE_STORE[fName] = file
    #     if count>2:
    #         break
    print()
    print("Loaded ", count, " Files")
    print("Loaded ", len(data), " options")

    opt = 354
    opt2 = 888

    eT = 30

    cvs = rt.TCanvas("cvs")
    cvs.cd()
    print(opt, opt2, eT)
    print(data.keys())
    print(data[opt].keys())
    print(data[opt]["FixedEtMetrics"].keys())
    data[opt]["FixedEtMetrics"][eT]["turnOn"].Draw()
    data[opt2]["FixedEtMetrics"][eT]["turnOn"].Draw("same")
    cvs.Draw()

    # isGoodTurnOn

    t1 = data[opt]["FixedEtMetrics"][eT]["turnOn"]
    t2 = data[opt2]["FixedEtMetrics"][eT]["turnOn"]

    y = rt.isGoodTurnON(t2, t1, 10, 19)
    print(y)

    ### Integral

    print(data[opt]["FixedEtMetrics"][eT])
    turnOn = data[opt]["FixedEtMetrics"][eT]["turnOn"]

    area = rt.getIntegral(turnOn, eT, eT - 4.0, eT + 4.0, "pref", rt.nullptr, False)
    print("area = ", area)

    fixedRates = [15, 16, 16.7, 17, 18]
    fixedRates = [16.7]

    dataForFixedRates = {}

    LT = data[931]["rateProgression"]
    effMin = 0.85
    effMax = 0.95
    OtherInfo += (
        "Efficeiecy range in which isGoodTurnOn is ran : [ "
        + str(effMin)
        + ","
        + str(effMax)
        + " ]\n"
    )
    for fR in fixedRates:
        print("Doing for Fixed Rate : ", fR)
        #     Baseline=getBaselineForFixedRate(fR,Run2RateHist,Run2TurnOnMap,Run2AcceptanceMap,effMin,effMax)
        Baseline = getBaselineForFixedRate(
            fR,
            Run2RateHist,
            Run2TurnOnMap,
            Run2AcceptanceMap,
            effMin,
            effMax,
            doFixedEt=False,
            fixedEt=28.0,
        )
        dataForFixedRates[fR] = {}
        dataForFixedRates[fR]["Baseline"] = Baseline
        dataForFixedRates[fR]["scan"] = getAllOptionsDataForFixedRate(
            data, fR, Baseline, LT
        )

    for tag in dataForFixedRates:
        print(tag, dataForFixedRates[tag]["Baseline"])

    # metricStr='area'
    metricStr = "acceptance"
    acceptanceTol = 0.85
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
    summaryFile.write(DatsetInfoStr)
    summaryFile.write(OtherInfo)

    for fRate in dataForFixedRates:
        summaryFile.write("=" * 20 + "\n")
        summaryFile.write("Summary at Fixed Rate = " + str(fRate) + " kHz : \n \n")

        prefix = prefixBase + "/" + str(fRate) + "kHz"
        print("For Fixed rate ", fRate, " Saving results to ", prefix)

        if not os.path.exists(prefix):
            print("making : ", prefix)
            os.system("mkdir -p " + prefix)

        #################   Baseline Details to the summary ###########

        summaryFile.write("  Baseline details the fixed rate \n")
        BaselineTurnOn = dataForFixedRates[fRate]["Baseline"][0]
        #     if BaselineTurnOn==ROOT.nullptr:
        #         print(fRate)
        #         h=0
        eTBaseLine = dataForFixedRates[fRate]["Baseline"][3]
        print(eTBaseLine)
        areaBaseline = rt.getIntegral(
            BaselineTurnOn,
            eTBaseLine,
            eTBaseLine - dETForAreaL,
            eTBaseLine + dETForAreaR,
            prefix + "/Baseline_",
            rt.nullptr,
            True,
        )
        acceptanceBaseline = dataForFixedRates[fRate]["Baseline"][4]
        baselineString = "\tBaseline eT   : " + str(eTBaseLine)
        baselineString += "\tBaseline Area : " + str(areaBaseline)
        baselineString += "\tBaseline Acceptance   : " + str(acceptanceBaseline)
        baselineString += "\n" + OtherInfo + "\n"
        summaryFile.write(baselineString)
        summaryFile.write("\n\n")
        #####################################################################

        fixedRateData = dataForFixedRates[fRate]["scan"]
        metric = np.asarray(fixedRateData[metricStr])
        sortedIdx = np.argsort(-1 * metric)
        isGoodTurnONMask = np.asarray(fixedRateData["isGoodTurnOn"]) == True

        allOptResults = open(prefix + "/allOptionsScan.txt", "w")
        allOptResults.write(DatsetInfoStr + "\n")
        allOptResults.write(OtherInfo + "\n")
        allOptResults.write("Full options scans at fixed Rate " + str(fRate) + " kHz\n")
        allOptResults.write("In decending order of the metric " + metricStr + "\n")
        allOptResults.write(baselineString + "\n")

        allOptResultsWithSelection = open(prefix + "/selectedOptions.txt", "w")
        allOptResultsWithSelection.write(
            "Full options scans at fixed Rate " + str(fRate) + " kHz\n"
        )
        allOptResultsWithSelection.write(
            "In decending order of the metric " + metricStr + "\n"
        )
        allOptResultsWithSelection.write(baselineString + "\n")

        sectedOptionsFile = rt.TFile(prefix + "/selectedOptions.root", "RECREATE")
        allOptionsFile = rt.TFile(prefix + "/allOptions.root", "RECREATE")
        baselineDir = allOptionsFile.mkdir("baselineDir")
        baselineDir.cd()

        keysToPrint = [
            "option",
            "option_pars",
            "eT",
            "rate",
            "isGoodTurnOn",
            "isTighterThanLT",
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

        eTsInTheScan = np.unique(np.asarray(fixedRateData["eT"])[isGoodTurnONMask])
        # eTsInTheScan=np.unique(np.asarray(fixedRateData['eT']))
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
            rateP = rt.TH1F(fixedRateData["rateProgression"][i])
            rateP.Write()
            allTurnOns = cDir.mkdir("allTurnons")
            allTurnOns.cd()
            #         for eT in data[opt]['FixedEtMetrics']:
            #             data[opt]['FixedEtMetrics'][eT]['turnOn'].Write()

            #         if not isGoodTurnONMask[i]:
            #             continue;
            if fixedRateData["eT"][i] > eTBaseLine + 2:
                continue
            if fixedRateData["area"][i] < acceptanceTol * areaBaseline:
                continue
            if fixedRateData["acceptance"][i] < acceptanceTol * acceptanceBaseline:
                continue

            count += 1
            eT = fixedRateData["eT"][i]
            if eT not in topOptionIdxForEachEt:
                topOptionIdxForEachEt[eT] = []
            if (fixedRateData["option_pars"][i][0] > 10.0) & (
                fixedRateData["option_pars"][i][1] < 0.3
            ):
                # & (fixedRateData['isTighterThanLT'][i] == True)): #NEW
                topOptionIdxForEachEt[eT].append(i)
            strToW = ""
            for j in range(len(keysToPrint)):
                ky = keysToPrint[j]
                strToW += str(fixedRateData[ky][i])
                strToW += "\t"
            strToW += "\n"
            allOptResultsWithSelection.write(strToW)
        allOptResults.close()
        print("Closing the selection file : ", allOptResultsWithSelection.name)
        allOptResultsWithSelection.close()
        allOptionsFile.Close()

        print("  Number of options passing the section ", count)
        summaryFile.write(
            "\n  Number of options passing the selection : " + str(count) + "\n\n"
        )

        #################   TOP n Options   to Summary ################

        summaryFile.write("  Best options for the fixed rate \n")
        baselineDir = sectedOptionsFile.mkdir("baselineDir")
        baselineDir.cd()
        BaselineTurnOn.Write()
        Run2RateHist.Write()

        top5opt_dir = sectedOptionsFile.mkdir("top_five_opions")
        top5opt_dir.cd()
        count = 0
        strToW = "\t"
        for j in range(len(keysToPrint)):
            strToW += keysToPrint[j] + "\t"
        summaryFile.write(strToW + "\n")
        turnOn = {}
        rateP = {}

        for i in sortedIdx:
            if fixedRateData["option_pars"][i][0] <= 10.0:
                continue
                # NEW
            #            if(fixedRateData['acceptance'][i] < 0.98 * maximum1): continue   #NEW
            if fixedRateData["eT"][i] > eTBaseLine + 2:
                #             print(fixedRateData['eT'][i] ," > ", eTBaseLine+2)
                continue
            if fixedRateData["area"][i] < acceptanceTol * areaBaseline:
                #             print( "area :",fixedRateData['area'][i] ,"  < ",1.0*areaBaseline)
                continue
            if fixedRateData["acceptance"][i] < acceptanceTol * acceptanceBaseline:
                #             print( "area :",fixedRateData['acceptance'][i]  ,"  < ",acceptanceTol*acceptanceBaseline)
                continue

            #         if not isGoodTurnONMask[i]:
            #             continue;
            count += 1

            eT = fixedRateData["eT"][i]
            area = rt.getIntegral(
                fixedRateData["turnOn"][i],
                eT,
                eT - dETForAreaL,
                eT + dETForAreaR,
                prefix + "/",
                BaselineTurnOn,
                True,
            )
            if count > 200:
                break
            strToW = "\t"
            for j in range(len(keysToPrint)):
                strToW += str(fixedRateData[keysToPrint[j]][i]) + "\t"
            summaryFile.write(strToW + "\n")

            cDir = top5opt_dir.mkdir(str(count))
            cDir.cd()
            turnOn[count] = fixedRateData["turnOn"][i]
            name = "top_" + str(count) + "_" + turnOn[count].GetName()
            turnOn[count].Write()

            rateP[count] = rt.TH1F(fixedRateData["rateProgression"][i])
            rateP[count].Write()

            allTurnOns = cDir.mkdir("allTurnons")
            allTurnOns.cd()
            opt = fixedRateData["option"][i]
            for eT in data[opt]["FixedEtMetrics"]:
                data[opt]["FixedEtMetrics"][eT]["turnOn"].Write()

        summaryFile.write("\n")

        #     continue

        #################   TOP n Options for unique eT  to Summary ################
        sectedOptionsFile.cd("/")
        top3OptInPt_dir = sectedOptionsFile.mkdir("top_three_options_in_Et")
        turnOn = {}
        rateP = {}
        # summaryFile.write("Hi theer")
        for eT in eTsInTheScan:
            summaryFile.write("  Best options for eT : " + str(eT) + "\n")
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
            # print(maximum)

            for i in topOptionIdxForEachEt[eT]:
                if fixedRateData["area"][i] < acceptanceTol * areaBaseline:
                    continue
                if fixedRateData["acceptance"][i] < acceptanceTol * acceptanceBaseline:
                    continue
                if fixedRateData["acceptance"][i] < 0.99 * maximum:
                    continue  # NEW
                #             if not isGoodTurnONMask[i]:
                #                 continue;
                count += 1
                area = rt.getIntegral(
                    fixedRateData["turnOn"][i],
                    eT,
                    eT - dETForAreaL,
                    eT + dETForAreaR,
                    prefix + "/",
                    BaselineTurnOn,
                    True,
                )
                if count > 5:
                    break
                strToW = "\t"
                for j in range(len(keysToPrint)):
                    strToW += str(fixedRateData[keysToPrint[j]][i]) + "\t"
                summaryFile.write(strToW + "\n")
                cDir = forCurrentPt_dir.mkdir(str(count))
                cDir.cd()
                turnOn[count] = fixedRateData["turnOn"][i]
                turnOn[count].Write()

                rateP[count] = fixedRateData["rateProgression"][i]
                rateP[count].Write()

                allTurnOns = cDir.mkdir("allTurnons")
                allTurnOns.cd()
                opt = fixedRateData["option"][i]
                for e in data[opt]["FixedEtMetrics"]:
                    data[opt]["FixedEtMetrics"][e]["turnOn"].Write()
            summaryFile.write("\n")
        sectedOptionsFile.Close()
    summaryFile.close()

    print("\nsummary written into ", summaryFileName)


if __name__ == "__main__":
    main()
