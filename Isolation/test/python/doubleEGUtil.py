import numpy as np
import os
import ROOT as rt


rt.gInterpreter.ProcessLine('#include "Util.h"')


def getTheObjectsFromFile(aFile):
    histStore = {}
    for key in aFile.GetListOfKeys():
        kk = key.GetName()
        curObj = aFile.Get(kk)
        if curObj.Class_Name() == "TDirectoryFile":
            histStore[kk] = getTheObjectsFromFile(curObj)
        else:
            histStore[kk] = curObj
    return histStore


dETForAreaL = 7.0
dETForAreaR = 7.0


# In[18]:


def getOptionDescFromRateProgressionName(rateProgressionName):
    #     print(rateProgressionName)
    items = rateProgressionName.replace("v2_", "").split("_")
    desc = {}
    desc["option"] = int(items[3].replace("option", ""))
    desc["option_str"] = "_".join(items[3:7])
    desc["optionParams"] = (
        float(items[4].replace("p", ".")),
        float(items[5].replace("p", ".")),
        float(items[6].replace("p", ".")),
    )
    desc["EtThreshold"] = float(items[7].replace("p", ".").replace("Et", ""))

    return desc


def getThresholdRate(rate):
    return float(int((rate + 0.1) * 10) / 10.0)


# ## Definie the Main data storege

# #### Load the Et and Rate metrics

# In[10]:


def updateDataStore(
    data,
    RateForTurnons,
    AcceptanceForTurnons,
    TurnOnProgressionDir,
    RateProgressionDir,
    SingleEGRateProgressionDir,
):
    xAxisForNames = AcceptanceForTurnons.GetXaxis()
    nBins = AcceptanceForTurnons.GetNbinsX()
    AcceptanceMap = {}
    for i in range(nBins):
        option_histName = xAxisForNames.GetBinLabel(i)
        if option_histName == "":
            continue
        AcceptanceMap[option_histName] = AcceptanceForTurnons.GetBinContent(i)

    TurnOnMap = {}
    for ky in TurnOnProgressionDir.GetListOfKeys():
        hName = ky.GetName()
        TurnOnMap[hName] = TurnOnProgressionDir.Get(hName)

    RateMap = {}
    for ky in RateProgressionDir.GetListOfKeys():
        hName = ky.GetName()
        opt = (
            hName.replace("rate_Progression_double", "")
            .replace("v2_", "")
            .split("_")[0]
        )
        if opt == "":
            opt = int(hName.replace("rate_Progression_doublev2_", "").split("_")[1])
        opt = int(opt)
        RateMap[opt] = RateProgressionDir.Get(hName)

    singleEGRateMap = {}
    for ky in SingleEGRateProgressionDir.GetListOfKeys():
        hName = ky.GetName()
        opt = hName.replace("rate_Progression", "").replace("v2_", "").split("_")[0]
        if opt == "":
            opt = int(hName.replace("rate_Progression_", "").split("_")[1])
        opt = int(opt)
        singleEGRateMap[opt] = SingleEGRateProgressionDir.Get(hName)

    xAxisForNames = RateForTurnons.GetXaxis()
    nBins = RateForTurnons.GetNbinsX()
    for i in range(nBins):
        option_histName = xAxisForNames.GetBinLabel(i)
        if option_histName == "":
            continue

        rate = RateForTurnons.GetBinContent(i)
        optionDesc = getOptionDescFromRateProgressionName(option_histName)

        opt = optionDesc["option"]
        eT = optionDesc["EtThreshold"]
        rateProgression = RateMap[opt]

        if opt not in data:
            data[opt] = {}
            data[opt]["params"] = optionDesc["optionParams"]
            data[opt]["option_str"] = optionDesc["option_str"]
            data[opt]["rateProgression"] = rateProgression
            data[opt]["singleEGRateProgression"] = singleEGRateMap[opt]

        acceptance = AcceptanceMap[option_histName]
        turnOn = TurnOnMap[option_histName]
        area = 9.0
        area = rt.getIntegral(turnOn, eT, eT - 7.0, eT + 7.0, "", rt.nullptr, False)
        isGoodTurnOn = True

        if "FixedEtMetrics" not in data[opt]:
            data[opt]["FixedEtMetrics"] = {}
        if eT not in data[opt]["FixedEtMetrics"]:
            data[opt]["FixedEtMetrics"][eT] = {}
        data[opt]["FixedEtMetrics"][eT]["rate"] = rate
        data[opt]["FixedEtMetrics"][eT]["acceptance"] = acceptance
        data[opt]["FixedEtMetrics"][eT]["turnOn"] = turnOn
        data[opt]["FixedEtMetrics"][eT]["area"] = area
        data[opt]["FixedEtMetrics"][eT]["isGoodTurnOn"] = isGoodTurnOn

        if "FixedRateMetrics" not in data[opt]:
            data[opt]["FixedRateMetrics"] = {}
        rate_rounded = getThresholdRate(rate)

        if rate_rounded in data[opt]["FixedRateMetrics"]:
            if data[opt]["FixedRateMetrics"][rate_rounded]["eT"] < eT:
                pass
        else:
            data[opt]["FixedRateMetrics"][rate_rounded] = {}
            data[opt]["FixedRateMetrics"][rate_rounded]["eT"] = eT
            data[opt]["FixedRateMetrics"][rate_rounded]["acceptance"] = acceptance
            data[opt]["FixedRateMetrics"][rate_rounded]["_rate"] = rate
            data[opt]["FixedRateMetrics"][rate_rounded]["turnOn"] = turnOn
            data[opt]["FixedRateMetrics"][rate_rounded]["area"] = area
            data[opt]["FixedRateMetrics"][rate_rounded]["isGoodTurnOn"] = isGoodTurnOn

    return data


def getBaselineForFixedRate(
    fRate,
    Run2RateHist,
    Run2TurnOnMap,
    Run2AcceptanceMap,
    effMin=0.43,
    effMax=0.82,
    doFixedEt=False,
    fixedEt=12,
):
    nBins = Run2RateHist.GetNbinsX()
    eT = -1
    idx = -1
    rtFixHere = -1e6
    for i in range(nBins):
        rate = Run2RateHist.GetBinContent(i)
        if rate > fRate:
            continue
        else:
            if (fRate - rate) < (fRate - rtFixHere):
                rtFixHere = rate
                idx = i
    if idx >= 0:
        eT = int(Run2RateHist.GetBinCenter(idx)) + 10
    if doFixedEt:
        eT = fixedEt

    if eT < 0:
        print(" No thresold found with 'rate < Fixed Rate ' (", fRate, "kHz in run 2 ")
        return (rt.nullptr, -1, -1, -1)
    if eT not in Run2TurnOnMap:
        print(" No Turnons found for thresold eT = ", eT, " in run2 turn on collection")
        return (rt.nullptr, -1, -1, -1)

    turnon = Run2TurnOnMap[eT]
    nBins = turnon.GetN()
    idxA = -10
    idxB = -10
    for i in range(nBins):
        eff = turnon.GetPointY(i)
        if eff > effMin and idxA < 0:
            idxA = i
        if eff > effMax and idxB < 0:
            idxB = i
    if idxB < 0 or idxA < 0:
        print(
            "Something wrong with the idx finding for fRate ",
            fRate,
            " which gave eT ",
            eT,
        )
    print("Et for fixed rate ", fRate, " -> ", int(Run2RateHist.GetBinCenter(idx)))
    print(
        "\t Baseline found to be at eT = ",
        eT,
        " for fixed Rate ",
        fRate,
        " indexes : ",
        idxA,
        "(",
        "{0:0.3f}".format(turnon.GetPointY(idxA)),
        "/ ",
        effMin,
        ")",
        idxB,
        "(",
        "{0:0.3f}".format(turnon.GetPointY(idxB)),
        " / ",
        effMax,
        ")",
    )
    if eT in Run2AcceptanceMap:
        accep = Run2AcceptanceMap[eT]
    else:
        accep = -1.0
    return (Run2TurnOnMap[eT], idxA, idxB, eT, accep)


dataForValidation = {"baseLine": [], "p1": [], "p2": [], "rslt": [], "turnOn": []}


def getAllOptionsDataForFixedRate(
    data, F_Rate, Baseline=None, doFixedEt=False, fixedEt=12
):
    # for key in dataForValidation:
    #    dataForValidation[key].clear()
    dataForFixedRate = {
        "option": [],
        "option_pars": [],
        "eT": [],
        "eTLowEtLeg": [],
        "fEt": [],
        "rate": [],
        "area": [],
        "acceptance": [],
        "turnOnLowEtLeg": [],
        "turnOn": [],
        "isGoodTurnOn": [],
        "rateProgression": [],
        "singleEGRateProgression": [],
    }
    for opt in data:
        idx = -1
        rateFixHere = -1e6
        eTmaxSearched = -1.0
        rateForETmaxSearched = -1.0
        for rate in data[opt]["FixedRateMetrics"]:
            if eTmaxSearched < data[opt]["FixedRateMetrics"][rate]["eT"]:
                eTmaxSearched = data[opt]["FixedRateMetrics"][rate]["eT"]
                rateForETmaxSearched = rate
            if rate > F_Rate:
                continue
            else:
                if (F_Rate - rate) < (F_Rate - rateFixHere):
                    rateFixHere = rate
        if rateFixHere < 0.0:
            print(
                "\tNo Et thresolds could be reached for option ",
                opt,
                " at fixed Rate",
                F_Rate,
                "( rate : ",
                rateForETmaxSearched,
                " for eT ",
                eTmaxSearched,
                " )",
            )
            continue
        vals = data[opt]["FixedRateMetrics"][rateFixHere]
        eT = vals["eT"] + 10
        eTx = eT
        if doFixedEt:
            eT = fixedEt
        if eT > 40:
            continue
        isoLegTurnOn = data[opt]["FixedEtMetrics"][eT]["turnOn"]
        area = rt.getIntegral(
            isoLegTurnOn, eT, eT - 7.0, eT + 7.0, "", rt.nullptr, False
        )
        isGoodTurnON = False  #

        if Baseline[0]:
            isGoodTurnON = rt.isGoodTurnON(
                Baseline[0], isoLegTurnOn, Baseline[1], Baseline[2]
            )

        dataForFixedRate["option"].append(opt)
        dataForFixedRate["eT"].append(eTx)
        dataForFixedRate["fEt"].append(fixedEt)
        dataForFixedRate["eTLowEtLeg"].append(vals["eT"])
        dataForFixedRate["rate"].append(vals["_rate"])
        dataForFixedRate["area"].append(area)
        dataForFixedRate["isGoodTurnOn"].append(isGoodTurnON)
        dataForFixedRate["acceptance"].append(vals["acceptance"])
        dataForFixedRate["turnOn"].append(isoLegTurnOn)
        dataForFixedRate["turnOnLowEtLeg"].append(vals["turnOn"])
        dataForFixedRate["rateProgression"].append(data[opt]["rateProgression"])
        dataForFixedRate["singleEGRateProgression"].append(
            data[opt]["singleEGRateProgression"]
        )
        dataForFixedRate["option_pars"].append(data[opt]["params"])

    return dataForFixedRate
