import numpy as np
import os
import ROOT as rt


rt.gInterpreter.ProcessLine('#include "Util.h"')
dETForAreaL = 7.0
dETForAreaR = 7.0


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


def getOptionDescFromRateProgressionName(rateProgressionName):
    items = rateProgressionName.replace("_v2", "").split("_")
    desc = {}
    desc["option"] = int(items[4].replace("option", ""))
    desc["option_str"] = "_".join(items[3:8])
    desc["optionParams"] = (
        float(items[5].replace("p", ".")),
        float(items[6].replace("p", ".")),
        float(items[7].replace("p", ".")),
    )
    desc["EtThreshold"] = float(items[8].replace("p", ".").replace("Et", ""))
    # print(rateProgressionName)
    # print(items)
    # print(desc)
    # exit(0)
    return desc


def getThresholdRate(rate):
    return float(int((rate + 0.1) * 10) / 10.0)


def updateDataStore(
    data, RateForTurnons, AcceptanceForTurnons, TurnOnProgressionDir, RateProgressionDir
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
        opt = int(hName.replace("rate_Progression", "").split("_")[1])
        RateMap[opt] = RateProgressionDir.Get(hName)

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

        acceptance = AcceptanceMap[option_histName]
        #        print(option_histName)
        #        print(TurnOnMap.keys())
        #       for keys, value in TurnOnMap.items():
        #           if keys == option_histName:
        #               print(keys)
        turnOn = TurnOnMap[option_histName]
        area = 9.0
        area = rt.getIntegral(
            turnOn, eT, eT - dETForAreaL, eT + dETForAreaR, "", rt.nullptr, False
        )

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
        #        print(rate)
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
    fixedEt=28.0,
):
    nBins = Run2RateHist.GetNbinsX()
    eT = -1
    idx = -1
    rtFixHere = -1e6
    if doFixedEt:
        for i in range(nBins):
            eT = int(Run2RateHist.GetBinCenter(i))
            if abs(eT - fixedEt) < 0.1:
                idx = i
                break

    else:
        for i in range(nBins):
            rate = Run2RateHist.GetBinContent(i)
            if rate > fRate:
                continue
            else:
                if (fRate - rate) < (fRate - rtFixHere):
                    rtFixHere = rate
                    idx = i
    if idx >= 0:
        eT = int(Run2RateHist.GetBinCenter(idx))
    if eT < 0:
        print(
            " No thresold found with 'rate < Fixed Rate ' (",
            fRate,
            "kHz in run 2  | eT<0",
        )
        return (rt.nullptr, -1, -1, -1)
    if eT not in Run2TurnOnMap:
        print(
            " No Turnons found for thresold eT = ",
            eT,
            " in run2 turn on collection | eT not in Run2TurnOnMap",
        )
        return (rt.nullptr, -1, -1, -1)

    turnon = Run2TurnOnMap[eT]
    nBins = turnon.GetN()
    idxA = -10
    idxB = -10
    print("N in turnons ", nBins)
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
    print(
        "\t Baseline found to be at eT = ",
        eT,
        " for fixed Rate ",
        fRate,
        " indexes : ",
        idxA,
        "(",
        "{0:0.3f}".format(turnon.GetPointY(idxA)),
        ")",
        idxB,
        "(",
        "{0:0.3f}".format(turnon.GetPointY(idxB)),
        ")",
    )
    return (Run2TurnOnMap[eT], idxA, idxB, eT, Run2AcceptanceMap[eT])


# dataForValidation={'baseLine':[],'p1':[],'p2':[],'rslt':[],'turnOn':[]}
def getAllOptionsDataForFixedRate(data, F_Rate, Baseline=None, LT=None):
    dataForFixedRate = {
        "option": [],
        "option_pars": [],
        "eT": [],
        "rate": [],
        "area": [],
        "acceptance": [],
        "turnOn": [],
        "isGoodTurnOn": [],
        "isTighterThanLT": [],
        "rateProgression": [],
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

        isGoodTurnON = False  #
        if Baseline[0]:
            isGoodTurnON = rt.isGoodTurnON(
                Baseline[0], vals["turnOn"], Baseline[1], Baseline[2]
            )
        isTighterThanLT = True
        if LT != None:
            for i in range(LT.GetNbinsX()):
                if LT.GetBinLowEdge(i + 1) < 24.0:
                    continue
                if LT.GetBinLowEdge(i + 1) > 38.0:
                    break

                if data[opt]["rateProgression"].GetBinContent(
                    i + 1
                ) > 0.9 * LT.GetBinContent(i + 1):
                    isTighterThanLT = False
        if isTighterThanLT:
            print(" Got A tight opt on as opt ", opt)

        #             dataForValidation['baseLine'].append(Baseline[0])
        #             dataForValidation['p1'].append(Baseline[1])
        #             dataForValidation['p2'].append(Baseline[2])
        #             dataForValidation['turnOn'].append(vals['turnOn'])
        #             dataForValidation['rslt'].append(isGoodTurnON)
        #             print("Got isGoodTurnOns as ",isGoodTurnON)

        dataForFixedRate["option"].append(opt)
        dataForFixedRate["eT"].append(vals["eT"])
        dataForFixedRate["rate"].append(vals["_rate"])
        dataForFixedRate["area"].append(vals["area"])
        dataForFixedRate["isGoodTurnOn"].append(isGoodTurnON)
        dataForFixedRate["isTighterThanLT"].append(isTighterThanLT)
        dataForFixedRate["acceptance"].append(vals["acceptance"])
        dataForFixedRate["turnOn"].append(vals["turnOn"])
        dataForFixedRate["rateProgression"].append(data[opt]["rateProgression"])
        dataForFixedRate["option_pars"].append(data[opt]["params"])
    return dataForFixedRate
