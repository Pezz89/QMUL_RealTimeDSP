#!/usr/bin/env python

from __future__ import division, print_function
import numpy as np
import glob
import os
import pdb

correctOfflinePeakCount = 0
incorrectOfflinePeakCount = 0
correctOfflineLabelCount = 0
incorrectOfflineLabelCount = 0
totalOfflineLabelCount = 0
totalOfflinePeakCount = 0
totalOfflinePeakEst = 0

fs = 1000

def evaluate(offlineData, springerData, offlineFilepath, springerFilepath):

    global correctOfflinePeakCount
    global incorrectOfflinePeakCount
    global correctOfflineLabelCount
    global incorrectOfflineLabelCount
    global totalOfflineLabelCount
    global totalOfflinePeakCount
    global totalOfflinePeakEst
    # Check which peaks were correctly identified in time by finding a
    # corresponding peak in springer data within 100ms
    springerTimes = springerData[:, 0]
    offlineTimes = offlineData[:, 0]
    # Remove peaks from first second (which is used for calibration in real
    # time analysis)
    springerTimes = springerTimes[springerTimes > fs]
    offlineTimes = offlineTimes[offlineTimes > fs]

    # Mark all offline peak times that are within range of each peak from
    # the springer data
    c = np.isclose(springerTimes, offlineTimes[np.newaxis].T, 0, 0.1*fs)
    # Any peaks that are within 100ms are marked as correct
    results = np.any(c, axis = 0)

    # Quantify results for offline data
    # Correct peak count:
    cpc = np.count_nonzero(results)
    correctOfflinePeakCount += cpc
    # Incorrect peak count:
    ipc = np.size(results) - np.count_nonzero(results)
    incorrectOfflinePeakCount += ipc

    # Total number of peak estimates:
    tpe = offlineTimes.size
    totalOfflinePeakEst += tpe
    # Total number of peaks actually present in the signal:
    tpc = springerTimes.size
    totalOfflinePeakCount += tpc

    # For peaks that were correctly identified check which ones were
    # correctly labeled
    springerLabels = springerData[:, 1]
    offlineLabels = offlineData[:, 1]

    # Get indexes where labels are correct
    inds = np.where(c)
    # Get correct labels and estimated labels for these indicies
    sL = springerLabels[inds[1]]
    oL = offlineLabels[inds[0]]

    # Create boolean mask to determine correct vs. incorrect classification
    correctLabels = sL == oL

    # Quantify results
    clc = np.count_nonzero(correctLabels)
    correctOfflineLabelCount += clc
    ilc = np.size(correctLabels) - np.count_nonzero(correctLabels)
    incorrectOfflineLabelCount += ilc
    tlc = springerLabels.size
    totalOfflineLabelCount += tlc
    print("-------------------------------------------------------------------------------")
    print("Filepath: {0}".format(offlineFilepath))
    print("Correct Peak Count:\t\t\t{0}\t({1:.01f}%)".format(cpc, (cpc/tpc)*100))
    print("Incorrect Peak Count:\t\t\t{0}\t({1:.01f}%)".format(ipc, (ipc/tpc)*100))
    print("Correct Label Count:\t\t\t{0}\t({1:.01f}%)".format(clc, (clc/tlc)*100))
    print("Incorrect Label Count:\t\t\t{0}\t({1:.01f}%)".format(ilc, (ilc/tlc)*100))
    print("Total Number of Peaks Estimate:\t\t{0}\t({1:+.0f})".format(tpe, tpe-tpc))
    print("Total Number of Peaks Count:\t\t{0}".format(tpc))
    print("-------------------------------------------------------------------------------")



def main():
    global correctOfflinePeakCount
    global incorrectOfflinePeakCount
    global correctOfflineLabelCount
    global incorrectOfflineLabelCount
    global totalOfflineLabelCount
    global totalOfflinePeakCount
    global totalOfflinePeakEst

    offlineResults = glob.glob("../OfflineResults/*.csv")
    realtimeResults = glob.glob("../RealtimeResults/*.csv")
    springerResults = glob.glob("../SpringerSegmentationData/*.csv")
    pairedResults = findMatchingResults(offlineResults, springerResults, realtimeResults)


    for (offlineFilepath, springerFilepath, realtimeFilepath) in pairedResults:
        offlineData = np.genfromtxt(offlineFilepath, delimiter=',').astype(int)
        springerData = np.genfromtxt(springerFilepath, delimiter=',').astype(int)
        realtimeData = np.genfromtxt(realtimeFilepath, delimiter=',').astype(int)
        evaluate(offlineData, springerData, offlineFilepath, springerFilepath)

    correctOfflinePeakCount = 0
    incorrectOfflinePeakCount = 0
    correctOfflineLabelCount = 0
    incorrectOfflineLabelCount = 0
    totalOfflineLabelCount = 0
    totalOfflinePeakCount = 0
    totalOfflinePeakEst = 0

    for (offlineFilepath, springerFilepath, realtimeFilepath) in pairedResults:
        springerData = np.genfromtxt(springerFilepath, delimiter=',').astype(int)
        realtimeData = np.genfromtxt(realtimeFilepath, delimiter=',').astype(int)
        evaluate(realtimeData, springerData, realtimeFilepath, springerFilepath)

    print("-------------------------------------------------------------------------------")
    print("")
    print("OVERALL PERFORMANCE")
    print("===================")
    print("")
    print("Correct Peak Count:\t\t\t{0}\t({1:.01f}%)".format(correctOfflinePeakCount, (correctOfflinePeakCount/totalOfflinePeakCount)*100))
    print("Incorrect Peak Count:\t\t\t{0}\t({1:.01f}%)".format(incorrectOfflinePeakCount, (incorrectOfflinePeakCount/totalOfflinePeakCount)*100))
    print("Correct Label Count:\t\t\t{0}\t({1:.01f}%)".format(correctOfflineLabelCount, (correctOfflineLabelCount/totalOfflineLabelCount)*100))
    print("Incorrect Label Count:\t\t\t{0}\t({1:.01f}%)".format(incorrectOfflineLabelCount, (incorrectOfflineLabelCount/totalOfflineLabelCount)*100))
    print("Total Number of Peaks Estimate:\t\t{0}\t({1:+.0f})".format(totalOfflinePeakEst, totalOfflinePeakEst-totalOfflinePeakCount))
    print("Total Number of Peaks Count:\t\t{0}".format(totalOfflinePeakCount))
    print("-------------------------------------------------------------------------------")





def findMatchingResults(resultsA, resultsB, resultsC):
    out = []
    for filepathA in resultsA:
        basenameA = os.path.basename(filepathA)
        for filepathB in resultsB:
            basenameB = os.path.basename(filepathB)
            for filepathC in resultsC:
                basenameC = os.path.basename(filepathC)
                if  basenameA == basenameB and basenameA == basenameC:
                    out.append((filepathA, filepathB, filepathC))
    return out

if __name__ == "__main__":
    main()
