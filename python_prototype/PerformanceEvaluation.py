#!/usr/bin/env python

from __future__ import division, print_function
import numpy as np
import glob
import os
import pdb

def main():

    offlineResults = glob.glob("../OfflineResults/*.csv")
    springerResults = glob.glob("../SpringerSegmentationData/*.csv")
    pairedResults = findMatchingResults(offlineResults, springerResults)

    fs = 1000

    correctPeakCount = 0
    incorrectPeakCount = 0
    correctLabelCount = 0
    incorrectLabelCount = 0
    totalLabelCount = 0
    totalPeakCount = 0
    totalPeakEst = 0
    for (offlineFilepath, springerFilepath) in pairedResults:
        offlineData = np.genfromtxt(offlineFilepath, delimiter=',').astype(int)
        springerData = np.genfromtxt(springerFilepath, delimiter=',').astype(int)

        # Check which peaks were correctly identified in time by finding a
        # corresponding peak in springer data within 100ms
        springerTimes = springerData[:, 0]
        offlineTimes = offlineData[:, 0]
        c = np.isclose(springerTimes, offlineTimes[np.newaxis].T, 0, 0.1*fs)
        results = np.any(c, axis = 0)
        cpc = np.count_nonzero(results)
        correctPeakCount += cpc
        ipc = np.size(results) - np.count_nonzero(results)
        incorrectPeakCount += ipc
        tpe = offlineTimes.size
        totalPeakEst += tpe
        tpc = springerTimes.size
        totalPeakCount += tpc

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
        correctLabelCount += clc
        ilc = np.size(correctLabels) - np.count_nonzero(correctLabels)
        incorrectLabelCount += ilc
        tlc = springerLabels.size
        totalLabelCount += tlc
        print("-------------------------------------------------------------------------------")
        print("Filepath: {0}".format(offlineFilepath))
        print("Correct Peak Count:\t\t\t{0}\t(%{1:.01f})".format(cpc, (cpc/tpc)*100))
        print("Incorrect Peak Count:\t\t\t{0}\t(%{1:.01f})".format(ipc, (ipc/tpc)*100))
        print("Correct Label Count:\t\t\t{0}\t(%{1:.01f})".format(clc, (clc/tlc)*100))
        print("Incorrect Peak Count:\t\t\t{0}\t(%{1:.01f})".format(ilc, (ilc/tlc)*100))
        print("Total Number of Peaks Estimate:\t\t{0}\t({1:+.0f})".format(tpe, tpe-tpc))
        print("Total Number of Peaks Count:\t\t{0}".format(tpc))
        print("-------------------------------------------------------------------------------")

    print("-------------------------------------------------------------------------------")
    print("")
    print("OVERALL PERFORMANCE")
    print("===================")
    print("")
    print("Correct Peak Count:\t\t\t{0}\t(%{1:.01f})".format(correctPeakCount, (correctPeakCount/totalPeakCount)*100))
    print("Incorrect Peak Count:\t\t\t{0}\t(%{1:.01f})".format(incorrectPeakCount, (incorrectPeakCount/totalPeakCount)*100))
    print("Correct Label Count:\t\t\t{0}\t(%{1:.01f})".format(correctLabelCount, (correctLabelCount/totalLabelCount)*100))
    print("Incorrect Peak Count:\t\t\t{0}\t(%{1:.01f})".format(incorrectLabelCount, (incorrectLabelCount/totalLabelCount)*100))
    print("Total Number of Peaks Estimate:\t\t{0}\t({1:+.0f})".format(totalPeakEst, totalPeakEst-totalPeakCount))
    print("Total Number of Peaks Count:\t\t{0}".format(totalPeakCount))
    print("-------------------------------------------------------------------------------")





def findMatchingResults(resultsA, resultsB):
    out = []
    for filepathA in resultsA:
        for filepathB in resultsB:
            if os.path.basename(filepathA) == os.path.basename(filepathB):
                out.append((filepathA, filepathB))
    return out

if __name__ == "__main__":
    main()
