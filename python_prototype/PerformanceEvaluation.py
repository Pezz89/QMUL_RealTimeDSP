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
        correctPeakCount += np.count_nonzero(results)
        incorrectPeakCount += np.size(results) - np.count_nonzero(results)
        totalPeakEst += offlineTimes.size
        totalPeakCount += springerTimes.size

        # For peaks that were correctly identified check which ones were
        # correctly labeled
        springerLabels = springerData[:, 1]
        offlineLabels = offlineData[:, 1]

        inds = np.where(c)
        sL = springerLabels[inds[1]]
        oL = offlineLabels[inds[0]]
        correctLabels = sL == oL
        correctLabelCount += np.count_nonzero(correctLabels)
        incorrectLabelCount += np.size(correctLabels) - np.count_nonzero(correctLabels)
        totalLabelCount += springerLabels.size

    print("Correct Peak Count:\t\t\t{0}\t(%{1:.01f})".format(correctPeakCount, (correctPeakCount/totalPeakCount)*100))
    print("Incorrect Peak Count:\t\t\t{0}\t(%{1:.01f})".format(incorrectPeakCount, (incorrectPeakCount/totalPeakCount)*100))
    print("Correct Label Count:\t\t\t{0}\t(%{1:.01f})".format(correctLabelCount, (correctLabelCount/totalLabelCount)*100))
    print("Incorrect Peak Count:\t\t\t{0}\t(%{1:.01f})".format(incorrectLabelCount, (incorrectLabelCount/totalLabelCount)*100))
    print("Total Number of Peaks Estimate:\t\t{0}\t({1:+.0f})".format(totalPeakEst, totalPeakEst-totalPeakCount))
    print("Total Number of Peaks Count:\t\t{0}".format(totalPeakCount))





def findMatchingResults(resultsA, resultsB):
    out = []
    for filepathA in resultsA:
        for filepathB in resultsB:
            if os.path.basename(filepathA) == os.path.basename(filepathB):
                out.append((filepathA, filepathB))
    return out

if __name__ == "__main__":
    main()
