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
    print("Correct Peak Count: {0}".format(correctPeakCount))
    print("Incorrect Peak Count: {0}".format(incorrectPeakCount))


    # For peaks that were correctly identified check which ones were
    # correctly labeled



def findMatchingResults(resultsA, resultsB):
    out = []
    for filepathA in resultsA:
        for filepathB in resultsB:
            if os.path.basename(filepathA) == os.path.basename(filepathB):
                out.append((filepathA, filepathB))
    return out

if __name__ == "__main__":
    main()
