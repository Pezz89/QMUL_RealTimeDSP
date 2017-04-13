#!/usr/bin/env python

from __future__ import division, print_function
import pysndfile
import scipy.signal as signal
import scipy.stats as stats
import glob
import pdb
import numpy as np
import matplotlib.pyplot as plt
import peakutils
from peakutils.plot import plot as pplot
from numpy.lib.stride_tricks import as_strided
import numpy.ma as ma
import os
import errno

plotFigures1 = False
plotFigures2 = False
plotFigures3 = False
plotFigures4 = False
plotFigures5 = False

# Define the number of iterations to use when reducing the threshold for lost
# peaks search
reductionIterations = 20
# Define tolerances for heart sound classifications
c1 = 0.15
c2 = 0.3


class system:

    def __init__(self):
        self.calibratedMax = 0
        # Pretend that we know the sample rate...
        self.fs = 1000
        # Calculate window and hop size for windowing input signal
        self.hopSize = int(0.01 * self.fs)
        self.winSize = int(0.02 * self.fs)
        self.overlapFactor = np.ceil(self.winSize/self.hopSize)

        self.n = 103
        self.E_s = np.zeros(self.n, dtype=float)
        self.E_sPtr = 0
        PaBufferSize = np.round((5.0*self.fs) / self.hopSize)
        self.Pa = np.zeros(PaBufferSize)
        self.PaPtr = 0

        # Peak finding members
        # Boolean mask of peaks found in the last 5 seconds
        self.PaPeaks = np.zeros(PaBufferSize, dtype=bool)
        self.peakFound = False
        self.overshoot = False
        self.overshootCounter = False

        # Used for storing all peaks for offline analysis.
        self.allPeaks = np.array([])



    @staticmethod
    def moving_average(a, n=3) :
        ret = np.cumsum(a, dtype=float)
        ret[n:] = ret[n:] - ret[:-n]
        return ret[n - 1:] / n

    def calibrateSystem(self, sample):
        if sample > self.calibratedMax:
            self.calibratedMax = sample

    def calculateAvrShannonEnergy(self):
        newVal = False
        for ptr, pos in zip(self.bufferReadPtrs, self.bufferReadPtrsPos):
            if pos == self.winSize-1:
                inds = np.arange(ptr-self.winSize, ptr)%self.bufferArray.size
                sig = self.bufferArray[inds]
                a = (-1/self.winSize) * np.sum((sig**2)*np.log(sig**2))
                # Set nan values created by a log(0) to 0
                if not np.isfinite(a):
                    a = 0
                self.E_s[self.E_sPtr] = a
                # Increment pointer and wrap around
                self.E_sPtr += 1
                self.E_sPtr %= self.n

                # Normalise average shannon energy
                mE_s = np.mean(self.E_s)
                sE_s = np.std(self.E_s)
                a = (a - mE_s)/sE_s
                if not np.isfinite(a):
                    a = 0
                # Increment pointer location
                self.PaPtr += 1
                self.PaPtr %= self.Pa.size
                self.Pa[self.PaPtr] = a
                newVal = True
        return newVal

    def main(self):
        # Get file names for all PCG data
        PCGFiles = glob.glob("../validation_dataset/*.wav")

        # For each PCG data file...
        for filepath in PCGFiles:
            print("Analysing: {0}".format(filepath))
            # Read audio in
            sndfile = pysndfile.PySndfile(filepath, 'r')
            fs = sndfile.samplerate()
            sig = data = sndfile.read_frames()

            # Simulate chunks of audio as input to real-time system
            sig = self.rolling_window(sig, 2048, 2048)


            # Allocate circular buffer for calculating shannon energy using
            # overlapping windows
            self.bufferArray = np.zeros(self.winSize+((self.overlapFactor-1)*self.hopSize))
            self.bufferWritePtr = 0
            self.bufferReadPtrs = np.arange(self.overlapFactor, dtype=int)
            self.bufferReadPtrs *= -self.hopSize
            self.bufferReadPtrsPos = np.arange(self.overlapFactor, dtype=int)
            self.bufferReadPtrsPos *= -self.hopSize

            # Allocate interger to count the number os samples that have passed
            sampleCount = 0
            # For each chunk of audio in the signal
            for audio in sig:
                # For every sample in current audio block. Iterating over a numpy
                # array with a python for loop... There will be no good coding
                # practices from this point forward.
                for sample in audio:
                    # Increment sample counter by the number of samples that have been
                    # processed
                    sampleCount += 1

                    # If less than a second's worth of samples have been processed...
                    if sampleCount < fs:
                        # Use the current sample as part of system calibration
                        self.calibrateSystem(sample)
                        continue

                    # Normalise sample
                    sample *= self.calibratedMax

                    self.bufferArray[self.bufferWritePtr] = sample
                    newVal = self.calculateAvrShannonEnergy()

                    self.bufferWritePtr += 1
                    self.bufferWritePtr %= self.bufferArray.size
                    self.bufferReadPtrs += 1
                    self.bufferReadPtrs %= self.bufferArray.size
                    self.bufferReadPtrsPos += 1
                    self.bufferReadPtrsPos %= self.winSize

                    # Every time a new shannon energy value is calculated...
                    if newVal:
                        threshold = 0.5
                        # Calculate if it is a peak above the threshold
                        self.findNewPeaks(threshold)


            plt.plot(self.Pa)
            plt.plot(self.PaPeaks)
            plt.show()
            pdb.set_trace()

    def findNewPeaks(self, threshold):
        # If Pa values are above the threshold and a peak hasn't previously
        # been found
        if (self.Pa[self.PaPtr] >= threshold):
            if not self.peakFound:
                # If a start index has not been set for the current overshoot...
                if not self.overshoot:
                    # Store the start index of the overshoot
                    self.overshootCounter = 0
                    self.overshoot = True
                else:
                    self.overshootCounter += 1

                # search for a peak since the Pa values first raised above the
                # threshold

                # compute first order difference
                x = np.arange(self.PaPtr-self.overshootCounter, self.PaPtr+1)%self.Pa.size
                dy = np.diff(self.Pa[x])

                # propagate left and right values successively to fill all plateau pixels (0-value)
                zeros,=np.where(dy == 0)

                while len(zeros):
                    # add pixels 2 by 2 to propagate left and right value onto the zero-value pixel
                    zerosr = np.hstack([dy[1:], 0.])
                    zerosl = np.hstack([0., dy[:-1]])

                    # replace 0 with right value if non zero
                    dy[zeros]=zerosr[zeros]
                    zeros,=np.where(dy == 0)

                    # replace 0 with left value if non zero
                    dy[zeros]=zerosl[zeros]
                    zeros,=np.where(dy == 0)

                # find the peaks by using the first order difference
                peaks = (np.hstack([dy, 0.]) < 0.) & (np.hstack([0., dy]) > 0.)
                # If a peak is found, store location in boolean mask
                if np.any(peaks):
                    self.peakFound = True
                    self.PaPeaks[x] = peaks
                    self.overshoot = False
                else:
                    self.PaPeaks[self.PaPtr] = False
            else:
                self.PaPeaks[self.PaPtr] = False
        else:
            self.PaPeaks[self.PaPtr] = False
            self.peakFound = False
            self.overshoot = False
        return self.peakFound



    @staticmethod
    def saveResults(inputFile, results):
        path = '../OfflineResults/'

        # Create folder if it doesn't already exist
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

        inputBasename = os.path.basename(os.path.splitext(inputFile)[0])
        outputPath = os.path.join(path, inputBasename+"_segs.csv")
        np.savetxt(outputPath, results, fmt='%i', delimiter=",")

    @staticmethod
    def filterExtraPeaks(rejectionCandidates, candidatePeaks, allPeaks, x, fs, Pa, offset):
        if not np.ma.isMaskedArray(candidatePeaks):
            candidatePeaks = ma.array(candidatePeaks)
            candidatePeaks.mask = np.zeros(candidatePeaks.size, dtype=bool)
        for i, inds in enumerate(rejectionCandidates):
            # Get index location of peaks to potentially be rejected
            # Use data member to bypass the mask
            peakIndex1 = candidatePeaks.data[inds[0]]
            peakIndex2 = candidatePeaks.data[inds[1]]
            # Calculate time difference between indexes
            indexDiff = (x[peakIndex2]-x[peakIndex1])/fs
            # Calculate the ratio between first and second peaks amplitude
            # If time diference is less than 50ms...
            if indexDiff < 0.07:
                # If first peak is more than half the amplitude of the second,
                # reject the second peak, else reject the first
                if Pa[peakIndex1] > (Pa[peakIndex2] / 2):
                    candidatePeaks[inds[1]] = ma.masked
                else:
                    candidatePeaks[inds[0]] = ma.masked
            else:
                # If the first peak's energy is higher than that of the second
                # peak
                if Pa[peakIndex1] > Pa[peakIndex2]:
                    # Calculate mean and variance of every 2nd interval
                    # before the current interval
                    prevPeaksMask = allPeaks.mask.copy()
                    # TODO: Deal with this edge case
                    # if not isinstance(xx, np.ndarray):
                    #     prevPeaksMask = np.array([prevPeaksMask])
                    # Mask any peak indexes beyond and including the current interval
                    prevPeaksMask[offset+inds[1]:] = True
                    newPeakDiff = np.diff(allPeaks[~prevPeaksMask])
                    # Create array of all previous second intervals
                    secondIntervals = newPeakDiff[1-(newPeakDiff.size % 2)::2]
                    # Get last calculated interval
                    if not np.any(secondIntervals.data):
                        candidatePeaks[inds[1]] = ma.masked
                        continue

                    lastInterval = secondIntervals[-1]

                    # Sperate last calculated interval from all other intervals
                    secondIntervals = secondIntervals[:-1]
                    pDMean = np.mean(secondIntervals)
                    pDVar = np.var(secondIntervals)

                    # If current interval is more or less than the mean +/- the
                    # variance, remove first peak, else remove the second peak
                    if (lastInterval > pDMean + pDVar) or (lastInterval < pDMean - pDVar):
                        candidatePeaks[inds[0]] = ma.masked
                    else:
                        candidatePeaks[inds[1]] = ma.masked

                else:
                    # Else, reject the first peak
                    candidatePeaks[inds[0]] = ma.masked
        return candidatePeaks


    @staticmethod
    def getLastUnmaskedCandidate(rejectionCandidates, peaks, i):
        i -= 1
        while i > -1:
            if np.any(~peaks[rejectionCandidates[i]].mask):
                a = peaks[rejectionCandidates[i]]
                b = rejectionCandidates[i][~a.mask]
                return b[-1]
            i -= 1
        return -1

    @staticmethod
    def rolling_window(a, window, hopSize):
        # This function was adapted from: http://stackoverflow.com/questions/4936620/using-strides-for-an-efficient-moving-average-filter
        if window < 1:
            raise ValueError, "`window` must be at least 1."
        if window > a.shape[-1]:
            raise ValueError, "`window` is too long."
        numWindows = np.ceil((a.shape[-1]-hopSize)/hopSize)
        shape = a.shape[:-1] + (numWindows, window)
        strides = (a.strides[0]*hopSize,) + a.strides[1:-1] + (a.strides[-1],)

        return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)


    # This function was modified from the following repository: https://bitbucket.org/lucashnegri/peakutils
    @staticmethod
    def indexes(y, thres=0.3, min_dist=1):
        """Peak detection routine.

        Finds the numeric index of the peaks in *y* by taking its first order difference. By using
        *thres* and *min_dist* parameters, it is possible to reduce the number of
        detected peaks. *y* must be signed.

        Parameters
        ----------
        y : ndarray (signed)
            1D amplitude data to search for peaks.
        thres : float between [0., 1.]
            Normalized threshold. Only the peaks with amplitude higher than the
            threshold will be detected.
        min_dist : int
            Minimum distance between each detected peak. The peak with the highest
            amplitude is preferred to satisfy this constraint.

        Returns
        -------
        ndarray
            Array containing the numeric indexes of the peaks that were detected
        """
        if isinstance(y, np.ndarray) and np.issubdtype(y.dtype, np.unsignedinteger):
            raise ValueError("y must be signed")

        min_dist = int(min_dist)

        # compute first order difference
        dy = np.diff(y)

        # propagate left and right values successively to fill all plateau pixels (0-value)
        zeros,=np.where(dy == 0)

        while len(zeros):
            # add pixels 2 by 2 to propagate left and right value onto the zero-value pixel
            zerosr = np.hstack([dy[1:], 0.])
            zerosl = np.hstack([0., dy[:-1]])

            # replace 0 with right value if non zero
            dy[zeros]=zerosr[zeros]
            zeros,=np.where(dy == 0)

            # replace 0 with left value if non zero
            dy[zeros]=zerosl[zeros]
            zeros,=np.where(dy == 0)

        # find the peaks by using the first order difference
        peaks = np.where((np.hstack([dy, 0.]) < 0.)
                        & (np.hstack([0., dy]) > 0.)
                        & (y > thres))[0]

        if peaks.size > 1 and min_dist > 1:
            highest = peaks[np.argsort(y[peaks])][::-1]
            rem = np.ones(y.size, dtype=bool)
            rem[peaks] = False

            for peak in highest:
                if not rem[peak]:
                    sl = slice(max(0, peak - min_dist), peak + min_dist + 1)
                    rem[sl] = True
                    rem[peak] = False

            peaks = np.arange(y.size)[~rem]

        # Find threshold boundaries, between which, y is over the threshold
        boundaries = np.where(np.diff(np.array(y > thres, dtype=int))>0)[0]
        # If there are no boundaries then y is never over the threshold
        if boundaries.shape[0] < 2:
            return np.array([], dtype=int)
        # Make start-end pairs of boundaries
        boundaries = rolling_window(boundaries, 2, 1)
        out = np.zeros(boundaries.shape[0], dtype=int)
        # For each threshold boundary, select only the first peak
        for ind, indexes in enumerate(boundaries):
            peaks[(peaks > indexes[0]) & (peaks < indexes[1])]
            out[ind] = np.min(peaks[(peaks > indexes[0]) & (peaks < indexes[1])])

        return out

if __name__ == "__main__":
    a = system()
    a.main()

    #x = ((np.arange(Pa.size)*hopSize)+np.round(winSize/2)).astype(int)

    '''
            # Calculate timings for each analysis window
            # x = np.linspace(np.round(winSize/2), data.shape[-1]-np.round(winSize/2), Pa.size)

            threshold = 0.5

            if plotFigures1:
                plt.subplot(2, 1, 1)
                plt.plot(x, Pa)
                plt.xlim([fs*1, fs*4])
                plt.axhline(threshold, linestyle='--', color='g')

                plt.subplot(2, 1, 2)
                plt.plot(data)
                plt.xlim([fs*1, fs*4])
                plt.xlabel('Time (samples)')
                plt.show()

            # Create a masked array of indexes for each peak found in Pa
            peaks = ma.array(indexes(Pa, thres=threshold))
            peaks.mask = np.zeros(peaks.size, dtype=bool)

            if plotFigures2:
                pplot(x, Pa, peaks[~peaks.mask])
                plt.xlim([fs*1, fs*4])
                plt.xlabel('Time (samples)')
                plt.axhline(threshold, linestyle='--', color='g')
                plt.show()

            # Calculate the interval between adjacent peaks
            peakDiff = np.diff(peaks)
            pDMean = np.mean(peakDiff)
            pDStd = np.std(peakDiff)
            # Calculate high and low interval limits using mean and standard
            # deviation
            lowIntervalLim = pDMean - pDStd
            highIntervalLim = pDMean + pDStd
            rejectionCandidates = np.where(peakDiff < lowIntervalLim)[0]
            # Flip array vertially
            rejectionCandidates = rejectionCandidates[np.newaxis].T
            # Create pairs of indexes for peaks to be compared
            rejectionCandidates = np.hstack((rejectionCandidates, rejectionCandidates+1))

            peaks = filterExtraPeaks(rejectionCandidates, peaks, peaks, x, fs, Pa, 0)


            if plotFigures3:
                pplot(x, Pa, peaks[~peaks.mask])
                plt.xlim([fs*1, fs*4])
                plt.xlabel('Time (samples)')
                plt.axhline(threshold, linestyle='--', color='g')
                plt.show()

            inclusionCandidates = np.where(peakDiff > highIntervalLim)[0]
            # Flip array vertially
            inclusionCandidates = inclusionCandidates[np.newaxis].T
            inclusionCandidates = np.hstack((inclusionCandidates, inclusionCandidates+1))


            newPeaks = []
            for i, inds in enumerate(inclusionCandidates):
                # Get index location of boundaries in which peaks may have been
                # lost. Use data member to bypass the mask
                peakIndex1 = peaks.data[inds[0]]
                peakIndex2 = peaks.data[inds[1]]

                reductionAmount = (threshold - np.min(Pa)) / reductionIterations
                reducedThreshold = threshold - reductionAmount
                while reducedThreshold > np.min(Pa):
                    # Create threshold to be iterativey reduced until peaks are found
                    reducedThreshold -= reductionAmount
                    # Get all peaks that aren't currently masked
                    foundPeakInds = indexes(Pa[peakIndex1:peakIndex2], thres=reducedThreshold)
                    if np.any(foundPeakInds):
                        # Offset peaks by starting search index to find their global index
                        foundPeakInds += peakIndex1
                        foundPeakInds = np.append(peakIndex1, foundPeakInds)
                        foundPeakDiff = np.diff(foundPeakInds)
                        foundRejectionCandidates = np.where(foundPeakDiff < lowIntervalLim)[0]
                        # Flip array vertially
                        foundRejectionCandidates = foundRejectionCandidates[np.newaxis].T
                        # Create pairs of indexes for peaks to be compared
                        foundRejectionCandidates = np.hstack((foundRejectionCandidates, foundRejectionCandidates+1))
                        newPeakInds = filterExtraPeaks(foundRejectionCandidates, foundPeakInds, peaks, x, fs, Pa, inds[0])
                        newPeaks.append(newPeakInds)
                        break
            for inds in newPeaks:
                peaks = ma.append(peaks, inds)
                peaks = ma.unique(peaks)

            # # Calculate the interval between adjacent peaks
            # peakDiff = np.diff(peaks)
            # pDMean = np.mean(peakDiff)
            # pDStd = np.std(peakDiff)
            # # Calculate high and low interval limits using mean and standard
            # # deviation
            # lowIntervalLim = pDMean - pDStd
            # highIntervalLim = pDMean + pDStd
            # rejectionCandidates = np.where(peakDiff < lowIntervalLim)[0]
            # # Flip array vertially
            # rejectionCandidates = rejectionCandidates[np.newaxis].T
            # # Create pairs of indexes for peaks to be compared
            # rejectionCandidates = np.hstack((rejectionCandidates, rejectionCandidates+1))

            # peaks = filterExtraPeaks(rejectionCandidates, peaks, peaks, x, fs, Pa, 0)

            # Get all valid peaks
            if np.any(peaks.mask):
                peaks = peaks[~peaks.mask]
            # Calculate the difference between all peaks
                peakDiff = np.diff(peaks)
            # Find the largest interval between all peaks in the last 20 seconds
            diastolicPeriod = np.max(peakDiff)
            # Get index value of last peak
            endPeak = peaks[-1]
            # Find all peaks within the last 20 seconds
            peaksInRangeMask = x[np.abs(peaks-endPeak)] < 20*fs
            peaksInRange = peaks[peaksInRangeMask]

            # Calculate the difference between all peaks
            peakDiff = np.diff(peaks)
            # Find the largest interval between all peaks in the last 20 seconds
            peaksInRangeDiff = np.diff(peaksInRange)
            diastolicPeriod = np.max(peaksInRangeDiff)

            # Find all diastolic periods
            diastolicPeriodsMask = (peakDiff <= diastolicPeriod + diastolicPeriod * c2) & (peakDiff >= diastolicPeriod - diastolicPeriod * c2)
            otherPeriods = peakDiff[~diastolicPeriodsMask]
            systolicPeriod = np.median(otherPeriods)
            classification = np.zeros(peakDiff.shape)
            # Handle situation where no systolic segments are detected. Happens in
            # poor recordings...
            if systolicPeriod:
                systolicPeriodsMask = (peakDiff <= systolicPeriod + systolicPeriod * c1) & (peakDiff >= systolicPeriod - systolicPeriod * c1)
                classification[systolicPeriodsMask] = 1
            classification[diastolicPeriodsMask] = 2

            classification = classification[np.newaxis].T
            times = np.round(x[peaks]/2)[np.newaxis].T
            results = np.hstack((times[:-1], classification)).astype(int)

            if plotFigures4:
                pplot(x, Pa, peaks[~peaks.mask])
                plt.xlim([fs*1, fs*4])
                plt.xlabel('Time (samples)')
                plt.axhline(threshold, linestyle='--', color='g')
                plt.show()

            if plotFigures5:
                p = peaks[:-1]
                plt.plot(x, Pa)
                plt.plot(x[p], classification, "x")
                plt.show()
            saveResults(filepath, results)
    '''
