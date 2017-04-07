#!/usr/bin/env python

from __future__ import division
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

plotFigures1 = False
plotFigures2 = False
plotFigures3 = True


def main():
    # Get file names for all PCG data
    PCGFiles = glob.glob("../validation_dataset/*.wav")
    # For each PCG data file...
    for filepath in PCGFiles:
        # Read audio in
        sndfile = pysndfile.PySndfile(filepath, 'r')
        data = sig = sndfile.read_frames()
        fs = sndfile.samplerate()
        # Downsampling is no longer neccesary as recording have already been
        # downsampled in the validation set.
        '''
        # Calculate decimation factor
        q = np.floor(fs/2205)
        # Decimate to samplerate of 2205Hz, applying zero-phase shift
        x = signal.decimate(x, q, zero_phase=True)
        '''

        # Normalise signal
        sig /= np.max(np.abs(sig))

        # Split signal into 0.02 second overlapping grains with 0.01 second
        # overlap
        hopSize = int(0.01 * fs)
        N = winSize = int(0.02 * fs)

        # Window signal
        sig = rolling_window(sig, winSize, hopSize)

        # Calculate average shannon energy
        E_s = (-1/N) * np.sum((sig**2)*np.log(sig**2), axis=1)
        # Set nan values created by a log(0) to 0
        E_s[np.isnan(E_s)] = 0

        # Normalise average shannon energy
        mE_s = np.mean(E_s)
        sE_s = np.std(E_s)
        Pa = (E_s - mE_s)/sE_s
        x = ((np.arange(Pa.size)*hopSize)+np.round(winSize/2)).astype(int)

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

        if plotFigures2:
            pplot(x, Pa, peaks)
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

        for i, inds in enumerate(rejectionCandidates):
            # If a candidate has been previously masked...

            '''
            if ma.is_masked(peaks[inds[0]]):
                # Find the last unmasked candidate
                inds[0] = peaks.data[inds[0]]
                getLastUnmaskedCandidate(rejectionCandidates, peaks, i)

                # If there are none less than the interval limit then continue
                if (inds[0] < 0) or (inds[0] > lowIntervalLim):
                    continue
            '''
            # Get index location of peaks to potentially be rejected
            peakIndex1 = peaks.data[inds[0]]
            peakIndex2 = peaks.data[inds[1]]
            # Calculate time difference between indexes
            indexDiff = (x[peakIndex2]-x[peakIndex1])/fs
            # Calculate the ratio between first and second peaks amplitude
            # If time diference is less than 50ms...
            if indexDiff < 0.05:
                # If first peak is more than half the amplitude of the second,
                # reject the second peak, else reject the first
                if Pa[peakIndex1] > (Pa[peakIndex2] / 2):
                    peaks[inds[1]] = ma.masked
                else:
                    peaks[inds[0]] = ma.masked
            else:
                # If the first peak's energy is higher than that of the second
                # peak
                if Pa[peakIndex1] > Pa[peakIndex2]:
                    # Calculate mean and variance of every 2nd interval
                    # before the current interval
                    prevPeaksMask = peaks.mask.copy()
                    # Mask any peak indexes beyond and including the current interval
                    prevPeaksMask[inds[1]:] = True
                    newPeakDiff = np.diff(peaks[~prevPeaksMask])
                    # Create array of all previous second intervals
                    secondIntervals = newPeakDiff[1-(newPeakDiff.size % 2)::2]
                    # Get last calculated interval
                    lastInterval = secondIntervals[-1]
                    # Sperate last calculated interval from all other intervals
                    secondIntervals = secondIntervals[:-1]
                    pDMean = np.mean(secondIntervals)
                    pDVar = np.var(secondIntervals)

                    # If current interval is more or less than the mean +/- the
                    # variance, remove first peak, else remove the second peak
                    if (lastInterval > pDMean + pDVar) or (lastInterval < pDMean - pDVar):
                        peaks[inds[0]] = ma.masked
                    else:
                        peaks[inds[1]] = ma.masked

                else:
                    # Else, reject the first peak
                    peaks[inds[0]] = ma.masked

        rejectionCandidates = np.where(peakDiff > highIntervalLim)[0]

        if plotFigures3:
            pplot(x, Pa, peaks[~peaks.mask])
            plt.xlim([fs*1, fs*4])
            plt.xlabel('Time (samples)')
            plt.axhline(threshold, linestyle='--', color='g')
            plt.show()
        pdb.set_trace()


        # Calculate avergae shannon energy of each segment

def getLastUnmaskedCandidate(rejectionCandidates, peaks, i):
    i -= 1
    while i > -1:
        if np.any(~peaks[rejectionCandidates[i]].mask):
            a = peaks[rejectionCandidates[i]]
            b = rejectionCandidates[i][~a.mask]
            return b[-1]
        i -= 1
    return -1

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

    # Find threshold boundaries
    boundaries = np.where(np.diff(np.array(y > thres, dtype=int))>0)[0]
    # Make start-end pairs of boundaries
    boundaries = rolling_window(boundaries, 2, 1)
    out = np.zeros(boundaries.shape[0], dtype=int)
    # For each threshold boundary, select only the first peak
    for ind, indexes in enumerate(boundaries):
        out[ind] = np.min(peaks[(peaks > indexes[0]) & (peaks < indexes[1])])

    return out

if __name__ == "__main__":
    main()
