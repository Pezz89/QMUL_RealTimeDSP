#!/usr/bin/env python

from __future__ import division
import pysndfile
import scipy.signal as signal
import scipy.stats as stats
import glob
import pdb
import numpy as np
import matplotlib.pyplot as plt
from numpy.lib.stride_tricks import as_strided

plotFigures = True


def main():
    # Get file names for all PCG data
    PCGFiles = glob.glob("../validation_dataset/*.wav")
    # For each PCG data file...
    for filepath in PCGFiles:
        # Read audio in
        sndfile = pysndfile.PySndfile(filepath, 'r')
        data = x = sndfile.read_frames()
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
        x /= np.max(np.abs(x))

        # Split signal into 0.02 second overlapping grains with 0.01 second
        # overlap
        hopSize = int(0.01 * fs)
        N = winSize = int(0.02 * fs)

        # Window signal
        x = rolling_window(x, winSize, hopSize)

        # Calculate average shannon energy
        E_s = (-1/N) * np.sum((x**2)*np.log(x**2), axis=1)
        # Set nan values created by a log(0) to 0
        E_s[np.isnan(E_s)] = 0

        # Normalise average shannon energy
        mE_s = np.mean(E_s)
        sE_s = np.std(E_s)
        Pa = (E_s - mE_s)/sE_s

        threshold = 0.5

        if plotFigures:
            plt.subplot(2, 1, 1)
            x = np.linspace(0, data.shape[-1], Pa.size)
            plt.plot(x, Pa)
            plt.xlim([fs*1, fs*4])
            plt.axhline(threshold, linestyle='--', color='g')

            plt.subplot(2, 1, 2)
            plt.plot(data)
            plt.xlim([fs*1, fs*4])
            plt.xlabel('Time (samples)')
            plt.show()
        pdb.set_trace()



        # Calculate avergae shannon energy of each segment


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


if __name__ == "__main__":
    main()
