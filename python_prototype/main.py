#!/usr/bin/env python

from __future__ import division
import pysndfile
import scipy.signal as signal
import glob
import pdb
import numpy as np
from numpy.lib.stride_tricks import as_strided


def main():
    # Get file names for all PCG data
    PCGFiles = glob.glob("../validation_dataset/*.wav")
    # For each PCG data file...
    for filepath in PCGFiles:
        # Read audio in
        sndfile = pysndfile.PySndfile(filepath, 'r')
        x = sndfile.read_frames()
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
        winSize = int(0.02 * fs)

        #rolling_window(x, winSize, hopSize)

        a = rolling_window(np.array([1,2,3,4,5,6,7,8,9]), 3, 2)
        pdb.set_trace()


        # Calculate avergae shannon energy of each segment


def rolling_window(a, window, hopSize):
    numWindows = np.ceil((a.shape[-1]-hopSize)/hopSize)
    shape = a.shape[:-1] + (numWindows, window)
    strides = (a.strides[0]*hopSize,) + a.strides[1:-1] + (a.strides[-1],)
    pdb.set_trace()

    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)


if __name__ == "__main__":
    main()
