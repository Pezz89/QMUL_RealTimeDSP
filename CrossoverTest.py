#!/usr/bin/env python

import numpy as np
import pdb
import scipy.signal as signal
import matplotlib.pyplot as plt

def main():
    # Sampling frequency
    fs = 44100
    # Cutoff frequency in hertz
    cutoff_freq = 5000.

    # Calculate ratio between cutoff frequency and sampling rate
    wc = cutoff_freq/fs

    # Deifine Q as the square root of 2
    q=np.sqrt(2.0)

    # Warp the frequency to convert from continuous to discrete time cutoff
    wd1 = 1.0 / np.tan(np.pi*wc)

    # Calculate coefficients from equation
    b0 = 1.0 / (1.0 + q*wd1 + wd1**2)
    b1 = 2*b0
    b2 = b0
    a1 = -2.0 * (wd1**2 - 1.0) * b0
    a2 = (1.0 - q*wd1 + wd1**2) * b0

    b_conv = np.convolve([b0, b1, b2], [b0, b1, b2])
    a_conv = np.convolve([1, a1, a2], [1, a1, a2])
    w, h = signal.freqz([b0, b1, b2], [1, a1, a2], plot = lambda w, h: plt.plot((fs * 0.5 / np.pi) * w, 20 * np.log10(abs(h))))
    w, h = signal.freqz(b_conv, a_conv, plot = lambda w, h: plt.plot((fs * 0.5 / np.pi) * w, 20 * np.log10(abs(h))))
    pdb.set_trace()
    b0 = b0*wd1**2;
    b1 = -b1*wd1**2;
    b2 = b2*wd1**2;
    b_conv = np.convolve([b0, b1, b2], [b0, b1, b2])
    a_conv = np.convolve([1, a1, a2], [1, a1, a2])
    w, h = signal.freqz([b0, b1, b2], [1, a1, a2], plot = lambda w, h: plt.plot((fs * 0.5 / np.pi) * w, 20 * np.log10(abs(h))))
    w, h = signal.freqz(b_conv, a_conv, plot = lambda w, h: plt.plot((fs * 0.5 / np.pi) * w, 20 * np.log10(abs(h))))


    plt.axvline(cutoff_freq, color='r')
    plt.show()


if __name__ == "__main__":
    main()
