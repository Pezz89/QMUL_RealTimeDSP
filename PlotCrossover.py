#!/usr/bin/env python

import numpy as np
import pdb
import scipy.signal as signal
import matplotlib.pyplot as plt

def main():
    # Sampling frequency
    fs = 44100
    # Cutoff frequency in hertz
    cutoff_freq = 1000.

    # Calculate ratio between cutoff frequency and sampling rate
    wc = cutoff_freq/fs

    # Deifine Q as the square root of 2
    q=np.sqrt(2.0)

    # Warp the frequency to convert from continuous to discrete time cutoff
    wd1 = 1.0 / np.tan(np.pi*wc)

    b = np.zeros(3)
    a = np.zeros(3)

    # Calculate coefficients from equation
    b[0] = 1.0 / (1.0 + q*wd1 + wd1**2)
    b[1] = 2*b[0]
    b[2] = b[0]
    a[0] = 1
    a[1] = -2.0 * (wd1**2 - 1.0) * b[0]
    a[2] = (1.0 - q*wd1 + wd1**2) * b[0]

    # Convolve coefficients with themselves to produce a 4th order
    # Linkwitz-Riley filter.
    b_conv = np.convolve(b, b)
    a_conv = np.convolve(a, a)

    dpi = 100
    fig = plt.figure(figsize=(1300/dpi, 1000/dpi), dpi=dpi)

    # Plot dB magnitude response of 2nd order Butterworth low-pass filter
    w, bwlp_h = signal.freqz(b, a, plot = lambda w, h: plt.plot((fs * 0.5 / np.pi) * w, 20 * np.log10(abs(h)), color='b', label = 'Exepected Low-pass Value'))

    # Convert low-pass to high pass
    b[0] = b[0]*wd1**2;
    b[1] = -b[1]*wd1**2;
    b[2] = b[2]*wd1**2;

    # Convolve coefficients with themselves to produce a 4th order
    # Linkwitz-Riley filter.
    b_conv = np.convolve(b, b)
    a_conv = np.convolve(a, a)

    # Plot dB magnitude response of 2nd order Butterworth high-pass filter
    w, bwhp_h = signal.freqz(b, a, plot = lambda w, h: plt.plot((fs * 0.5 / np.pi) * w, 20 * np.log10(abs(h)), color='g', label='Expected High-pass Value'))
    # Define values manually measured from ossiloscope for low-pass
    manual_plt_x = np.array([300, 500, 750, 1000, 1250, 1500, 1750, 2000, 2250, 2500, 3000, 4000])
    manual_plt_y = 20 * np.log10(abs(np.array([1.0, 0.957, 0.864, 0.707, 0.538, 0.414, 0.318, 0.256, 0.203, 0.169, 0.122, 0.076])))
    # Plot values against low pass curve to check accuracy
    plt.plot(manual_plt_x, manual_plt_y, marker='x',mew=2, ms=10, linestyle='None', color='r', label='Low-pass Manual Measurement')

    # Define values manually measured from ossiloscope for high-pass
    manual_plt_x = np.array([300, 500, 750, 1000, 1250, 1500, 1750, 2000, 2250, 2500, 3000, 4000])
    manual_plt_y = 20 * np.log10(abs(np.array([0.082, 0.231, 0.467, 0.706, 0.814, 0.879, 0.913, 0.938, 0.950, 1.000, 1.000, 1.000])))
    # Plot values against high pass curve to check accuracy
    plt.plot(manual_plt_x, manual_plt_y, marker='x',mew=2, ms=10, linestyle='None', color='r', label='High-pass Manual Measurement')

    # Display cutoff frequency
    plt.axvline(cutoff_freq, color='r', linestyle='--')
    plt.ylim(-25, 5)
    plt.xlim(0, manual_plt_x[-1]+500)

    # Get current tick locations and append 271 to this array
    x_ticks = np.append(plt.xticks()[0], cutoff_freq)

    # Set xtick locations to the values of the array `x_ticks`
    plt.xticks(x_ticks)
    plt.legend(
        loc='upper center',
        bbox_to_anchor=(0.5, 1.10),
        ncol=2,
        fancybox=True,
        shadow=True
    )
    plt.grid(True)
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Magnitude (log dB)")
    fig.savefig("./BWFreqResp.png")

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

    b = np.zeros(3)
    a = np.zeros(3)
    # Calculate coefficients from equation
    b[0] = 1.0 / (1.0 + q*wd1 + wd1**2)
    b[1] = 2*b[0]
    b[2] = b[0]
    a[0] = 1
    a[1] = -2.0 * (wd1**2 - 1.0) * b[0]
    a[2] = (1.0 - q*wd1 + wd1**2) * b[0]

    # Convolve coefficients with themselves to produce a 4th order
    # Linkwitz-Riley filter.
    b_conv = np.convolve(b, b)
    a_conv = np.convolve(a, a)
    dpi = 100
    fig = plt.figure(figsize=(1300/dpi, 1000/dpi), dpi=dpi)

    # Plot dB magnitude response of 2nd order Butterworth low-pass filter
    w, bwlp_h = signal.freqz(b, a, plot = lambda w, h: plt.plot((fs * 0.5 / np.pi) * w, 20 * np.log10(abs(h)), color='r', label = 'Butterworth'))
    # Plot dB magnitude response of 4th order Linkwitz-Riley low-pass filter
    w, lrlp_h = signal.freqz(b_conv, a_conv, plot = lambda w, h: plt.plot((fs * 0.5 / np.pi) * w, 20 * np.log10(abs(h)), color='b'))

    # Convert low-pass to high pass
    b[0] = b[0]*wd1**2;
    b[1] = -b[1]*wd1**2;
    b[2] = b[2]*wd1**2;

    # Convolve coefficients with themselves to produce a 4th order
    # Linkwitz-Riley filter.
    b_conv = np.convolve(b, b)
    a_conv = np.convolve(a, a)
    # Plot dB magnitude response of 2nd order Butterworth high-pass filter
    w, bwhp_h = signal.freqz(b, a, plot = lambda w, h: plt.plot((fs * 0.5 / np.pi) * w, 20 * np.log10(abs(h)), color='r'))
    # Plot dB magnitude response of 4th order Linkwitz-Riley high-pass filter
    w, lrhp_h = signal.freqz(b_conv, a_conv, plot = lambda w, h: plt.plot((fs * 0.5 / np.pi) * w, 20 * np.log10(abs(h)), color='b', label='Linkwitz-Riley'))

    plt.plot((fs * 0.5 / np.pi) * w, 20 * np.log10(abs(bwlp_h)+abs(bwhp_h)), linewidth=2.0, linestyle=':', color='r', label='Btrwrth Crossover Gain')
    plt.plot((fs * 0.5 / np.pi) * w, 20 * np.log10(abs(lrlp_h)+abs(lrhp_h)), linewidth=2.0, linestyle=':', color='b', label='Lw-Rl Crossover Gain')

    # Display cutoff frequency
    plt.axvline(cutoff_freq, color='r', linestyle='--')
    plt.ylim(-25, 5)
    plt.xlim(0, 15000)

    # Get current tick locations and append 271 to this array
    x_ticks = np.append(plt.xticks()[0], cutoff_freq)

    # Set xtick locations to the values of the array `x_ticks`
    plt.xticks(x_ticks)
    plt.legend(
        loc='upper center',
        bbox_to_anchor=(0.5, 1.10),
        ncol=2,
        fancybox=True,
        shadow=True
    )
    plt.annotate(
        '-3dB',
        xy=(5000, -3),
        xytext=(5000, -3),
        xycoords='data',
        textcoords='data'
    )
    plt.annotate(
        '-6dB',
        xy=(5000, -6),
        xytext=(5000, -6),
        xycoords='data',
        textcoords='data'
    )
    plt.annotate(
        '0dB',
        xy=(5000, 0),
        xytext=(5000, 0),
        xycoords='data',
        textcoords='data'
    )
    plt.annotate(
        '+3dB',
        xy=(5000, 3),
        xytext=(5000, 3),
        xycoords='data',
        textcoords='data'
    )
    plt.grid(True)
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Magnitude (log dB)")
    fig.savefig("./XOverFreqResp.png")


if __name__ == "__main__":
    main()
