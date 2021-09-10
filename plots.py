import numpy as np
import matplotlib.pyplot as plt
import os
import h5py
from scipy.signal import argrelmax


def plot_peak_freqs(path, axes, Nmaxes = 3):
    # Nmaxes number of maximum frequencies on plot

    stim_freq = []

    peak_freqs = []
    for _ in range(Nmaxes):
        peak_freqs.append([])

    counter = 1
    for file in os.listdir(path):
        if file.split(".")[-1] != "hdf5":
            continue

        #print(file)
        filepath = path + file
        with h5py.File(filepath, 'r') as h5file:
            stim_freq.append(counter) # read from h5file !!!!

            fft_amplitudes = h5file["extracellular/electrode_1/firing/processing/msach/fft_amplitudes"][:]
            fft_freqs = h5file["extracellular/electrode_1/firing/processing/msach/fft_freqs"][:]
            peak_freq_args = argrelmax(fft_amplitudes, mode="wrap")[0]

            fft_freqs_peak = fft_freqs[peak_freq_args]
            peak_freq_args = np.argsort( -fft_amplitudes[peak_freq_args] )

            #plt.plot(fft_freqs, fft_amplitudes)
            for peak_freq_idx in range(Nmaxes):
                peak_freqs[peak_freq_idx].append(fft_freqs_peak[peak_freq_args[peak_freq_idx]])

        counter += 1

    stim_freq = np.asarray(stim_freq)
    for idx in range(Nmaxes):
        axes.plot(stim_freq, peak_freqs[idx], label=str(idx))

    axes.legend()

def plot_mean_CV(path, axes, parameter="mean_cv"):
    stim_freq = []
    vals = []
    counter = 1
    for file in os.listdir(path):
        if file.split(".")[-1] != "hdf5":
            continue

        print(file)
        filepath = path + file
        with h5py.File(filepath, 'r') as h5file:
            stim_freq.append(counter)  # read from h5file !!!!
            # .format(parameter)
            val = h5file["extracellular/electrode_1/firing/processing/msach/{}".format(parameter)][()]
            vals.append(val)
            counter += 1

    print(len(vals))
    axes.plot(stim_freq, vals, label=parameter)
    axes.legend()

if __name__ == '__main__':
    path = "./Results/"
    fig, axes = plt.subplots(ncols=3, figsize=(10, 5))
    plot_peak_freqs(path, axes[0])
    plot_mean_CV(path, axes[1], parameter="mean_cv")
    plot_mean_CV(path, axes[2], parameter="mean_correlation")

    plt.show()


