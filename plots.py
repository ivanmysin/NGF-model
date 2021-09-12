import numpy as np
import matplotlib.pyplot as plt
import os
import h5py
from scipy.signal import argrelmax #Calculate the relative maxima of data.


def plot_peak_freqs(path, axes, Nmaxes = 1):
    # Nmaxes number of maximum frequencies on plot

    stim_freq = []

    peak_freqs = []
    for _ in range(Nmaxes):
        peak_freqs.append([])

    counter = 1
    for file in os.listdir(path):
        if file.split(".")[-1] != "hdf5":
            continue

        #Переписать для файлов с более чем 1 параметрами
        freq = file.split("_")[1]
        freq = float(freq) #Получение значения частоты из названия файла
        
        filepath = path + file
        with h5py.File(filepath, 'r') as h5file:
            
            
            stim_freq.append(freq) # read from h5file !!!!

            fft_amplitudes = h5file["extracellular/electrode_1/firing/processing/ngf/fft_amplitudes"][:]
            fft_freqs = h5file["extracellular/electrode_1/firing/processing/ngf/fft_freqs"][:]
            peak_freq_args = argrelmax(fft_amplitudes, mode="wrap")[0]
            
            

            fft_freqs_peak = fft_freqs[peak_freq_args]
            peak_freq_args = np.argsort( -fft_amplitudes[peak_freq_args] )

           
            for peak_freq_idx in range(Nmaxes):
                if peak_freq_args.size == 0:
                    peak_freqs[peak_freq_idx].append(0)
                else:
                    peak_freqs[peak_freq_idx].append(fft_freqs_peak[peak_freq_args[peak_freq_idx]])

        counter += 1

    stim_freq = np.asarray(stim_freq)
    for idx in range(Nmaxes):
        axes.scatter(stim_freq, peak_freqs[idx], label=str(idx))

    axes.legend()
    #axes.set_xlabel("Stim freq, Hz")
    axes.set_xlabel("Stim strength")
    axes.set_ylabel("Maximal freq of spike rate, Hz")

def plot_mean_CV(path, axes, parameter="mean_cv"):
    stim_freq = []
    vals = []
    counter = 1
    

    
    for file in os.listdir(path):
        if file.split(".")[-1] != "hdf5":
            continue

        print(file)
        freq = file.split("_")[1]
        freq = float(freq)
        
        filepath = path + file
        with h5py.File(filepath, 'r') as h5file:
            stim_freq.append(freq)  # read from h5file !!!!
            # .format(parameter)
            val = h5file["extracellular/electrode_1/firing/processing/ngf/{}".format(parameter)][()]
            vals.append(val)
            counter += 1

    print(len(vals))
    axes.scatter(stim_freq, vals, label=parameter)
    axes.legend()
    
def plot_Freq_Box(path, axes):
    
    CellFreqslist = []
    x_p = []
    for file in sorted(os.listdir(path)):
        if file.split(".")[-1] != "hdf5":
            continue

        print(file)
        freq = file.split("_")[1]
        x_p.append(freq)
        
        filepath = path + file
        with h5py.File(filepath, 'r') as h5file:
            CellFreqs = h5file["extracellular/electrode_1/firing/processing/ngf/Spike_Freq_Set"][:]
            CellFreqslist.append(CellFreqs)
    
    #CellFreqslist = sorted(CellFreqslist, key = lambda)
    
    axes.boxplot(CellFreqslist, labels=x_p)

if __name__ == '__main__':
    path = "./Results/20Hz/"
    fig, axes = plt.subplots(ncols=4, figsize=(10, 5))
    plot_peak_freqs(path, axes[0])
    plot_mean_CV(path, axes[1], parameter="mean_cv")
    plot_mean_CV(path, axes[2], parameter="mean_correlation")
    plot_Freq_Box(path, axes[3])
    
    #fig.savefig(path + "Figures/with_gap.png")

    plt.show()


