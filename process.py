import numpy as np
import os
from elephant import signal_processing as sigp 
import matplotlib.pyplot as plt
import h5py
import processingLib as plib
from scipy.signal import get_window



processing_param = {

}

def processing_and_save(filepath):
    with h5py.File(filepath, 'a') as h5file:
        tsim = h5file["time"]
        delta_t = tsim[1] - tsim[0]
        duration = tsim[-1]

        firing_group = h5file["extracellular/electrode_1/firing"]

        try:
            firing_process_group = firing_group.create_group("processing")
        except ValueError:
            del h5file["extracellular/electrode_1/firing/processing"]
            firing_process_group = firing_group.create_group("processing")

        firing_origin = firing_group["origin_data"]

        parzen_win = get_window("parzen", 51)
        for celltype in firing_origin.keys():
            if celltype != "ngf": continue # ngf !!!!!!

            firing_process_group_celltype = firing_process_group.create_group(celltype)
            spike_rate_pop = np.zeros( tsim.size, dtype=np.float)
            spike_rate_cells = []
            mean_cv = 0

            counter = 1
            for dsetname, firings_dset in firing_origin[celltype].items():
                firing_process_group_celltype_number = firing_process_group_celltype.create_group(dsetname)
                if firings_dset.size == 0:
                    firing_process_group_celltype_number.create_dataset("cv", data=np.nan)

                    continue
                isi = np.diff(firings_dset[:])
                cv = np.std(isi) / np.mean(isi)
                mean_cv += cv
                firing_process_group_celltype_number.create_dataset("cv", data=cv)


                spike_rate_cell = np.zeros( tsim.size, dtype=np.float)
                spike_rate_cell[ np.floor(delta_t * firings_dset[:]).astype(np.int) ] += 1
                spike_rate_cell = np.convolve(parzen_win, spike_rate_cell, mode="same")
                spike_rate_pop += spike_rate_cell
                spike_rate_cells.append(spike_rate_cell)

                counter += 1

            mean_cv /= counter
            spike_rate_pop /= counter

            fft_amplitudes = 2 * np.abs( np.fft.rfft(spike_rate_pop) ) / spike_rate_pop.size
            fft_amplitudes = fft_amplitudes[1:]
            fft_freqs = np.fft.rfftfreq(spike_rate_pop.size, d=0.001*delta_t)
            fft_freqs = fft_freqs[1:]

            fft_amplitudes = fft_amplitudes[fft_freqs < 500]
            fft_freqs = fft_freqs[fft_freqs < 500]

            firing_process_group_celltype.create_dataset("fft_amplitudes", data = fft_amplitudes)
            firing_process_group_celltype.create_dataset("fft_freqs", data = fft_freqs)
            firing_process_group_celltype.create_dataset("mean_cv", data = mean_cv)

            if len(spike_rate_cells) > 2:
                corrMatrix = np.corrcoef( np.vstack(spike_rate_cells) )
                mean_correlation = np.mean(corrMatrix)
            else:
                mean_correlation = np.nan

            firing_process_group_celltype.create_dataset("mean_correlation", data=mean_correlation)




    return
        
        
if __name__ == "__main__":
    # from basic_parameters import get_basic_params
    #
    # basic_params = get_basic_params()
    # filepath = basic_params["file_results"]
    path = "./Results/"
    for file in os.listdir(path):

        if file.split(".")[-1] != "hdf5":
            continue
        filepath = path + file
        print(filepath)
        processing_and_save(filepath)






