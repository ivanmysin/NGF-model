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

        parzen_win = get_window("parzen", 151)
        for celltype in firing_origin.keys():
            if celltype != "ngf": continue # ngf !!!!!!

            firing_process_group_celltype = firing_process_group.create_group(celltype)
            spike_rate_pop = np.zeros( tsim.size, dtype=np.float) #Заготовка для вычисления популяционной частоты разряда!
            spike_rate_cells = []
            mean_cv = 0
            Spike_Freq_Set = []
            
            #Вычисление мгновенной частоты разряда для популяции и для отдельных нейронов,
            #коэффициентов вариации межспайковых интервалов,
            #Распределениz средних частот отдельных нейронов
            counter = 1
            for dsetname, firings_dset in firing_origin[celltype].items(): #Получение данных от конкретных нейронов
                firing_process_group_celltype_number = firing_process_group_celltype.create_group(dsetname) #Создание папок под обработку каждого нейрона
                if firings_dset.size == 0:
                    firing_process_group_celltype_number.create_dataset("cv", data=np.nan)
                    continue
                isi = np.diff(firings_dset[:])
                cv = np.std(isi) / np.mean(isi)
                mean_cv += cv
                firing_process_group_celltype_number.create_dataset("cv", data=cv)
                
                cell_spike_freq = firings_dset.size
                firing_process_group_celltype_number.create_dataset("Cell_Spike_Freq", data = cell_spike_freq)
                Spike_Freq_Set.append(cell_spike_freq)
                
                spike_rate_cell = np.zeros( tsim.size, dtype=np.float)
                spike_rate_cell[ np.floor( firings_dset[:] / delta_t).astype(np.int) ] += 1 #Заполняет единицами элементы
                #массива, номера которых соответствуют единичному временному отрезку delta_t, в течение которого зарегистрировался спайк
                #массив не должен быть типа list, а np.array, например
                spike_rate_cell = np.convolve(parzen_win, spike_rate_cell, mode="same")
                
                firing_process_group_celltype_number.create_dataset("spike_rate", data=spike_rate_cell)
                spike_rate_pop += spike_rate_cell
                spike_rate_cells.append(spike_rate_cell)

                counter += 1

            mean_cv /= counter
            spike_rate_pop /= counter
            
            #Вычисление спектров популяционной мгновенной частоты нейронов
            fft_amplitudes = 2 * np.abs( np.fft.rfft(spike_rate_pop) ) / spike_rate_pop.size
            fft_amplitudes = fft_amplitudes[1:]
            fft_freqs = np.fft.rfftfreq(spike_rate_pop.size, d=0.001*delta_t)
            fft_freqs = fft_freqs[1:]

            fft_amplitudes = fft_amplitudes[fft_freqs < 500]
            fft_freqs = fft_freqs[fft_freqs < 500]
            #Вывод групповых значений
            firing_process_group_celltype.create_dataset("fft_amplitudes", data = fft_amplitudes)
            firing_process_group_celltype.create_dataset("fft_freqs", data = fft_freqs)
            firing_process_group_celltype.create_dataset("mean_cv", data = mean_cv)
            firing_process_group_celltype.create_dataset("pop_spike_rate", data = spike_rate_pop)
            firing_process_group_celltype.create_dataset("Spike_Freq_Set", data = Spike_Freq_Set)
            
            #Вычисление ковариационной матрицы для мгновенных частот разряда отдельных нейронов 
            if len(spike_rate_cells) > 2:
                corrMatrix = np.corrcoef( np.vstack(spike_rate_cells) )
                mean_correlation = np.mean(corrMatrix)
            else:
                mean_correlation = np.nan

            firing_process_group_celltype.create_dataset("mean_correlation", data=mean_correlation)
    return
        
def Sort_To_Folders(path):
    import os
    import shutil
    print(path)
    for file in os.listdir(path): #Названия папок тоже считываются
        file_name = file.split(".")
        if file_name[-1] != "hdf5":
            continue
        file_name = file_name[0].split("_") 
        freq_of_file = file_name[0]
        print("freq_of_file: ", freq_of_file)
        destination = path + "/" + freq_of_file
        try:
            os.mkdir(destination)
        except FileExistsError:
            pass
        shutil.move(file, destination)
        
if __name__ == "__main__":
    # from basic_parameters import get_basic_params
    #
    # basic_params = get_basic_params()
    # filepath = basic_params["file_results"]
    path = "./Results/"
    Sort_To_Folders(path)
    for file in os.listdir(path):
        if file.split(".")[-1] == "pickle":
            continue
        path_to_files = path + "/" + file
        print("path_to_files: ", path_to_files)
        for f in os.listdir(path_to_files):
            processing_and_save(filepath)
    '''
    for file in os.listdir(path):
        if file.split(".")[-1] != "hdf5":
            continue
        filepath = path + file
        print(filepath)
        processing_and_save(filepath)
    '''





