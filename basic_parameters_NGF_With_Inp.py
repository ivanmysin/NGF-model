import numpy as np
import presimulation_lib as prelib
from copy import deepcopy
import mpi4py
import pickle
"""
pyr - CA1 piramidal cells
olm - Oriens-Lacunosum/Moleculare (OLM) Cells
pvbas - parvalbumin positive basket cells
cckbas - cholecystokinin positive basket cells
ivy - Ivy cells
ngf - neurogliaform cells
bis - bistratifered cells
aac - axo-axonal cells
csa - Schaffer collateral-associated cells
ca3_spatial - CA3 pyramidal cells with place cells dynamics
ca3_non_spatial - CA3 pyramidal cells without place cells dynamics
mec - stellate cells of the medial entorhinal cortex
msach - medial septal cholinergic neurons
msteevracells  - medial septal GABAergic neurons (Teevra cells)
"""

def get_basic_params(freq_param, ampl_param):
    basic_params = {
        #"file_results":  "./Results/theta_state.hdf5", # non_theta_state_ripples.hdf5",   #file for saving results
        #"file_params" : "./Results/params.pickle",
        "file_results":  "./Results/" + str(freq_param) + "Hz_" + str(ampl_param/0.016) + "_Mult" + ".hdf5", 
        "file_params" : "./Results/params_" + str(freq_param) + "Hz_" + str(ampl_param/0.016) + "_Mult" + ".pickle",
        #"file_results":  "./Results/theta_state" + str(freq_param) + ".hdf5", # non_theta_state_ripples.hdf5",   #file for saving results
        #"file_params" : "./Results/params" + str(freq_param) + ".pickle",
        "duration" : 5200, # 10 sec simulation time

        "del_start_time" : 200, # time after start for remove

        # number of cells
        "CellNumbers" : {
            "Nngf" : 200, # 130,

            "Nnon_spatial" : 50,
        },
        "CellParameters" : {
            "non_spatial" : {
                "cellclass" : "ArtifitialCell",
                "R" : 0.3,     # gamma R
                "mu" : np.pi, #1.3,    # theta phase (В радианах?)
                "freqs" : freq_param, # theta frequency
                "spike_rate" : 1.0,
                "latency" : 0, #10.0,
                "delta_t" : 0.1, #?
            },
            
            "ngf" : {
                "cellclass" : "ngfcell",
                "iext" : 0.00,
                "iext_std" : 0.00,
            },
        },
        # indexes of neurons for saving soma potentials
        # Сколько нейронов данного типа записывать
        "save_soma_v" : {
            "ngf" : np.arange(200), #НЕ ЗАБЫТЬ!!!
        },

        # parameters of connections
        "connections" : {
            # connection to pyramidal neurons            
            "non_spatial2ngf": {
                "gmax": ampl_param, #0.016 * 300, #
                "gmax_std" : 0, #0.002, # 0.5, #
                
                "Erev": 0,
                "tau_rise": 0.5,
                "tau_decay": 3,

                "prob": 1, # 0.3, 
                
                "delay": 1.5,
                "delay_std" : 0.5,
                

                "sourse_compartment" : "acell",
                "target_compartment" : "dendrite_list",

                #"NMDA" : {
                #    "gNMDAmax" : 0.02, # mS
                #    "gmax_std" : 0.001,
                #    "tcon" : 2.3,   # ms
                #    "tcoff" : 95.0, # ms
                #    "enmda" : 0, 
                #},
            },         
            
           "ngf2ngf": {
                "gmax": 0.75,
                "gmax_std" : 0, #0.7, #Из-за разброса даже выигравшие существование синапсы могут занулиться!
                
                "Erev": -75,
                "tau_rise": 3.1,
                "tau_decay": 42,

                "prob": 1, #17/200,
                
                "delay": 1.2,
                "delay_std" : 0.2,
                

                "sourse_compartment" : "soma",
                "target_compartment" : "dendrite_list",
            },
        }, # end of connetion settings

        # parameters of gap junctions
        "gap_junctions_params" : {           
            "ngf2ngf" : {
                "r" : 1e6,
                "r_std" : 10,
                "prob": 1,
                "compartment1" : "dendrite_list",
                "compartment2" : "dendrite_list",
            },
        },
    }    
    return basic_params

def get_object_params(Nthreads=1, freq_param=0, ampl_param = 1):
    """
    Function return list.
    Each element of the list is dictionary of parameters for one thread of CPU
    """

    basic_params = get_basic_params(freq_param, ampl_param)
    
    OBJECTS_PARAMS = []
    for _ in range(Nthreads):
        thread_param = {
            "neurons" : [],
            "save_soma_v" : None,
            "gids_of_celltypes" : [],
            "synapses_params" : [],
            "gap_junctions" : [],
            "common_params" : {},
        }

        OBJECTS_PARAMS.append(thread_param)

    cell_types_in_model = [] #Содержит и нейроны, и генераторы. Позже входит в состав OBJECTS_PARAMS.
    gids_of_celltypes = {} #Общий список гидов, соответствующий общему списку нейронов до распределения по потокам

    for celltype, numbers in sorted(basic_params["CellNumbers"].items()):

        celltype = celltype[1:] #Название типа (отбрасываем первый символ - N)
        start_idx = len(cell_types_in_model) #Индекс первого элемента данного типа в общем массиве
        cell_types_in_model.extend( [celltype, ] * numbers ) #Заполнение участка общего массива однотипными элементами
        end_idx = len(cell_types_in_model) #Индекс последнего элемента данного типа в общем массиве
        gids_of_celltypes[celltype] = np.arange(start_idx, end_idx) #Гиды, видимо, соответствуют индексам в общем массиве
        #gids_of_celltypes здесь локальная переменная, которая никуда не сохраняется?
        #Амонимия: ниже gids_of_celltypes еще и ключ от словаря, который заполняется без участия локальной 
        #переменной gids_of_celltypes
    #print("cell_types_in_model", cell_types_in_model)
    #print("gids_of_celltypes", gids_of_celltypes["ngf"])
    #print("gids_of_celltypes", gids_of_celltypes["ivy"])
    #print("gids_of_celltypes", gids_of_celltypes["aac"])
    
    Ncells = len(cell_types_in_model) #Генераторы также входят в это число
    
    #Видимо, берет общее число нейронов, растаскивает их по потокам и нумерует от 0 до нейроны_на_поток
    neurons_by_threads = np.tile(np.arange(Nthreads), int(np.ceil(Ncells/Nthreads)) )
    neurons_by_threads = neurons_by_threads[:Ncells]
    #print("neurons_by_threads", neurons_by_threads)
    #Выбор подмножества нейронов для записи МП:
    save_soma_v_idx = np.empty(shape=0, dtype=np.int)#одна на все потоки

    for celltype, list_idx in basic_params["save_soma_v"].items():
        #print("celltype", celltype, "list_idx", list_idx)
        #Видимо, сбрасывает нумерацию элементов подмассива, индексы первых элементов которых отличны от нуля
        #Перенумеровывает подмножества из общего массива.
        indices = [i for i, x in enumerate(cell_types_in_model) if x == celltype]
        #print("indices", indices)
        if len(indices) == 0:
            continue

        indices = np.asarray(indices) #Список индексов, соответствующий общему массиву нейронов
        list_idx = np.asarray(list_idx) #"Обнуленный" список индексов
        #print("indices", indices)
        #print("list_idx", list_idx)
        
        save_soma_v_idx = np.append(save_soma_v_idx, indices[list_idx[list_idx<len(indices)]])
        #print("save_soma_v_idx", save_soma_v_idx)

    for th_idx in range(Nthreads):
        if th_idx == 0:
            #Образец для сборки после симуляции (каждый нейрон помещается на свою начальную позицию и исходя из 
            #нее строятся графики и ведутся обсчеты. Также (!) отфильровываются нейроны, которые должны обсчитываться 
            #в нулевом потоке, т е то же, что делается после else для других потоков.)
            OBJECTS_PARAMS[th_idx]["save_soma_v"] = save_soma_v_idx
            #print("for zero thread: ", save_soma_v_idx)
        else:
            #Сохранение индексов общего массива, которые соответствуют номеру того или иного потока
            save_soma_v_idx_tmp = save_soma_v_idx[neurons_by_threads[save_soma_v_idx] == th_idx]
            OBJECTS_PARAMS[th_idx]["save_soma_v"] = save_soma_v_idx_tmp # !!!!!!!
            #print("neurons_by_threads", neurons_by_threads)
            #print("save_soma_v_idx_tmp", save_soma_v_idx_tmp)
        #Гиды (индивидуален для каждого нейрона) разбрасываются по потокам по тому же принципу, что и индексы.
        OBJECTS_PARAMS[th_idx]["gids_of_celltypes"] = np.arange(len(cell_types_in_model))[neurons_by_threads == th_idx]
        #print("gids_of_celltypes", gids_of_celltypes)
        #print("OBJECTS_PARAMS[th_idx][gids_of_celltypes] ", OBJECTS_PARAMS[th_idx]["gids_of_celltypes"])
    OBJECTS_PARAMS[0]["cell_types_in_model"] = cell_types_in_model

    #Задание конкретных параметров для разных типов генераторов    
    for celltypename, cellparam in basic_params["CellParameters"].items():
        #print("celltypename:", celltypename, "cellparam:", cellparam)
        if cellparam["cellclass"] == "ArtifitialCell": #Здесь задаются параметры привязки для всего класса генераторов.
            Rgen = cellparam["R"]
            kappa, I0 = prelib.r2kappa(Rgen)
            cellparam["kappa"] = kappa #этих ключей нет в исходном словаре параметров
            cellparam["I0"] = I0     
        else:
            continue
    ArtCell_idx = 1 
    for cell_idx, celltype in enumerate(cell_types_in_model):
        #print("cell_idx:", cell_idx, "celltype:", celltype)
        cell_param = basic_params["CellParameters"][celltype] #"Размножение" словарей параметров (полные клоны базового словаря)
        #Задание базовых свойств каждого нейрона (включая генераторы)
        neuron = {#Создается отдельный словарь для каждого нейрона.
            "celltype" : celltype, 
            "cellclass" : cell_param["cellclass"],
            "cellparams" : {},
            "gid" : cell_idx, #Гид задается дважды?
        }
        neuron["cellparams"] = deepcopy(cell_param)
        #Модификация параметров для каждого гененратора            
        if cell_param["cellclass"] == "ArtifitialCell":
            #Вариация силы фазовой привязки:
            Rgen = ArtCell_idx/10*0.2
            if Rgen >= 1: Rgen = 0.99
            kappa, I0 = prelib.r2kappa(Rgen)
            neuron["cellparams"]["kappa"] = kappa 
            neuron["cellparams"]["I0"] = I0  
            #Вариация частоты разряда:
            neuron["cellparams"]["spike_rate"] = 5 + 2 * (ArtCell_idx/10)
            #Генератор случайных чисел: нет параметра!
            ArtCell_idx += 1
            #print("kappa for ArtCell ", cell_idx, ":", neuron["cellparams"]["kappa"]) 
        else:
            if cell_param["iext"] > 0:
                neuron["cellparams"]["iext"] = np.random.lognormal( np.log(cell_param["iext"]), cell_param["iext_std"]   )
            else:
                neuron["cellparams"]["iext"] = np.random.normal(cell_param["iext"], cell_param["iext_std"] )
            

        th_idx = int(neurons_by_threads[cell_idx])
        OBJECTS_PARAMS[th_idx]["neurons"].append(neuron)
        
    #Задание параметров конкретных синапсов между конкретными нейронами
    for presynaptic_cell_idx, pre_celltype in enumerate(cell_types_in_model):
        for postsynaptic_cell_idx, post_celltype in enumerate(cell_types_in_model):
            dist_normalizer = 0
            conn_name = pre_celltype + "2" + post_celltype


            if presynaptic_cell_idx == postsynaptic_cell_idx: continue

            try:
                conn_name = pre_celltype + "2" + post_celltype
                conn_data = basic_params["connections"][conn_name]

            except AttributeError:
                continue
            except KeyError:
                continue
            
            #Для масштабирования при отладке ("вероятность" может стать больше 1, если мы хотим уменьшить 
            #количество нейронов, но сохранить плотность связей на них. Можно заменить масштабированием весов):
            number_connections = int( np.floor(conn_data["prob"]) )
            #print("number_connections", number_connections)
            if (np.random.rand() < (conn_data["prob"] - number_connections) ):
                number_connections += 1 #либо 0 либо 1, если вероятность от 0 до 1 и больше единицы, ели вероятность больше 1.
                #Для случая болльше 1 ниже идет пробежка по циклу.

            gmax = deepcopy(conn_data["gmax"])
            #Преобразование "условной" дальности в индексах в реальную дальность в координатах между пирамидами
            #и искусственными генераторами 

            #Задание задержек и прочих параметров:
            for _ in range(number_connections): #

                delay = np.random.lognormal(mean=np.log(conn_data["delay"]), sigma=conn_data["delay_std"]) 
                if delay <= 0.5:
                    delay = 0.5
                
                gmax_syn =  np.random.normal(loc=gmax, scale=conn_data["gmax_std"])
                #np.random.lognormal(mean=np.log(gmax), sigma=conn_data["gmax_std"])

                if gmax_syn < 0: #0.000001:
                    continue
                    
                connection = { #Здесь только соединения нейрон с нейроном, не генератор на нейрон?
                    "pre_gid" : presynaptic_cell_idx,
                    "post_gid" : postsynaptic_cell_idx,
                    
                    "gmax" : gmax_syn,
                    "Erev" : conn_data["Erev"],
                    "tau_rise" : conn_data["tau_rise"],
                    "tau_decay" : conn_data["tau_decay"],
                    "delay" : delay,
                    
                    "sourse_compartment" : conn_data["sourse_compartment"], #В случает генератора 'acell'
                    "target_compartment" : conn_data["target_compartment"],
                }
                #if pre_celltype == 'non_spatial':
                    #print(pre_celltype + "2" + post_celltype, presynaptic_cell_idx, "to", postsynaptic_cell_idx, "gmax", connection["gmax"])
                #НМДА пока отключены
                try: 
                    gmax_nmda = conn_data["NMDA"]["gNMDAmax"]
                    gmax_nmda =  np.random.lognormal(mean=np.log(gmax_nmda), sigma=conn_data["NMDA"]["gmax_std"])
                    
                    connection["NMDA"] = {
                        "gNMDAmax" : gmax_nmda,
                        "tcon" : conn_data["NMDA"]["tcon"],   
                        "tcoff" : conn_data["NMDA"]["tcoff"], 
                        "enmda" : conn_data["NMDA"]["enmda"], 
                    }

                except KeyError:
                    pass

                connection["gmax"] *= 0.001  # recalulate nS to micromhos
                connection["delay"] += 1.5  # add delay on spike generation
                try:
                    connection["NMDA"]["delay"] += 1.5
                    connection["NMDA"]["gNMDAmax"] *= 0.001
                except KeyError:
                    pass
                # synapses.append(connection)

                th_idx = int(neurons_by_threads[postsynaptic_cell_idx])
                OBJECTS_PARAMS[th_idx]["synapses_params"].append(connection)

    #counter = 0
    #count_of_doubles = 0
    gap_juncs = []
    sgid_gap = 0 #Для щелевых контактов нужны свои гиды.
    #Создание щелевых контактов между нейронами:
    for cell1_idx, celltype1 in enumerate(cell_types_in_model):
        for cell2_idx, celltype2 in enumerate(cell_types_in_model):

            if cell1_idx == cell2_idx: continue

            try:
                conn_name = celltype1 + "2" + celltype2
                conn_data = basic_params["gap_junctions_params"][conn_name]
            #Поскольку лишь по некоторым ключам есть параметры щелевых контактов, при переборе 
            #всего массива контактов будут вылетать исключения, которые можно игнорировать и 
            #продолжать перебор.
            except AttributeError:
                continue
            except KeyError:
                continue
            #print(celltype1, cell1_idx, "2", celltype2, cell2_idx)
            #counter +=1
            if (np.random.rand() > conn_data["prob"]): continue
           
            gap = {
                "gid1" : cell1_idx,
                "gid2" : cell2_idx,
                "r" : np.random.normal(conn_data["r"], conn_data["r_std"], 1),
                
                "compartment1" : conn_data["compartment1"],
                "compartment2" : conn_data["compartment2"],
                
                "sgid_gap" : sgid_gap,
            }
            
            gap_juncs.append(gap)
            sgid_gap += 2

            th_idx = int(neurons_by_threads[cell1_idx])
            OBJECTS_PARAMS[th_idx]["gap_junctions"].append(gap)
            
            th_idx2 = int(neurons_by_threads[cell2_idx])
            if th_idx2 != th_idx:
                OBJECTS_PARAMS[th_idx2]["gap_junctions"].append(gap)
                #count_of_doubles += 1
    #print("gap_junctions from vector", len(OBJECTS_PARAMS[th_idx]["gap_junctions"]))
    #print("gap_junctions from counter", counter)
    #print("gaps count_of_doubles", count_of_doubles)
    
    for th_idx in range(Nthreads):
        OBJECTS_PARAMS[th_idx]["duration"] = deepcopy(basic_params["duration"])
        OBJECTS_PARAMS[th_idx]["file_results"] = deepcopy(basic_params["file_results"])
        OBJECTS_PARAMS[th_idx]["del_start_time"] = deepcopy(basic_params["del_start_time"])

    if not(basic_params["file_params"] is None):
        with open(basic_params["file_params"], 'wb') as file_param:
            pickle.dump(OBJECTS_PARAMS, file_param)
        
        
    return OBJECTS_PARAMS


if __name__ == "__main__":
    # This test code. You can run and see format of list and dictionaries genereated by functions in this file
    Nthreads = 4
    objc_p_list = []
    ArtFreqList = np.arange(5, 105, 5)
    #ArtFreqList = [100,]
    ArtAmplList = 0.016 * np.arange(100, 1000, 100)
    print(ArtAmplList)
    for ampl_param in ArtAmplList:
        for freq_param in ArtFreqList:
            print("ampl_param", ampl_param, "freq_param", freq_param)
            objc_p_list.append(get_object_params(Nthreads=Nthreads, freq_param=freq_param, ampl_param=ampl_param))
   