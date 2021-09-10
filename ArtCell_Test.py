"""

"""
GenNum = 50
Freq = 5
SimulTime = 1000
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

import neuron as n
import presimulation_lib as prelib
import math
import matplotlib.pyplot as ppl
n.h.load_file("stdrun.hoc")
n.h.load_file("stdgui.hoc")
n.h.load_file("import3d.hoc")

n.load_mechanisms("./mods/") #(r"C:\Users\Povarov\Mechs\ArtCells")
pc = n.h.ParallelContext()
s = n.h.Section(name='soma')
s.insert('pas')
SynList = []

for i in range(GenNum): #Для каждого генератора свой синапс
    syn = n.h.Exp2Syn(0.5, sec = s)
    SynList.append(syn)
    #По Bezaire et al 2016
    SynList[i].tau1 = 2
    SynList[i].tau2 = 6.3
    SynList[i].e = 0
    SynList[i].g = 0.0035 #Настройка этого параметра ни на что не влияет
    #G_max - как выставить? В mod файле есть просто g, измеряется в мкС

Arts = n.h.List()
GenTimes = n.h.List()
ConnList = n.h.List()
for i in range(GenNum):
    GenTimes.append(n.h.Vector())
    pre_R = i/10*0.2 #math.floor(i/10) * 0.2
    if pre_R >= 1: pre_R = 0.99
    theta_kappa, theta_i0 = prelib.r2kappa(pre_R) #Положим разную степень фазовой привязки у разных генераторов
    Arts.append(n.h.ArtificialRhytmicCell(0, 0))
    Arts[i].mu = math.pi #фаза привязки
    Arts[i].latency = 0 #Регулирует рефрактерный период
    Arts[i].kappa = theta_kappa #* (Freq/10)
    #print("theta_kappa: ", theta_kappa)
    Arts[i].I0 = theta_i0  
    Arts[i].spike_rate = 5 + 2 * (i/10) #* (Freq/5)*0.3*(Freq != 5) #Если хотим получить одинаковое кол-во спайков на период любой длительности, 
    #частоту разряда делаем пропорциональной несущей частоте.
    Arts[i].freqs = Freq #регулирует частоту несущей волны, к фазам которой осуществляется привязка 
    Arts[i].myseed = i*2
    pc.set_gid2node(i, pc.id())
    ProtoCon = n.h.NetCon(Arts[i], None)
    ProtoCon.record(GenTimes[i])
    pc.cell(i, ProtoCon)
    conn = pc.gid_connect(i, SynList[i])
    conn.weight[0] = 0.1 #+ Freq/50
    ConnList.append(conn)

t = n.h.Vector()
t.record(n.h._ref_t)
sv = n.h.Vector()
sv.record(s(0.5)._ref_v)
sp = n.h.Vector()
#conn.record(sp)

n.h.tstop = SimulTime
n.h.celsius = 37

pc.set_maxstep(5)
n.h.finitialize()
pc.psolve(SimulTime) #замена функции run()
pc.barrier()

#n.h.run()
import numpy as np
########## Растерограмма разрядов всех генераторов + Мембранный ответ ##########
fig, axes = ppl.subplots(nrows=2, sharex=True)  #ppl.figure()
for i in range(GenNum):
    times = np.asarray(GenTimes[i])
    points = np.ones(len(times)) * (i + 1) 
    axes[0].scatter(times, points, color='k', s=0.5)          

axes[1].plot(t, sv,  color="black")
axes[1].set_yticks(np.arange(-80, 0, 5))
#axes[1].set_ylim(-80, 40)
ppl.show()

########### Активность полной популяции за маленькие интервалы ##################
fig, axes_1 = ppl.subplots(nrows=1, sharex=True) 
bin_size = 10 #ms
PopFreq = np.zeros(int(SimulTime/bin_size)) 
border = 0
for bin_num in range(len(PopFreq)):
    border += bin_size
    for TVec in GenTimes:
        for t in np.asarray(TVec):
            if t <= border and t > (border - bin_size):
                PopFreq[bin_num] += 1
PopFreq /= (SimulTime/1000) 
BinTimes = np.arange(int(SimulTime/bin_size)) * bin_size
axes_1.plot(BinTimes, PopFreq,  color="black") 

############# Дополнительные параметры ####################
AP_sum_min = len(np.asarray(GenTimes[0]))
AP_sum_max = len(np.asarray(GenTimes[GenNum - 1]))
print("AP_sum_min: ", AP_sum_min, " AP_sum_max: ", AP_sum_max)
MP_mean = np.mean(np.array(sv))
print("MP_mean: ", MP_mean) #Возможно, не лучший показатель, хотя, может помочь в оценке тонического компонента.
AmpVec = []
T = int(1000/Freq)
sv = np.array(sv)
#WARNUNG!!! 40 значений на 1 мс
for i in range(Freq): #Вычисление амплитуды колебаний
    amp_per_T = abs(np.min(sv[(i*T*40):((i+1)*T*40)])) - abs(np.max(sv[(i*T*40):((i+1)*T*40)]))
    #print("pare_", i, "-", i+1, ": ", i*T, (i+1)*T, "min:", min(sv[(i*T*40):((i+1)*T*40)]), "max:", max(sv[(i*T*40):((i+1)*T*40)]))
    AmpVec.append(amp_per_T)
Amp = np.mean(np.array(AmpVec))
print("Amp: ", Amp)

######## Спектры популяционной частоты разряда генераторов и колебаний МП ############
fig_2, axes_2 = ppl.subplots(nrows=1, sharex=True)
delta_t = 1000/40000
fft_ampl_MP = 2 * np.abs( np.fft.rfft(sv) ) / sv.size
fft_ampl_MP = fft_ampl_MP[1:] 
fft_freqs_MP = np.fft.rfftfreq(sv.size, d=0.001*delta_t) #Домножаем для получения Гц
fft_freqs_MP = fft_freqs_MP[1:]
axes_2.plot(fft_freqs_MP, fft_ampl_MP, color = 'black')
fft_ampl_PopFreq = 2 * np.abs( np.fft.rfft(PopFreq) ) / PopFreq.size
fft_ampl_PopFreq = fft_ampl_PopFreq[1:] 
fft_freqs_PopFreq = np.fft.rfftfreq(PopFreq.size, d=0.001*bin_size) #В качестве длительности сэмпла берем длительность бина
fft_freqs_PopFreq = fft_freqs_PopFreq[1:]
axes_2.plot(fft_freqs_PopFreq, fft_ampl_PopFreq, color = 'red')
#axes_2.set_xticks(np.arange(0, 2500, 100))
axes_2.set_xlim(0, 200)
ppl.show()