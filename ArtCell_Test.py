"""
Для того, чтобы ввести несколько генераторов не в фазе нужно разнообразить
myseed
"""
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

import neuron as n
import presimulation_lib as prelib
import math
n.h.load_file("stdrun.hoc")
n.h.load_file("stdgui.hoc")
n.h.load_file("import3d.hoc")
#import presimulation_lib as presim

n.load_mechanisms("./mods/") #(r"C:\Users\Povarov\Mechs\ArtCells")
pc = n.h.ParallelContext()
s = n.h.Section(name='soma')
s.insert('pas')
syn = n.h.Exp2Syn(0.5, sec = s)

theta_kappa, theta_i0 = prelib.r2kappa(0.2)
Arts = n.h.List()
GenTimes = n.h.List()
GenTimes.append(n.h.Vector())
GenTimes.append(n.h.Vector())
for i in range(2):
    Arts.append(n.h.ArtificialRhytmicCell(0, 0))
    Arts[i].mu = math.pi #фаза привязки
    Arts[i].latency = 0 #5 + i*10 #Видимо, регулирует рефрактерный период
    Arts[i].kappa = theta_kappa
    Arts[i].I0 = theta_i0
    Arts[i].spike_rate = 20
    Arts[i].freqs = 1 #регулирует частоту несущей волны, к фазам которой осуществляется привязка 
    Arts[i].myseed = i*2
    pc.set_gid2node(i, pc.id())
    ProtoCon = n.h.NetCon(Arts[i], None)
    ProtoCon.record(GenTimes[i])
    pc.cell(i, ProtoCon)
    conn = pc.gid_connect(i, syn)
    conn.weight[0] = 0.1
"""
import presimulation_lib as prelib
art = n.h.ArtificialRhytmicCell(0, 0)
theta_kappa, theta_i0 = prelib.r2kappa(0.9)
art.latency = 5
art.kappa = theta_kappa
art.I0 = theta_i0
art.spike_rate = 20
art.freqs = 5
art.myseed = 0
"""

t = n.h.Vector()
t.record(n.h._ref_t)
sv = n.h.Vector()
sv.record(s(0.5)._ref_v)
sp = n.h.Vector()
#conn.record(sp)

n.h.tstop = 10000
n.h.celsius = 37

pc.set_maxstep(5)
n.h.finitialize()
pc.psolve(1000) #замена функции run()
pc.barrier()

#n.h.run()
import numpy as np
SpikeTimes = np.array(GenTimes[1])
times = np.array(t)
print('Spike_times ', SpikeTimes)
print('size of SpikeTimes   ', np.shape(SpikeTimes))
print(np.shape(SpikeTimes)[0])
tmp = np.zeros((1, np.shape(times)[0])) 
############
'''
a = 0  
for i in range(np.shape(times)[0]):
    if times[i] > SpikeTimes[a]:
        print('AAAAAAAAAAAA')
        tmp[i] = 50
        a = a + 1  
        if a >= np.shape(SpikeTimes):
            break
#################
print('Size of tmp ', np.shape(tmp))
print('size of times ', np.shape(times))
print('delta_11 ', times[1] - times[0], 'delta_t2 ', times[2]-times[1])
#T = np.arange(0, 10000, 0.2)
'''
import matplotlib.pyplot as ppl
ppl.plot(t, sv, color="black")
ppl.ylim(-80, -40)
ppl.show()
#ppl.scatter(times, tmp)
#ppl.show()

