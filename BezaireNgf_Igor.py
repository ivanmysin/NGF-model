# -*- coding: utf-8 -*-
"""
Created on Fri Jan  8 18:54:01 2021

@author: Povarov
"""

import numpy as np
import neuron as n
#import math
        
n.h.load_file("stdrun.hoc")
n.h.load_file("stdgui.hoc")
n.h.load_file("import3d.hoc")

a= n.h.SectionList()

#dend = [n.h.Section(name='dend[%d]' % i) for i in range(3)]

#n.load_mechanisms(r"C:\Users\Povarov\Mechs\NGF_Bezaire")
#n.load_mechanisms(r"C:\Users\Povarov\Mechs\Hoc_Mod")
#n.load_mechanisms(r"C:\Users\Povarov\Mechs\MyExp2Syn")
#n.load_mechanisms(r"C:\Users\Povarov\Mechs\gap")

class BezaireNgf:
    def __init__(self, gid):
        #self.cv = n.h.CVode()
        #self.cv.active(1)
        #self.cv.atol(0.005)
        
        self.gid = gid
        self.SegLen = 10
        self.AllSections = n.h.List()
        self.SynList = n.h.List()
        self.SynList_inh = n.h.List()
        self.SynList_ex = n.h.List()
        
        self._MorphConstruct()
        self._MechInit()
        self._InsertMechs()
        self._SynConstruct()
        #n.h.topology()
        #n.h.Shape()
         
    def _MorphConstruct(self):
        self.soma = n.h.Section(name='soma', cell = self)
        n.h.pt3dclear(sec = self.soma)
        xvec = n.h.Vector([0, 0]); yvec = n.h.Vector([0, 20])
        zvec = n.h.Vector([0, 0]); dvec = n.h.Vector([10, 10])
        n.h.pt3dadd(xvec, yvec, zvec, dvec, sec = self.soma)
        self.AllSections.append(self.soma)
        
        #Apical right:
        self.dend_1 = n.h.Section(name='dend_1', cell = self)
        n.h.pt3dclear(sec = self.dend_1)
        xvec = n.h.Vector([0, 38]); yvec = n.h.Vector([20, 112])
        zvec = n.h.Vector([0, 0]); dvec = n.h.Vector([4, 4])
        n.h.pt3dadd(xvec, yvec, zvec, dvec, sec = self.dend_1)
        self.AllSections.append(self.dend_1)
        self.dend_1.connect(self.soma(1), 0)
        
        self.dend_2 = n.h.Section(name='dend_2', cell = self)
        n.h.pt3dclear(sec = self.dend_2)
        xvec = n.h.Vector([38, 77]); yvec = n.h.Vector([112, 204])
        zvec = n.h.Vector([0, 0]); dvec = n.h.Vector([3, 3])
        n.h.pt3dadd(xvec, yvec, zvec, dvec, sec = self.dend_2)
        self.AllSections.append(self.dend_2)
        self.dend_2.connect(self.dend_1(1), 0)
        
        self.dend_3 = n.h.Section(name='dend_3', cell = self)
        n.h.pt3dclear(sec = self.dend_3)
        xvec = n.h.Vector([77, 155]); yvec = n.h.Vector([204, 388])
        zvec = n.h.Vector([0, 0]); dvec = n.h.Vector([2, 2])
        n.h.pt3dadd(xvec, yvec, zvec, dvec, sec = self.dend_3)
        self.AllSections.append(self.dend_3)
        self.dend_3.connect(self.dend_2(1), 0)
        
        self.dend_4 = n.h.Section(name='dend_4', cell = self)
        n.h.pt3dclear(sec = self.dend_4)
        xvec = n.h.Vector([155, 194]); yvec = n.h.Vector([388, 480])
        zvec = n.h.Vector([0, 0]); dvec = n.h.Vector([1, 1])
        n.h.pt3dadd(xvec, yvec, zvec, dvec, sec = self.dend_4)
        self.AllSections.append(self.dend_4)
        self.dend_4.connect(self.dend_3(1), 0)
        
        self.dend_5 = n.h.Section(name='dend_5', cell = self)
        n.h.pt3dclear(sec = self.dend_5)
        xvec = n.h.Vector([194, 233]); yvec = n.h.Vector([480, 572])
        zvec = n.h.Vector([0, 0]); dvec = n.h.Vector([1, 1])
        n.h.pt3dadd(xvec, yvec, zvec, dvec, sec = self.dend_5)
        self.AllSections.append(self.dend_5)
        self.dend_5.connect(self.dend_4(1), 0)
        
        #Apical left:
        self.dend_6 = n.h.Section(name='dend_6', cell = self)
        n.h.pt3dclear(sec = self.dend_6)
        xvec = n.h.Vector([0, -38]); yvec = n.h.Vector([20, 112])
        zvec = n.h.Vector([0, 0]); dvec = n.h.Vector([4, 4])
        n.h.pt3dadd(xvec, yvec, zvec, dvec, sec = self.dend_6)
        self.AllSections.append(self.dend_6)
        self.dend_6.connect(self.soma(1), 0)
        
        self.dend_7 = n.h.Section(name='dend_7', cell = self)
        n.h.pt3dclear(sec = self.dend_7)
        xvec = n.h.Vector([-38, -77]); yvec = n.h.Vector([112, 204])
        zvec = n.h.Vector([0, 0]); dvec = n.h.Vector([3, 3])
        n.h.pt3dadd(xvec, yvec, zvec, dvec, sec = self.dend_7)
        self.AllSections.append(self.dend_7)
        self.dend_7.connect(self.dend_6(1), 0)
        
        self.dend_8 = n.h.Section(name='dend_8', cell = self)
        n.h.pt3dclear(sec = self.dend_8)
        xvec = n.h.Vector([-77, -155]); yvec = n.h.Vector([204, 388])
        zvec = n.h.Vector([0, 0]); dvec = n.h.Vector([2, 2])
        n.h.pt3dadd(xvec, yvec, zvec, dvec, sec = self.dend_8)
        self.AllSections.append(self.dend_8)
        self.dend_8.connect(self.dend_7(1), 0)
        
        self.dend_9 = n.h.Section(name='dend_9', cell = self)
        n.h.pt3dclear(sec = self.dend_9)
        xvec = n.h.Vector([-155, -194]); yvec = n.h.Vector([388, 480])
        zvec = n.h.Vector([0, 0]); dvec = n.h.Vector([1.5, 1.5])
        n.h.pt3dadd(xvec, yvec, zvec, dvec, sec = self.dend_9)
        self.AllSections.append(self.dend_9)
        self.dend_9.connect(self.dend_8(1), 0)
        
        self.dend_10 = n.h.Section(name='dend_10', cell = self)
        n.h.pt3dclear(sec = self.dend_10)
        xvec = n.h.Vector([-194, -233]); yvec = n.h.Vector([480, 572])
        zvec = n.h.Vector([0, 0]); dvec = n.h.Vector([1, 1])
        n.h.pt3dadd(xvec, yvec, zvec, dvec, sec = self.dend_10)
        self.AllSections.append(self.dend_10)
        self.dend_10.connect(self.dend_9(1), 0)
        
        #Basal left:
        self.dend_11 = n.h.Section(name='dend_11', cell = self)
        n.h.pt3dclear(sec = self.dend_11)
        xvec = n.h.Vector([0, -38]); yvec = n.h.Vector([0, -92])
        zvec = n.h.Vector([0, 0]); dvec = n.h.Vector([2, 2])
        n.h.pt3dadd(xvec, yvec, zvec, dvec, sec = self.dend_11)
        self.AllSections.append(self.dend_11)
        self.dend_11.connect(self.soma(0), 0)
        
        self.dend_12 = n.h.Section(name='dend_12', cell = self)
        n.h.pt3dclear(sec = self.dend_12)
        xvec = n.h.Vector([-38, -77]); yvec = n.h.Vector([-92, -184])
        zvec = n.h.Vector([0, 0]); dvec = n.h.Vector([1.5, 1.5])
        n.h.pt3dadd(xvec, yvec, zvec, dvec, sec = self.dend_12)
        self.AllSections.append(self.dend_12)
        self.dend_12.connect(self.dend_11(1), 0)
        
        self.dend_13 = n.h.Section(name='dend_13', cell = self)
        n.h.pt3dclear(sec = self.dend_13)
        xvec = n.h.Vector([-77, -116]); yvec = n.h.Vector([-184, -276])
        zvec = n.h.Vector([0, 0]); dvec = n.h.Vector([1, 1])
        n.h.pt3dadd(xvec, yvec, zvec, dvec, sec = self.dend_13)
        self.AllSections.append(self.dend_13)
        self.dend_13.connect(self.dend_12(1), 0)
        
        #Basal right:
        self.dend_14 = n.h.Section(name='dend_14', cell = self)
        n.h.pt3dclear(sec = self.dend_14)
        xvec = n.h.Vector([0, 38]); yvec = n.h.Vector([0, -92])
        zvec = n.h.Vector([0, 0]); dvec = n.h.Vector([2, 2])
        n.h.pt3dadd(xvec, yvec, zvec, dvec, sec = self.dend_14)
        self.AllSections.append(self.dend_14)
        self.dend_14.connect(self.soma(0), 0)
        
        self.dend_15 = n.h.Section(name='dend_15', cell = self)
        n.h.pt3dclear(sec = self.dend_15)
        xvec = n.h.Vector([38, 77]); yvec = n.h.Vector([-92, -184])
        zvec = n.h.Vector([0, 0]); dvec = n.h.Vector([1.5, 1.5])
        n.h.pt3dadd(xvec, yvec, zvec, dvec, sec = self.dend_15)
        self.AllSections.append(self.dend_15)
        self.dend_15.connect(self.dend_14(1), 0)
        
        self.dend_16 = n.h.Section(name='dend_16', cell = self)
        n.h.pt3dclear(sec = self.dend_16)
        xvec = n.h.Vector([77, 116]); yvec = n.h.Vector([-184, -276])
        zvec = n.h.Vector([0, 0]); dvec = n.h.Vector([1, 1])
        n.h.pt3dadd(xvec, yvec, zvec, dvec, sec = self.dend_16)
        self.AllSections.append(self.dend_16)
        self.dend_16.connect(self.dend_15(1), 0)
        
        for sec in range(len(self.AllSections)):
            self.AllSections[sec].nseg = round(n.h.distance(self.AllSections[sec](0), 
                                         self.AllSections[sec](1)) / self.SegLen)
        
    def _MechInit(self):
        import math
        Vrest=-65
        self.mechinit = {'celsius':34.0, 
                    'RmDend':5555, 'RmSoma':5555,
                    'CmSoma':1.4, 'CmDend':1.4, 
                    'RaDend':100, 'RaSoma':100,
                    'ca_outside':2, 'ca_inside':5.e-6, 'catau':10, #Ca conc in mM
                    'ekval':-90, 'enaval':55, 'eHCNval':-30, 'ecaval':130, 'eleakval':Vrest,            
                	'gKdr':0.013, 'kdrf':0.15, 'kdrfdend':0.002, #max ion channel conductances in mho/cm2
                    'gKvA_1':0.00015, 'gKvA_2':5.2e-06, 
                    'gCavN_1':0.0008, 'gCavN_2':0.0006, 'gCavL_1':0.005, 'gCavL_2':0.056, 
                    'gKvCaB_1':0.0000002, 'gKvCaB_2':1.02e-06, 'gKCaS_1':0.000002, 'gKCaS_2':4.5e-07,
                    'gNav':0.15, 'gNasoma':3.8, 'gNadend':0.26, 
                    'gleak':8.47e-05, 'gHCN':0.00002,
                    'myRa':14, 'CM':1.8,
                    'offset1':20, 'offset2':20, 'offset3':13.6, 'offset4':13.35,
                    'slope1':0.68, 'slope2':0.57, 'slope3':1.186, 'slope4':49.5,
                    'offset5':9, 'offset6':9, 'slope5':0.07, 'slope6':0.26}
        self.mechinit['ecaval'] = 8.314 * (273.15 + self.mechinit.get('celsius')) \
                                / (2 * 9.649e4) * math.log(self.mechinit.get('ca_outside') \
                                / self.mechinit.get('ca_inside')) * 1000 #about 170, otherwise set to 130
        #self.mechinit['ecaval'] = 130                           
    
    def _InsertMechs(self):
        for i in range(17):
            
            self.AllSections[i].insert('iconc_Ca')
            self.AllSections[i].catau_iconc_Ca = 10
            self.AllSections[i].caiinf_iconc_Ca = 5.e-6
            #self.AllSections[i].cao_iconc_Ca = 2 #ПРОВЕРИТЬ (откуда Z?)!
            self.AllSections[i].insert('ch_KvAngf')
            self.AllSections[i].gmax_ch_KvAngf = self.mechinit.get('gKvA_2') #gKvA
            self.AllSections[i].insert('ch_CavN')#HAV-N- Ca channel
            self.AllSections[i].gmax_ch_CavN = self.mechinit.get('gCavN_2') #gCavN
            self.AllSections[i].insert('ch_CavL')
            self.AllSections[i].gmax_ch_CavL = self.mechinit.get('gCavL_2') #gCavL
            self.AllSections[i].insert('ch_KCaS')
            self.AllSections[i].gmax_ch_KCaS = self.mechinit.get('gKCaS_2') #gKCaS
            self.AllSections[i].insert('ch_KvCaB')
            self.AllSections[i].gmax_ch_KvCaB = self.mechinit.get('gKvCaB_2') #gKvCaB
            
            self.AllSections[i].cm = self.mechinit.get('CM')
            self.AllSections[i].Ra = self.mechinit.get('myRa')
            #self.AllSections[i].ek = self.mechinit.get('ekval')
            self.AllSections[i].eca = self.mechinit.get('ecaval')
            #self.SectionList[i].e_ch_leak = -60 
            self.AllSections[i].insert('ch_leak')
            self.AllSections[i].gmax_ch_leak = self.mechinit.get('gleak')
            self.AllSections[i].e_ch_leak = self.mechinit.get('eleakval')
            #Натриевая и калиевая проводимость в соме и дендритах различны
            #Зависимость от расстояния пока не задана, пока все равномерно
            self.AllSections[i].insert('ch_Navngf')
            if i == 0: #soma
                self.AllSections[i].gmax_ch_Navngf = self.mechinit.get('gNasoma') #0.10*gna_scale
            else: #dend
                self.AllSections[i].gmax_ch_Navngf = self.mechinit.get('gNadend')
            #В mod файле нет параметров offset и slope
            #self.AllSections[i].offset1_ch_Navngf = self.mechinit.get('offset1')
            #self.AllSections[i].offset2_ch_Navngf = self.mechinit.get('offset2')
            #self.AllSections[i].offset3_ch_Navngf = self.mechinit.get('offset3')
            #self.AllSections[i].offset4_ch_Navngf = self.mechinit.get('offset4')
            #self.AllSections[i].slope1_ch_Navngf = self.mechinit.get('slope1')
            #self.AllSections[i].slope2_ch_Navngf = self.mechinit.get('slope2')
            #self.AllSections[i].slope3_ch_Navngf = self.mechinit.get('slope3')
            #self.AllSections[i].slope4_ch_Navngf = self.mechinit.get('slope4')
            self.AllSections[i].ena = 55
            self.AllSections[i].insert('ch_Kdrfastngf')
            if i == 0:
                self.AllSections[i].gmax_ch_Kdrfastngf = self.mechinit.get('kdrf') 
            else:
                self.AllSections[i].gmax_ch_Kdrfastngf = self.mechinit.get('kdrfdend')
            self.AllSections[i].offset5_ch_Kdrfastngf = self.mechinit.get('offset5') 
            self.AllSections[i].offset6_ch_Kdrfastngf = self.mechinit.get('offset6') 
            self.AllSections[i].slope6_ch_Kdrfastngf = self.mechinit.get('slope6') 
            self.AllSections[i].slope5_ch_Kdrfastngf = self.mechinit.get('slope5')	
            #insert ch_Kdrslow
            #gmax_ch_Kdrslow=kdrf
    
    def _SynConstruct(self):
        #Inhibitory from other Ngf cells
        for i in range(16): #последняя цифра не включается
            syn_ = n.h.MyExp2Syn(0.5, sec=self.AllSections[i + 1])
            self.SynList_inh.append(syn_)
            self.SynList_inh[i].tau1 = 2
            self.SynList_inh[i].tau2 = 10
            self.SynList_inh[i].e = -80
        
        #Excitatory inputs from EC3
        for i in range(16): #последняя цифра не включается
            syn_ = n.h.MyExp2Syn(0.5, sec=self.AllSections[i + 1])
            self.SynList_ex.append(syn_)
            self.SynList_ex[i].tau1 = 2
            self.SynList_ex[i].tau2 = 10
            self.SynList_ex[i].e = 40
             
###################################################################        
#pc = n.h.ParallelContext()		
#a=BezaireNgf(1)
#print('DONE')
"""
ngf = n.h.List() #ngf = [0] * (2)
for num in range(2):
    ngf.append(BezaireNgf(num))   
print(type(ngf[0])) 
"""         
#ngf_1 = BezaireNgf(1)
#ngf_2 = BezaireNgf(2)

#h_interp('setpointer ngf.object(0).GapList.object(0).vgap, ngf.object(1).dend_2.v(0.5)')
#h_interp('setpointer ngf.object(1).GapList.object(0).vgap, ngf.object(0).dend_2.v(0.5)')
"""
#pc.source_var(_ref_v, source_global_index, sec=section)
pc.source_var(ngf.object(1).dend_2(0.5)._ref_v, 1, sec = ngf.object(1).dend_2)
pc.source_var(ngf.object(0).dend_2(0.5)._ref_v, 0, sec = ngf.object(0).dend_2)
gap_1 = n.h.gap3(ngf[0].dend_2(0.5))
gap_1.g = 0.17
gap_1.vgap = -70
gap_2 = n.h.gap3(ngf[1].dend_2(0.5))
gap_2.g = 0.17
gap_2.vgap = -70
#pc.target_var(targetPointProcess, _ref_target_variable, source_global_index)
pc.target_var(gap_1, ngf.object(0).dend_2(0.5)._ref_v, 0)
pc.target_var(gap_2, ngf.object(1).dend_2(0.5)._ref_v, 1)
pc.setup_transfer()
"""
#n.h.topology()
#n.h.Shape()
#ConList = n.h.List()
#a = n.h.Section(name='a')
#gap_0 = n.h.gapJ(0.5, sec=a)
"""
ConDict_Ngf2Ngf_syn = dict()

#con = n.h.NetCon(ngf_1.AllSections.object(0)(0.5)._ref_v, ngf_2.SynList.object(0))
con_1 = n.h.NetCon(ngf.object(0).soma(0.5)._ref_v, ngf.object(1).SynList_inh.object(1))
con_1.weight[0] = 0.5
con_1.delay = 10
con_1.threshold = -15
con = n.h.NetCon(ngf.object(0).soma(0.5)._ref_v, ngf.object(1).SynList_ex.object(0))
con.weight[0] = 0.5
con.delay = 100
con.threshold = -15
#print(ngf_1.SynList.count())
#print(ngf_1.soma.psection())
"""
"""
gap_0 = n.h.gapJ_m(ngf.object(0).dend_2(0.5))
gap_1 = n.h.gapJ_m(ngf.object(1).dend_2(0.5))
#n.h.setpointer(_ref_hocvar, 'POINTER_name', point_proces_object)
n.h.setpointer(ngf.object(1).dend_2(0.5)._ref_v, 'vgap', gap_0)
n.h.setpointer(ngf.object(0).dend_2(0.5)._ref_v, 'vgap', gap_1)
gap_0.g = 0.017 #0.0001
gap_1.g = 0.017 #0.0001
"""
"""
#nc = h.NetStim()
#ns = h.NetCon(nc, target...)
iclamp = n.h.IClamp(ngf[0].dend_2(0.5))
#iclamp = n.h.IClamp(ngf_1.soma(0.5))
iclamp.amp = 0.5 #3 # nA
iclamp.delay = 120 # ms
iclamp.dur = 5 # ms
"""
"""
iclamp_1 = n.h.IClamp(ngf[0].dend_2(0.5))
#iclamp = n.h.IClamp(ngf_1.soma(0.5))
iclamp_1.amp = 0.5 # nA
iclamp_1.delay = 220 # ms
iclamp_1.dur = 200 # ms


t = n.h.Vector()
t.record(n.h._ref_t)
ngf_2_dend_2 = n.h.Vector()
ngf_2_dend_2.record(ngf[1].dend_2(0.5)._ref_v)
ngf_1_dend_2 = n.h.Vector()
ngf_1_dend_2.record(ngf[0].dend_2(0.5)._ref_v)
ngf_1_soma = n.h.Vector()
ngf_1_soma.record(ngf[0].soma(0.5)._ref_v)
ngf_2_soma = n.h.Vector()
ngf_2_soma.record(ngf[1].soma(0.5)._ref_v)

n.h.tstop = 500
n.h.celsius = 37

n.h.run()

#pc.setup_transfer()

import matplotlib.pyplot as ppl
ppl.plot(t, ngf_2_dend_2, color="black")
ppl.plot(t, ngf_2_soma, color="blue")
ppl.show()
ppl.plot(t, ngf_1_dend_2, color="red")
ppl.plot(t, ngf_1_soma, color="green")
ppl.show()
"""