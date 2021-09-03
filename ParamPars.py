#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  3 14:46:56 2021

@author: user
"""
N_Ngf = 200
#Syn per Ngf = 10 realised as G_max *= 10
"Ngf_to_Ngf": {
     "gmax": 1.5 * 10, #nS 1 синапс идет за 10
     "gmax_std" : 0.7, # ?
     
     "Erev": -60, # in Bezaire -11
     "tau_rise": 4.8,
     "tau_decay": 32.03,

     "prob": 17.0/200.0, 
     
     "delay": 1.2,
     "delay_std" : 0.2,
     

     "sourse_compartment" : "soma",
     "target_compartment" : "dendrite_list",
 },

"gap_junctions_params" : {
    "ngf2ngf" : {
        "r" : 1e6,
        "r_std" : 10,
        "prob": 17/200,
        "compartment1" : "dendrite_list",
        "compartment2" : "dendrite_list",
    },
},
