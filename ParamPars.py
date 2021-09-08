#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  3 14:46:56 2021

@author: user
"""
'''
import numpy as np
pre_gid = objc_p[0]["synapses_params"][0]["pre_gid"]
post_gid = objc_p[0]["synapses_params"][0]["post_gid"]
connected = (objc_p[0]["synapses_params"][0]["gmax"]) > 0
print("for pre_gid ", pre_gid, "and post_gid ", post_gid, "connection is ", connected)

import matplotlib.pyplot as plt
for thr in range(4):
    for pare in range(len(objc_p[thr]["synapses_params"])):
        pre_gid = objc_p[thr]["synapses_params"][pare]["pre_gid"]
        post_gid = objc_p[thr]["synapses_params"][pare]["post_gid"]
        connected = objc_p[thr]["synapses_params"][pare]["gmax"]   
        if connected > 0:
            plt.scatter(pre_gid, post_gid, color = 'black')
        else:
            plt.scatter(pre_gid, post_gid, color = 'white')
plt.show()
'''
'''
import numpy as np
ArtParamList = np.arange(5, 105, 5)
ArtParamList = np.array(ArtParamList)
for freq_param in ArtParamList:
    print(freq_param)
print(ArtParamList)
'''