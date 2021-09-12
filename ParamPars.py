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
'''
import math
import numpy as np

a = [1.1, 2.3, 3.4, 4.6, 5.8]
a = np.array(a)
b = np.zeros(30)
print(a[-1])
print(np.diff(a))

b[np.floor(a * 2).astype(np.int)] +=1

print("a*2", np.array(a)*2)
print("floor(a*2)", np.floor(a[:]*2))
print("b", b)

print(a[a<4])

import os
path = "./Results/"
for file in os.listdir(path):
    if file.split(".")[-1] != "hdf5":
        continue

    freq = file.split("state")[1]
    print("freq 1", freq)
    freq = float(freq.split(".")[0])
    print("freq 2", freq)
'''
#Образец: создание полярных диаграмм для циркулярной статистики.
'''
import numpy as np
import matplotlib.pyplot as plt
rads = np.arange(0, (2*np.pi), 0.01)
radian = rads * np.random.rand(rads.size)
r = (np.zeros_like(radian)+0.5) * np.random.rand(rads.size)
fig = plt.figure()
axes = fig.add_subplot(1, 1, 1, polar=True)
axes.set_title('Circle in polar format:r=R')
for i in range(len(r)):
    axes.plot([0, radian[i]], [0, r[i]], "-") 
    axes.plot(radian[i], r[i], ".")
plt.show()
'''
#Работа с директориями:
import os
import shutil
path = "./Results/" #os.getcwd()
print(path)
for file in os.listdir(path):
    file_name = file.split(".")
    print("file_name: ", file_name)
    if file_name[-1] != "hdf5":
        continue
    file_name = file_name[0].split("_") 
    freq_of_file = file_name[0]
    print("freq_of_file: ", freq_of_file)
    destination = path + freq_of_file
    try:
        os.mkdir(destination)
    except FileExistsError:
        pass
    shutil.move(path+file, destination)
    


