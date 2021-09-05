import os

dir = "./Results/"
#files = [ "params.pickle", ]

for file in os.listdir(dir):
    if file.split(".")[-1] != "pickle":
         continue
    filepath = dir + file
    os.system('mpiexec -n 4 python3 base_model.py {}'.format(filepath) )