import os

path = "./Results/"
#files = [ "params.pickle", ]

for file in os.listdir(path):
    if file.split(".")[-1] != "pickle":
         continue
    filepath = path + file
    os.system('mpiexec -n 4 python3 base_model.py {}'.format(filepath) )
    