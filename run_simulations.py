import os

filepath = "./Results/params.pickle"
os.system('mpiexec -n 4 python3 base_model.py {}'.format(filepath) )