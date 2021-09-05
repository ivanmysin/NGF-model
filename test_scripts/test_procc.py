import numpy as np


binned_spiketrain = []
for _ in range(10):
    sp = np.random.rand(100)
    binned_spiketrain.append(sp)

binned_spiketrain = np.vstack(binned_spiketrain)

print(binned_spiketrain.shape)

corrM = np.corrcoef(binned_spiketrain)
print(corrM)