import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib as mtl
data_muic = np.genfromtxt('MuIC_W2.csv', delimiter='\t')
data_eic = np.genfromtxt('EIC_W2.csv', delimiter='\t')

data_muic_log = np.log10(data_muic)
data_eic_log = np.log10(data_eic)

print("*****Load MuIC Successfully!*****")
print("size:", len(data_muic))
print("*****Load EIC Successfully!*****")
print("size:", len(data_eic))
plt.hist(data_muic_log, bins=100, edgecolor="b", histtype="step")
plt.hist(data_eic_log, bins=100, edgecolor="r", histtype="step")
plt.yscale('log')
plt.xlabel("$W^2 [GeV^2]$")
plt.ylabel('counts')
xtick = [1,2,3,4,5,6]
xname = ['$10^1$','$10^2$','$10^3$','$10^4$','$10^5$','$10^6$']
plt.xticks(xtick, xname)
plt.legend(['MuIC', 'EIC'])
plt.show()
