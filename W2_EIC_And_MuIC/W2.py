import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib as mtl
data_muic = np.genfromtxt('MuIC_W2.csv', delimiter='\t')
data_eic = np.genfromtxt('EIC_W2.csv', delimiter='\t')

data_muic_log = np.log10(np.sqrt(data_muic))
data_eic_log = np.log10(np.sqrt(data_eic))

#data_muic_log = np.sqrt(data_muic)
#data_eic_log = np.sqrt(data_eic)

print("*****Load MuIC Successfully!*****")
print("size:", len(data_muic))
print("*****Load EIC Successfully!*****")
print("size:", len(data_eic))
plt.hist(data_muic_log, bins=100, edgecolor="b", histtype="step")
plt.hist(data_eic_log, bins=100, edgecolor="r", histtype="step")
plt.yscale('log')
plt.xlabel("$W [GeV]$")
plt.ylabel('counts')
xtick = np.log10([10,20,30,40,50,60,70,80,90,100,200,300,400,500,600,700,800,900,1000])
#xtick = [1,2,3,4,5,6]
xname = ['$10^1$','','','','','','','','','$10^2$','','','','','','','','','$10^3$']
#xname = ['$10^1$','$10^2$','$10^3$','$10^4$','$10^5$','$10^6$']
plt.xticks(xtick,xname)
plt.legend(['MuIC', 'EIC'])
plt.show()
