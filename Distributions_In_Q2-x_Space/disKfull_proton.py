import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib as mtl
data_proton = np.genfromtxt('protonKfull.csv', delimiter='\t')
length_proton = int(len(data_proton)/4.0)
print("*****Load proton Successfully!*****")
print("size:", length_proton, "* 4 (eta, pt, Q^2, x)")

proton_theta = data_proton[0: length_proton-1]
proton_pt = data_proton[length_proton: 2*length_proton-1]
proton_Q2 = data_proton[2*length_proton: 3*length_proton-1]
proton_x = data_proton[3*length_proton: 4*length_proton-1]
print("szie of theta:", len(proton_theta), data_proton[2])
print("szie of pt:", len(proton_theta), proton_pt.max())
print("szie of Q^2:", len(proton_theta))
print("szie of x:", len(proton_theta))
df = pd.DataFrame()

proton_theta_rescale = np.pi/2.0 + (np.pi/8.0)*(np.log(np.tan(proton_theta/2.0)))

pt_max = 400

rbins = np.linspace(0, np.log10(pt_max), 120)
abins = np.linspace(0, 2*np.pi, 120)
hist, _, _ = np.histogram2d(proton_theta_rescale, proton_pt, bins=(abins, np.power(10, rbins)))
A, R = np.meshgrid(abins, rbins)
fig, ax = plt.subplots(subplot_kw=dict(polar=True))

pc = ax.pcolor(A, R, hist.T, norm=mtl.colors.LogNorm(vmin=1, vmax=20))
plt.grid(True, ls=':')

rtick = [0, 1, 2]
rname = ['1GeV/c', '10GeV/c', '$p_T$=100GeV/c']
etaname = [4, 3, 2, 1, 0, -1, -2, -3, -4, -5, -6, -7, -8, 7, 6, 5]
etatick = [0, np.pi/8, 2*np.pi/8, 3*np.pi/8, 4*np.pi/8, 5*np.pi/8, 6*np.pi/8, 7*np.pi/8, 8*np.pi/8, 9*np.pi/8, 10*np.pi/8, 11*np.pi/8, 12*np.pi/8, 13*np.pi/8, 14*np.pi/8, 15*np.pi/8]

plt.yticks(rtick, rname)
plt.xlabel("$\eta$")
plt.xticks(etatick, etaname)
#ax[0, 0].set_rgrids([1, 10, 100, 400])
fig.colorbar(pc, ticks=[1, 10, 100, 1000, 10000, 100000, 1000000])
plt.show()


