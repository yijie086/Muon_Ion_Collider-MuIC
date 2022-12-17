import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib as mtl
data_onejet = np.genfromtxt('onejetKppp.csv', delimiter='\t')
length_onejet = int(len(data_onejet)/4.0)
print("*****Load onejet Successfully!*****")
print("size:", length_onejet, "* 4 (eta, pt, Q^2, x)")

onejet_theta = data_onejet[0: length_onejet-1]
onejet_pt = data_onejet[length_onejet: 2*length_onejet-1]
onejet_Q2 = data_onejet[2*length_onejet: 3*length_onejet-1]
onejet_x = data_onejet[3*length_onejet: 4*length_onejet-1]
print("szie of theta:", len(onejet_theta), data_onejet[2])
print("szie of pt:", len(onejet_theta), onejet_pt.max())
print("szie of Q^2:", len(onejet_theta))
print("szie of x:", len(onejet_theta))
df = pd.DataFrame()

onejet_theta_rescale = np.pi/2.0 + (np.pi/8.0)*(np.log(np.tan(onejet_theta/2.0)))

onejet_eta_015w = []
onejet_eta_0154 = []
onejet_eta_0143 = []
onejet_eta_01w3 = []
onejet_eta_124w = []
onejet_eta_1243 = []
onejet_eta_1232 = []
onejet_eta_12w2 = []
onejet_eta_233w = []
onejet_eta_2332 = []
onejet_eta_2321 = []
onejet_eta_23w1 = []
onejet_eta_342w = []
onejet_eta_3421 = []
onejet_eta_3410 = []
onejet_eta_451w = []
onejet_eta_4510 = []
onejet_eta_w5w1 = []
onejet_pt_015w = []
onejet_pt_0154 = []
onejet_pt_0143 = []
onejet_pt_01w3 = []
onejet_pt_124w = []
onejet_pt_1243 = []
onejet_pt_1232 = []
onejet_pt_12w2 = []
onejet_pt_233w = []
onejet_pt_2332 = []
onejet_pt_2321 = []
onejet_pt_23w1 = []
onejet_pt_342w = []
onejet_pt_3421 = []
onejet_pt_3410 = []
onejet_pt_451w = []
onejet_pt_4510 = []
onejet_pt_w5w1 = []

pt_max = 400

for j in range(length_onejet-1):
    if onejet_Q2[j]>1 and onejet_Q2[j]<10:
        if onejet_x[j]<0.00001:
            onejet_eta_015w.append(onejet_theta_rescale[j])
            onejet_pt_015w.append(onejet_pt[j])
        elif onejet_x[j]<0.0001 and onejet_x[j]>0.00001:
            onejet_eta_0154.append(onejet_theta_rescale[j])
            onejet_pt_0154.append(onejet_pt[j])
        elif onejet_x[j]<0.001 and onejet_x[j]>0.0001:
            onejet_eta_0143.append(onejet_theta_rescale[j])
            onejet_pt_0143.append(onejet_pt[j])
        else:
            onejet_eta_01w3.append(onejet_theta_rescale[j])
            onejet_pt_01w3.append(onejet_pt[j])

    elif onejet_Q2[j]>10 and onejet_Q2[j]<100:
        if onejet_x[j] < 0.0001:
            onejet_eta_124w.append(onejet_theta_rescale[j])
            onejet_pt_124w.append(onejet_pt[j])
        elif onejet_x[j] < 0.001 and onejet_x[j] > 0.0001:
            onejet_eta_1243.append(onejet_theta_rescale[j])
            onejet_pt_1243.append(onejet_pt[j])
        elif onejet_x[j] < 0.01 and onejet_x[j] > 0.001:
            onejet_eta_1232.append(onejet_theta_rescale[j])
            onejet_pt_1232.append(onejet_pt[j])
        else:
            onejet_eta_12w2.append(onejet_theta_rescale[j])
            onejet_pt_12w2.append(onejet_pt[j])

    elif onejet_Q2[j]>100 and onejet_Q2[j]<1000:
        if onejet_x[j] < 0.001:
            onejet_eta_233w.append(onejet_theta_rescale[j])
            onejet_pt_233w.append(onejet_pt[j])
        elif onejet_x[j] < 0.01 and onejet_x[j] > 0.001:
            onejet_eta_2332.append(onejet_theta_rescale[j])
            onejet_pt_2332.append(onejet_pt[j])
        elif onejet_x[j] < 0.1 and onejet_x[j] > 0.01:
            onejet_eta_2321.append(onejet_theta_rescale[j])
            onejet_pt_2321.append(onejet_pt[j])
        else:
            onejet_eta_23w1.append(onejet_theta_rescale[j])
            onejet_pt_23w1.append(onejet_pt[j])

    elif onejet_Q2[j]>1000 and onejet_Q2[j]<10000:
        if onejet_x[j] < 0.01:
            onejet_eta_342w.append(onejet_theta_rescale[j])
            onejet_pt_342w.append(onejet_pt[j])
        elif onejet_x[j] < 0.1 and onejet_x[j] > 0.01:
            onejet_eta_3421.append(onejet_theta_rescale[j])
            onejet_pt_3421.append(onejet_pt[j])
        else:
            onejet_eta_3410.append(onejet_theta_rescale[j])
            onejet_pt_3410.append(onejet_pt[j])

    elif onejet_Q2[j]>10000 and onejet_Q2[j]<100000:
        if onejet_x[j] < 0.1:
            onejet_eta_451w.append(onejet_theta_rescale[j])
            onejet_pt_451w.append(onejet_pt[j])
        else:
            onejet_eta_4510.append(onejet_theta_rescale[j])
            onejet_pt_4510.append(onejet_pt[j])
    elif onejet_Q2[j]>100000:
        onejet_eta_w5w1.append(onejet_theta_rescale[j])
        onejet_pt_w5w1.append(onejet_pt[j])

onejet_theta_select = []
onejet_pt_select = []

rbins = np.linspace(0, np.log10(pt_max), 120)
abins = np.linspace(0, 2*np.pi, 120)

h2d015w, _, _ = np.histogram2d(onejet_eta_015w, onejet_pt_015w, bins=(abins, np.power(10, rbins)))
h2d0154, _, _ = np.histogram2d(onejet_eta_0154, onejet_pt_0154, bins=(abins, np.power(10, rbins)))
h2d0143, _, _ = np.histogram2d(onejet_eta_0143, onejet_pt_0143, bins=(abins, np.power(10, rbins)))
h2d01w3, _, _ = np.histogram2d(onejet_eta_01w3, onejet_pt_01w3, bins=(abins, np.power(10, rbins)))
h2d124w, _, _ = np.histogram2d(onejet_eta_124w, onejet_pt_124w, bins=(abins, np.power(10, rbins)))
h2d1243, _, _ = np.histogram2d(onejet_eta_1243, onejet_pt_1243, bins=(abins, np.power(10, rbins)))
h2d1232, _, _ = np.histogram2d(onejet_eta_1232, onejet_pt_1232, bins=(abins, np.power(10, rbins)))
h2d12w2, _, _ = np.histogram2d(onejet_eta_12w2, onejet_pt_12w2, bins=(abins, np.power(10, rbins)))
h2d233w, _, _ = np.histogram2d(onejet_eta_233w, onejet_pt_233w, bins=(abins, np.power(10, rbins)))
h2d2332, _, _ = np.histogram2d(onejet_eta_2332, onejet_pt_2332, bins=(abins, np.power(10, rbins)))
h2d2321, _, _ = np.histogram2d(onejet_eta_2321, onejet_pt_2321, bins=(abins, np.power(10, rbins)))
h2d23w1, _, _ = np.histogram2d(onejet_eta_23w1, onejet_pt_23w1, bins=(abins, np.power(10, rbins)))
h2d342w, _, _ = np.histogram2d(onejet_eta_342w, onejet_pt_342w, bins=(abins, np.power(10, rbins)))
h2d3421, _, _ = np.histogram2d(onejet_eta_3421, onejet_pt_3421, bins=(abins, np.power(10, rbins)))
h2d3410, _, _ = np.histogram2d(onejet_eta_3410, onejet_pt_3410, bins=(abins, np.power(10, rbins)))
h2d451w, _, _ = np.histogram2d(onejet_eta_451w, onejet_pt_451w, bins=(abins, np.power(10, rbins)))
h2d4510, _, _ = np.histogram2d(onejet_eta_4510, onejet_pt_4510, bins=(abins, np.power(10, rbins)))
h2dw5w1, _, _ = np.histogram2d(onejet_eta_w5w1, onejet_pt_w5w1, bins=(abins, np.power(10, rbins)))

rtick = [0, 1, 2]
rname = [1, 10, 100]
etaname = [4, 3, 2, 1, 0, -1, -2, -3, -4, -5, -6, -7, -8, 7, 6, 5]
etatick = [0, np.pi/8, 2*np.pi/8, 3*np.pi/8, 4*np.pi/8, 5*np.pi/8, 6*np.pi/8, 7*np.pi/8, 8*np.pi/8, 9*np.pi/8, 10*np.pi/8, 11*np.pi/8, 12*np.pi/8, 13*np.pi/8, 14*np.pi/8, 15*np.pi/8]
nonename = []
A, R = np.meshgrid(abins, rbins)
fig, ax = plt.subplots(6, 6, subplot_kw=dict(polar=True))
plt.subplots_adjust(wspace=0, hspace=0.1)
maxcolor = 100
alphapara = 0.2

ax[5, 0].pcolor(A, R, h2d015w.T, norm=mtl.colors.LogNorm(vmin=1, vmax=maxcolor))
plt.sca(ax[5, 0])
ax[5, 0].grid(True, alpha=alphapara)
plt.yticks(rtick, nonename)
plt.xticks(etatick, nonename)
ax[5, 1].pcolor(A, R, h2d0154.T, norm=mtl.colors.LogNorm(vmin=1, vmax=maxcolor))
plt.sca(ax[5, 1])
ax[5, 1].grid(True, alpha=alphapara)
plt.yticks(rtick, nonename)
plt.xticks(etatick, nonename)
ax[5, 2].pcolor(A, R, h2d0143.T, norm=mtl.colors.LogNorm(vmin=1, vmax=maxcolor))
plt.sca(ax[5, 2])
ax[5, 2].grid(True, alpha=alphapara)
plt.yticks(rtick, nonename)
plt.xticks(etatick, nonename)
ax[5, 3].pcolor(A, R, h2d01w3.T, norm=mtl.colors.LogNorm(vmin=1, vmax=maxcolor))
plt.sca(ax[5, 3])
ax[5, 3].grid(True, alpha=alphapara)
plt.yticks(rtick, nonename)
plt.xticks(etatick, nonename)
ax[4, 1].pcolor(A, R, h2d124w.T, norm=mtl.colors.LogNorm(vmin=1, vmax=maxcolor))
plt.sca(ax[4, 1])
ax[4, 1].grid(True, alpha=alphapara)
plt.yticks(rtick, nonename)
plt.xticks(etatick, nonename)
ax[4, 2].pcolor(A, R, h2d1243.T, norm=mtl.colors.LogNorm(vmin=1, vmax=maxcolor))
plt.sca(ax[4, 2])
ax[4, 2].grid(True, alpha=alphapara)
plt.yticks(rtick, nonename)
plt.xticks(etatick, nonename)
ax[4, 3].pcolor(A, R, h2d1232.T, norm=mtl.colors.LogNorm(vmin=1, vmax=maxcolor))
plt.sca(ax[4, 3])
ax[4, 3].grid(True, alpha=alphapara)
plt.yticks(rtick, nonename)
plt.xticks(etatick, nonename)
ax[4, 4].pcolor(A, R, h2d12w2.T, norm=mtl.colors.LogNorm(vmin=1, vmax=maxcolor))
plt.sca(ax[4, 4])
ax[4, 4].grid(True, alpha=alphapara)
plt.yticks(rtick, nonename)
plt.xticks(etatick, nonename)

ax[3, 2].pcolor(A, R, h2d233w.T, norm=mtl.colors.LogNorm(vmin=1, vmax=maxcolor))
plt.sca(ax[3, 2])
ax[3, 2].grid(True, alpha=alphapara)
plt.yticks(rtick, nonename)
plt.xticks(etatick, nonename)
ax[3, 3].pcolor(A, R, h2d2332.T, norm=mtl.colors.LogNorm(vmin=1, vmax=maxcolor))
plt.sca(ax[3, 3])
ax[3, 3].grid(True, alpha=alphapara)
plt.yticks(rtick, nonename)
plt.xticks(etatick, nonename)
ax[3, 4].pcolor(A, R, h2d2321.T, norm=mtl.colors.LogNorm(vmin=1, vmax=maxcolor))
plt.sca(ax[3, 4])
ax[3, 4].grid(True, alpha=alphapara)
plt.yticks(rtick, nonename)
plt.xticks(etatick, nonename)
ax[3, 5].pcolor(A, R, h2d23w1.T, norm=mtl.colors.LogNorm(vmin=1, vmax=maxcolor))
plt.sca(ax[3, 5])
ax[3, 5].grid(True, alpha=alphapara)
plt.yticks(rtick, nonename)
plt.xticks(etatick, nonename)

ax[2, 3].pcolor(A, R, h2d342w.T, norm=mtl.colors.LogNorm(vmin=1, vmax=maxcolor))
plt.sca(ax[2, 3])
ax[2, 3].grid(True, alpha=alphapara)
plt.yticks(rtick, nonename)
plt.xticks(etatick, nonename)
ax[2, 4].pcolor(A, R, h2d3421.T, norm=mtl.colors.LogNorm(vmin=1, vmax=maxcolor))
plt.sca(ax[2, 4])
ax[2, 4].grid(True, alpha=alphapara)
plt.yticks(rtick, nonename)
plt.xticks(etatick, nonename)
ax[2, 5].pcolor(A, R, h2d3410.T, norm=mtl.colors.LogNorm(vmin=1, vmax=maxcolor))
plt.sca(ax[2, 5])
ax[2, 5].grid(True, alpha=alphapara)
plt.yticks(rtick, nonename)
plt.xticks(etatick, nonename)

ax[1, 4].pcolor(A, R, h2d451w.T, norm=mtl.colors.LogNorm(vmin=1, vmax=maxcolor))
plt.sca(ax[1, 4])
ax[1, 4].grid(True, alpha=alphapara)
plt.yticks(rtick, nonename)
plt.xticks(etatick, nonename)
ax[1, 5].pcolor(A, R, h2d4510.T, norm=mtl.colors.LogNorm(vmin=1, vmax=maxcolor))
plt.sca(ax[1, 5])
ax[1, 5].grid(True, alpha=alphapara)
plt.yticks(rtick, nonename)
plt.xticks(etatick, nonename)

ax[0, 5].pcolor(A, R, h2dw5w1.T, norm=mtl.colors.LogNorm(vmin=1, vmax=maxcolor))
plt.sca(ax[0, 5])
ax[0, 5].grid(True, alpha=alphapara)
plt.yticks(rtick, nonename)
plt.xticks(etatick, nonename)

plt.sca(ax[0, 0])
plt.axis('off')
plt.sca(ax[0, 1])
plt.axis('off')
plt.sca(ax[0, 2])
plt.axis('off')
plt.sca(ax[0, 3])
plt.axis('off')
plt.sca(ax[0, 4])
plt.axis('off')
plt.sca(ax[1, 0])
plt.axis('off')
plt.sca(ax[1, 1])
plt.axis('off')
plt.sca(ax[1, 2])
plt.axis('off')
plt.sca(ax[1, 3])
plt.axis('off')
plt.sca(ax[2, 0])
plt.axis('off')
plt.sca(ax[2, 1])
plt.axis('off')
plt.sca(ax[2, 2])
plt.axis('off')
plt.sca(ax[3, 0])
plt.axis('off')
plt.sca(ax[3, 1])
plt.axis('off')
plt.sca(ax[4, 0])
plt.axis('off')
plt.sca(ax[5, 4])
plt.axis('off')
plt.sca(ax[5, 5])
plt.axis('off')
plt.sca(ax[4, 5])
plt.axis('off')

plt.show()