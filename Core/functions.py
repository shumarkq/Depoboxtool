#include functions that can be used for box deposition analysis.
#More functions will be added for different purpose.

from numpy import *
import matplotlib.pyplot as plt
from math import erf

#create number concentration lognorm size distribution
def bins(num, sigmag, dpg, plot=False):
    num = num
    sigmag=sigmag
    dpg = dpg
    plot = plot
    #dps is very important to decide bins lower and upper bound
    #here we choose -4 to 3 to represent 10**-4 and 10**3 mu dp range
    dps = logspace(-4,3,num)
    lower = dps[:-1]
    upper = dps[1:]
    interval = [i for i in zip(lower,upper)]
    Ndp = [num/2. * (erf(log(dp[1]/dpg)/(sqrt(2)*log(sigmag))) - erf(log(dp[0]/dpg)/(sqrt(2)*log(sigmag)))) for dp in interval]
    wt = [i/num for i in Ndp]
    dpmedian = [median(interval[i]) for i in arange(0,len(interval))]

    if plot:
        plt.close()
        fig = plt.figure(figsize=(8,6))
        ax = fig.add_axes([.15,.15,.75,.75])
        ax.plot(dpmedian,wt,marker='o')
        ax.set_xscale('log')
        plt.show()

    return [(dpmedian[i], wt[i]) for i in arange(0,len(interval))]

#create mass concentration lognorm size distribution
def mbins(num, sigmag, dpg, dens, plot=False):
    num = num
    sigmag=sigmag
    dpg = dpg
    dens = dens
    plot = plot
    dps = logspace(-4,3,num)
    lower = dps[:-1]
    upper = dps[1:]
    interval = [i for i in zip(lower,upper)]
    Ndp = [num/2. * (erf(log(dp[1]/dpg)/(sqrt(2)*log(sigmag))) - erf(log(dp[0]/dpg)/(sqrt(2)*log(sigmag)))) for dp in interval]
    dpmedian = [median(interval[i]) for i in arange(0,len(interval))]
    massmedian = [dpmedian[i] ** 3 * pi / 6 * dens * Ndp[i] for i in arange(0, len(Ndp))]
    totalmass = sum(massmedian)
    wt = [i/totalmass for i in massmedian]

    if plot:
        plt.close()
        fig = plt.figure(figsize=(8,6))
        ax = fig.add_axes([.15,.15,.75,.75])
        ax.plot(dpmedian,wt, marker='o')
        ax.set_xscale('log')
        plt.show()

    return [(dpmedian[i], wt[i]) for i in arange(0,len(interval))]





#get medians for df data
def getmedian(data, name):
    newdata = data[data[name] > 0]
    ti3 = (newdata.loctime == '3:00:00') | (newdata.loctime == '4:00:00') | (newdata.loctime == '5:00:00') | (newdata.loctime == '6:00:00')
    ti7 = (newdata.loctime == '7:00:00') | (newdata.loctime == '8:00:00') | (newdata.loctime == '9:00:00') | (newdata.loctime == '10:00:00')
    ti11 = (newdata.loctime == '11:00:00') | (newdata.loctime == '12:00:00') | (newdata.loctime == '13:00:00') | (newdata.loctime == '14:00:00')
    ti15 = (newdata.loctime == '15:00:00') | (newdata.loctime == '16:00:00') | (newdata.loctime == '17:00:00') | (newdata.loctime == '18:00:00')
    ti19 = (newdata.loctime == '19:00:00') | (newdata.loctime == '20:00:00') | (newdata.loctime == '21:00:00') | (newdata.loctime == '22:00:00')
    ti23 = (newdata.loctime == '23:00:00') | (newdata.loctime == '0:00:00') | (newdata.loctime == '1:00:00') | (newdata.loctime == '2:00:00')
    box3 = newdata[ti3]
    box7 = newdata[ti7]
    box11 = newdata[ti11]
    box15 = newdata[ti15]
    box19 = newdata[ti19]
    box23 = newdata[ti23]
    out = [median(i) for i in [box3[name],box7[name],box11[name],box15[name],box19[name],box23[name]]]
    return out

#"""
#get NMBS
forNMBF = {'Z01': Z01, 'ZS14': ZS14, 'CMAQ': CMAQ, 'CMAQstd': CMAQstd}
NMBFs = {}
for k, v in forNMBF.items():
    if mean(v) >= mean(Vd_obs):
        NMBFs[k] = mean(v) / mean(Vd_obs) - 1
    else:
        NMBFs[k] = 1 - mean(Vd_obs) / mean(v)

for k, v in NMBFs.items():
    print(k, v)
#"""
