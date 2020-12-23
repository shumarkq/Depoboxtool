#This is an example to investigate Vd difference as the function of (diamater and standard dievation matrix) for different deposition algorithm.
#Details can be found in Shu et al. 2021 (GMD).

from Models import *
from eval_luc import *
import pandas as pd
pd.set_option('display.max_rows', 1000)
import matplotlib.pylab as plt
from functions import bins
import seaborn as sns

dps = [0.01, 0.05, 0.1, 0.5, 1, 1.5, 2, 2.5, 3, 5, 10, 50]
sigmags = [1.01, 1.2, 1.4, 1.7, 2.0, 2.5]
num = 100

dataset = pd.read_csv('../Data/vong_grass.csv')
for data in [dataset.loc[11]]:
    dens = data.density
    temp = data.temp
    press = data.press
    RH = data.RH / 100.
    LAI = data.LAI
    Uh = data.Uh
    ustar = data.ustar
    wstar = data.wstarloc
    h = data.h
    d = data.d
    z0 = data.z0
    z = data.z
    Lo = data.Lo
    Vd_obs = data.Vd_cm
    index = data.index

    z01s = []
    qbases = []
    qimps = []
    qlais = []
    for dpg in dps:
        for sigmag in sigmags:

            #Z01 bw
            onebin = bins(num, sigmag, dpg)
            vds = []
            for mi in arange(0, len(onebin)):
                dpi, weighti = onebin[mi][0], onebin[mi][1]
                vdi = 100 * weighti * run(dens, dpi * 1e-6, temp, press, RH, ustar, Lo,
                    luc(wstar, ustar, Uh, z0, z, d, 'Z01').grass(), 'Z01').Z01()[0]
                vds.append(vdi)
            vdw = sum(vds)
            z01s.append([dpg, sigmag, vdw])

            #CMAQ modes
            qbase = 100 * mode(dens, dpg * 1e-6, temp, press, RH, ustar, Lo, luc(wstar, ustar, Uh, z0, z, d, 'CMAQmode').grass(), sigmag, 'CMAQmode', True, True, True).CMAQmode()[0]
            qbases.append(qbase)

            qimp = 100 * mode(dens, dpg * 1e-6, temp, press, RH, ustar, Lo, luc(wstar, ustar, Uh, z0, z, d, 'CMAQmode').grass(), sigmag, 'CMAQmode', True, True, False).CMAQmode()[0]
            qimps.append(qimp)

            qlai = 100 * laimode(dens, dpg * 1e-6, temp, press, RH, ustar, Lo, luc(wstar, ustar, Uh, z0, z, d, 'CMAQmodeLAI').grass(), sigmag, 'CMAQmodeLAI', True, True).CMAQmodeLAI()[0]
            qlais.append(qlai)



    df = pd.DataFrame({'dp':[i[0] for i in z01s],
                       'sigma':[i[1] for i in z01s],
                       'z01':[i[2] for i in z01s],
                       'base': qbases,
                       'imp': qimps,
                       'lai': qlais
                       }
                      )
    df.eval('z01_base = (z01-base)/base', inplace=True)
    df.eval('imp_base = (imp-base)/base', inplace=True)
    df.eval('lai_base = (lai-base)/base', inplace=True)

    def getplots(ax, data, label, title, xticklabels=True, yticklabels=False, cbar=False, axcb=None, axcblabel=None):
        mask = abs(data) > 0
        table = df[mask].pivot_table(values=label, index='sigma', columns='dp')
        heatmap = sns.heatmap(\
        table, vmin=-1, vmax=1, annot=True, fmt= '1.2f' ,linewidths=0.2, xticklabels=xticklabels, yticklabels=yticklabels,\
        cmap='RdBu_r', cbar=cbar, center=0, annot_kws={"size": 6}, ax=ax, cbar_ax=axcb,cbar_kws={'label':axcblabel})
        heatmap.tick_params(labelsize=6)
        cax = plt.gcf().axes[-1]
        cax.tick_params(labelsize=6)
        ax.tick_params(axis=u'both', which=u'both', length=0)
        ax.set_title(title, fontsize=8)
        ax.invert_yaxis()

    fig, axes = plt.subplots(nrows=1,ncols=4,figsize=(12,3.5),gridspec_kw={'wspace':0.05, 'hspace':0.5, 'width_ratios':[1,1,1,0.08]})
    p1=getplots(axes[0], df.z01_base, 'z01_base', '(Z01-PR11)/PR11', True, True, False)
    p2=getplots(axes[1], df.imp_base, 'imp_base', '(OFF-PR11)/PR11', True, False, False)
    p3=getplots(axes[2], df.lai_base, 'lai_base', '(VGLAI-PR11)/PR11', True, False, True, axes[3],'Grass')

    axes[0].set_ylabel(r'$\sigma_g$', fontsize=7)
    axes[1].set_ylabel('')
    axes[2].set_ylabel('')
    axes[1].set_xlabel('dp(%sm)'% r'$\mu$', fontsize=7)
    axes[0].set_xlabel('')
    axes[2].set_xlabel('')
    plt.savefig('../Fig/variance_grass.png',dpi=300)




