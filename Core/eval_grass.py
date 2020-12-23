#This is an example to evaluate deposition schemes using site measurements.
#Details can be found in Shu et al. 2021 (GMD).

from Models import *
from eval_luc import *
import pandas as pd
pd.set_option('display.max_rows', 1000)
from datetime import datetime
import matplotlib.pylab as plt
from functions import bins


dataset = pd.read_csv('../Data/vong_grass.csv')
time = dataset.loctime
for dimi in [0.54]:

    bins1 = bins(100, 1.01, dimi, plot=False)
    bins2 = bins(100, 1.7, dimi, plot = False)
    bins3 = bins(100, 2.5, dimi, plot=False)

    for data in [dataset]:
        dens = data.density
        dim = dimi * 1e-6
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

        #one dp
        dpone = {}
        dpone['Z01'] = [100 * run(dens[i], dim, temp[i], press[i], RH[i], ustar[i], Lo[i],
               luc(wstar[i], ustar[i], Uh[i], z0[i], z[i], d[i], 'Z01').grass(), 'Z01').Z01()[0] for i in index]
        dpone['CMAQ'] = [100 * run(dens[i], dim, temp[i], press[i], RH[i], ustar[i], Lo[i],
                luc(wstar[i], ustar[i], Uh[i], z0[i], z[i], d[i], 'CMAQ').grass(), 'CMAQ').CMAQ()[0] for i in index]
        dpone['CMAQLAI'] = [100 * run(dens[i], dim, temp[i], press[i], RH[i], ustar[i], Lo[i],
               luc(wstar[i], ustar[i], Uh[i], z0[i], z[i], d[i], 'CMAQLAI').grass(), 'CMAQLAI').CMAQLAI()[0] for i in index]



        #dpbins
        sigmagbins = {'1.7':bins2, '2.5':bins3}
        dpbins = {}

        for depmodel in ['Z01','CMAQ','CMAQLAI']:
            dpbins[depmodel] = {}
            for modesigma, modeinterval in sigmagbins.items():
                weights = []
                for di in data.index:
                    vds = []
                    for mi in arange(0, len(modeinterval)):
                        dpi, weighti = modeinterval[mi][0], modeinterval[mi][1]

                        if depmodel in ['Z01']:
                            vdi = 100 * weighti * run(dens[di], dpi * 1e-6, temp[di], press[di], RH[di], ustar[di], Lo[di],
                                  luc(wstar[di], ustar[di], Uh[di], z0[di], z[di], d[di], depmodel).grass(), depmodel).Z01()[0]

                        elif depmodel in ['CMAQ']:
                            vdi = 100 * weighti * run(dens[di], dpi * 1e-6, temp[di], press[di], RH[di], ustar[di], Lo[di],
                                  luc(wstar[di], ustar[di], Uh[di], z0[di], z[di], d[di], depmodel).grass(), depmodel).CMAQ()[0]

                        elif depmodel in ['CMAQLAI']:
                            vdi = 100 * weighti * run(dens[di], dpi * 1e-6, temp[di], press[di], RH[di], ustar[di], Lo[di],
                                  luc(wstar[di], ustar[di], Uh[di], z0[di], z[di], d[di], depmodel).grass(),depmodel).CMAQLAI()[0]

                        else:
                            print('no model was selected')

                        vds.append(vdi)
                    weights.append(sum(vds))
                dpbins[depmodel][modesigma] = weights

        #dpmode
        dpmode = {}
        dpmode['CMAQ_base'] = {}
        for sigmag in [1.7, 2.5]:
            dpmode['CMAQ_base'][sigmag] = [100 * mode(dens[i], dim, temp[i], press[i], RH[i], ustar[i], Lo[i],
                luc(wstar[i], ustar[i], Uh[i], z0[i], z[i], d[i], 'CMAQmode').grass(), sigmag, 'CMAQmode',True, True, True).CMAQmode()[0] for i in index]

        dpmode['CMAQ_imp'] = {}
        for sigmag in [1.7, 2.5]:
            dpmode['CMAQ_imp'][sigmag] = [100 * mode(dens[i], dim, temp[i], press[i], RH[i], ustar[i], Lo[i],
                luc(wstar[i], ustar[i], Uh[i], z0[i], z[i], d[i], 'CMAQmode').grass(), sigmag, 'CMAQmode',True, True, False).CMAQmode()[0] for i in index]

        dpmode['CMAQ_lai'] = {}
        for sigmag in [1.7, 2.5]:
            dpmode['CMAQ_lai'][sigmag] = [100 * laimode(dens[i], dim, temp[i], press[i], RH[i], ustar[i], Lo[i],
                luc(wstar[i], ustar[i], Uh[i], z0[i], z[i], d[i], 'CMAQmodeLAI').grass(), sigmag, 'CMAQmodeLAI',True, True).CMAQmodeLAI()[0] for i in index]


        plt.close()
        fig = plt.figure(figsize=(8, 4))
        ax = fig.add_axes([0.15, 0.12, 0.5, 0.8])
        markersize = 2
        ax.plot(time, Vd_obs, color = 'grey', lw =1 , ls = '-', marker='o', markersize = markersize, mfc='white', label='Obs (Vong et al 2004)')


        ax.plot(time, dpone['Z01'], color = 'red', lw = 1, ls = '-', marker='o', markersize = markersize, mfc='white',  label = 'Z01 single diameter' )
        ax.plot(time, dpbins['Z01']['1.7'], color = 'red', lw = 1, ls = '-', marker = '^', markersize = markersize, mfc='white',  label = 'Z01 sectional, %s=1.7' % r'$\sigma$' )
        ax.plot(time, dpbins['Z01']['2.5'], color = 'red', lw = 1, ls = '-', marker = 'v', markersize = markersize, mfc='white',  label = 'Z01 sectional, %s=2.5' % r'$\sigma$' )

        ax.plot(time, dpone['CMAQ'], color = 'blue', lw = 1, ls = '-', marker='o', markersize = markersize, mfc='white',  label = 'PX(base&IMP) single diameter' )
        ax.plot(time, dpbins['CMAQ']['1.7'], color = 'blue', lw = 1, ls = '-', marker = '^', markersize = markersize, mfc='white',  label = 'PX(base&IMP) sectional, %s=1.7' % r'$\sigma$')
        ax.plot(time, dpbins['CMAQ']['2.5'], color = 'blue', lw = 1, ls = '-', marker = 'v', markersize = markersize, mfc='white',  label = 'PX(base&IMP) sectional, %s=2.5' % r'$\sigma$')
        ax.plot(time, [i / 1. for i in dpmode['CMAQ_base'][1.7]], color='blue', lw=1, ls='--', marker='^', markersize=markersize, mfc='white', label='PX(base) modal, %s=1.7' % r'$\sigma$')
        ax.plot(time, [i / 1. for i in dpmode['CMAQ_base'][2.5]], color='blue', lw=1, ls='--', marker='v', markersize=markersize, mfc='white', label='PX(base) modal, %s=2.5' % r'$\sigma$')
        ax.plot(time, [i / 1. for i in dpmode['CMAQ_imp'][1.7]], color='blue', lw=1, ls='--', marker='<', markersize=markersize, mfc='white', label='PX(IMP) modal, %s=1.7' % r'$\sigma$')
        ax.plot(time, [i / 1. for i in dpmode['CMAQ_imp'][2.5]], color='blue', lw=1, ls='--', marker='>', markersize=markersize, mfc='white', label='PX(IMP) modal, %s=2.5' % r'$\sigma$')

        ax.plot(time, dpone['CMAQLAI'], color = 'green', lw = 1, ls = '-', marker='o', markersize = markersize, mfc='white',  label = 'PX(LAI) single diameter' )
        ax.plot(time, dpbins['CMAQLAI']['1.7'], color = 'green', lw = 1, ls = '-', marker = '^', markersize = markersize, mfc='white',  label = 'PX(LAI) sectional, %s=1.7' % r'$\sigma$')
        ax.plot(time, dpbins['CMAQLAI']['2.5'], color = 'green', lw = 1, ls = '-', marker = 'v', markersize = markersize, mfc='white',  label = 'PX(LAI) sectional, %s=2.5' % r'$\sigma$')
        ax.plot(time, [i / 1. for i in dpmode['CMAQ_lai'][1.7]], color='green', lw=1, ls='--', marker='^', markersize=markersize, mfc='white', label='PX(LAI) modal, %s=1.7'% r'$\sigma$')
        ax.plot(time, [i / 1. for i in dpmode['CMAQ_lai'][2.5]], color='green', lw=1, ls='--', marker='v', markersize=markersize, mfc='white', label='PX(LAI) modal, %s=2.5'% r'$\sigma$')

        locs, labels = plt.xticks()
        ax.set_ylabel('Vd (cm/s)')
        ax.set_xlabel('Hour of day (local time at Oregon)')
        ax.set_title('           Grass (Density = 1500 kg/m3, Diameter = 0.52 %sm)' % r'$\mu$', fontsize=9)
        ax.legend(ncol=1, bbox_to_anchor=(1.55,.87), fontsize=8)
        fig.savefig('../Fig/Eval_grass.png', dpi=300)








