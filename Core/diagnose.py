#This is an example to diagnose the influence of different impaction term in different papers.
#Details can be found in Shu et al. 2021 (GMD).
from Models import *
from eval_luc import *
import pandas as pd
pd.set_option('display.max_rows', 1000)
from datetime import datetime
import matplotlib.pylab as plt
from functions import bins, getmedian


dataset = pd.read_csv('../Data/matsuda_df.csv')
date = dataset.Date
time = dataset.loctime
thetime = []
for i in dataset.index:
    newdate = datetime.strptime(date[i], '%m/%d')
    newtime = datetime.strptime(time[i], '%X')
    newdatetime = datetime.combine(newdate.date(),newtime.time())
    newdatetime = newdatetime.replace(year=2010)
    thetime.append(newdatetime)
dataset['thetime'] = thetime


for dimi in [0.48]:

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
        index = data.index

        #dpbins
        sigmagbins = {'1.01':bins1, '1.7':bins2, '2.5':bins3}
        for soustokes in ['Slinn', 'Binkowski', 'Giorgi', 'Peters']:
            for soueim in ['Slinn', 'Wiman', 'Giorgi', 'Peters', 'Pleim']:
                combine = soustokes + '&' + soueim
                for modesigma, modeinterval in sigmagbins.items():
                    weights = []
                    for di in data.index:
                        vds = []
                        for mi in arange(0, len(modeinterval)):
                            dpi, weighti = modeinterval[mi][0], modeinterval[mi][1]
                            vdi = 100 * weighti * diag(dens[di], dpi * 1e-6, temp[di], press[di], RH[di], ustar[di], Lo[di],
                                luc(wstar[di], ustar[di], Uh[di], z0[di], z[di], d[di], 'CMAQdiag').dforest(), 'CMAQdiag', soustokes, soueim).CMAQdiag()[0]
                            vds.append(vdi)
                        weights.append(sum(vds))
                    data['%s&%s' % (combine,modesigma)] = weights


        for modesigma, modeinterval in sigmagbins.items():
            weights = []
            for di in data.index:
                vds = []
                for mi in arange(0, len(modeinterval)):
                    dpi, weighti = modeinterval[mi][0], modeinterval[mi][1]
                    vdi = 100 * weighti * run(dens[di], dpi * 1e-6, temp[di], press[di], RH[di], ustar[di], Lo[di],
                                luc(wstar[di], ustar[di], Uh[di], z0[di], z[di], d[di], 'CMAQLAI').dforest(), 'CMAQLAI').CMAQLAI()[0]
                    vds.append(vdi)
                weights.append(sum(vds))
            data['Slinn&Slinn&LAI&%s' % (modesigma)] = weights
            #"""

        medians = {}
        for soustokes in ['Slinn', 'Binkowski', 'Giorgi', 'Peters']:
            for soueim in ['Slinn', 'Wiman', 'Giorgi', 'Peters', 'Pleim']:
                for sigma in ['1.01', '1.7', '2.5']:
                #for sigma in ['1.01']:
                    name = soustokes + '&' + soueim + '&' + sigma
                    medians[name] = getmedian(data, name)

        #add LAI one
        for sigma in ['1.01', '1.7', '2.5']:
            # for sigma in ['1.01']:
            name = 'Slinn&Slinn&LAI&' + sigma
            medians[name] = getmedian(data, name)


        Vd_obs = [0.18, 0.44, 0.95, 0.85, 0.28, 0.25] #cm/s
        timeinterval = [0.5, 1, 1.5, 2, 2.5, 3]

        plt.close()
        fig = plt.figure(figsize=(8, 4))
        ax = fig.add_axes([0.15, 0.12, 0.5, 0.8])
        markersize = 2

        diagls = ['-', ':', '-.', '--']
        diagcolors = ['blue', 'green', 'red', 'cyan', 'brown']
        markers = ['x', '+']
        for lsi, soustokes in enumerate(['Slinn', 'Binkowski', 'Giorgi', 'Peters']):
            for colori, soueim in enumerate(['Slinn', 'Wiman', 'Giorgi', 'Peters', 'Pleim']):
                for si, sigma in enumerate(['1.01']):
                    name = soustokes + '&' + soueim + '&' + sigma
                    combine = soustokes + '&' + soueim
                    if combine in ['Binkowski&Wiman']:
                        ax.plot(timeinterval, [i/4 for i in medians[name]], color=diagcolors[colori], lw=1, ls=diagls[lsi], marker=markers[si], markersize=markersize, mfc='white', label='0.25x ' + name.split('&')[0] + '&' + name.split('&')[1])

                    else:
                        ax.plot(timeinterval, medians[name], color=diagcolors[colori], lw = 1, ls = diagls[lsi], marker= markers[si], markersize=markersize, mfc='white', label = '1x '+ name.split('&')[0] + '&' + name.split('&')[1])

        locs, labels = plt.xticks()
        ax.set_xlim(0.3,3.2)
        ax.set_ylabel('Median Vd of each time interval (cm/s)')
        ax.set_xticklabels(['','2-6', '6-10', '10-14', '14-18', '18-22', '22-2',''])
        ax.set_xlabel('Time intervals of day (local time at central Japan)')
        ax.set_title('           Deciduous forest (%s = 1500 kg/m3, Dpg = 0.48 %sm)' % (r'$\rho$', r'$\mu$'),fontsize=9)
        ax.legend(ncol=1, bbox_to_anchor=(1.47,1),fontsize=8)
        fig.savefig('../Fig/Diag_imp_df_monodisperse.png', dpi=150)


