#This is core function to include all tested deposition algorithms.
#Details can be found in Shu et al., 2021 (GMD).

from numpy import *

defaultkey = dict(
    # (values, units, physical meaning)
    # constant
    pi=(3.1416, '/', 'constant'),
    g=(9.8, 'm/s2', 'gravity accelaration'),
    kvon=(0.4, '/', 'Von Karman constant'),
    kbo=(1.38e-23, 'J/K', 'Boltzmann constant'),


    # varaibles
    dens=(1500, 'Kg/m3', 'particle density'),
    dim=(2.5e-6, 'm', 'particle diameter'),
    temp=(298, 'K', 'temperature'),
    RH=(0.9, '/', 'relative humidity'),
    press=(101325, 'Pasca', 'pressure'),
    ustar=(0.5, 'm/s', 'friction velocity'),

    # calculate using T later
    kvair=(1.57e-5, 'm2/s', 'kinematic viscosity of air, as function of Temp'),
    dvair=(1.89e-5, 'kg/m/s', 'Dynamic viscosity of air, for Vg, as function of Temp'),
    mfpath=(0.067e-6, 'm', 'air molecules mean free path'), #PZ10

    # refine in Landuse.py
    z=(0, 'm', 'measurement height'),
    Lo=(50, 'm', 'Monin-Obukhov length'),
    zr=(1, 'm', 'height at which the Vd is evaluated'),
    z0=(0.03, 'm', 'roughness length'),
    surface=([], '', 'vegetation or else')  # for different models, use different keywords
    )


class __models(object):
    def __init__(self, dens, dim, temp, press, RH, ustar, Lo, par=dict()):
        # variables
        self.dens = dens
        self.dim = dim
        self.temp = temp
        self.press = press
        self.RH = RH
        self.ustar = ustar
        self.Lo = Lo
        self.par = par


        # constant
        self.g = 9.8
        self.kvon = 0.4
        self.kbo = 1.38e-23

        # ACP page 909, eq(19.18)
        self.dvair = 1.8e-5 * (self.temp / 298) ** 0.85
        #self.dvair2 = 1.458e-6 * (self.temp) * sqrt(self.temp) / (self.temp + 110.4)
        #print(self.dvair, self.dvair2)

        airdens = 1.225  # kg/m3
        self.kvair = self.dvair / airdens
        # ACP page 401, assume z= Ma/Mb = 1
        self.mfpath = self.kbo * self.temp / (sqrt(2) * pi * self.press * self.dim ** 2)
        self.C = 1 + 2 * (self.mfpath / self.dim) * (1.257 + 0.4 * exp((-0.55) * self.dim / self.mfpath))
        self.Vg = (self.dens * (self.dim) ** 2 * self.g) / (18 * self.dvair) * self.C


    def CMAQ(self):
        # phiH is the stability function for heat.
        self.zr = self.par['z'] - self.par['d']

        if (1/self.Lo) < 0:
            phiH = 2 * log((sqrt(1 - 16 * self.zr / self.Lo)+1)/(sqrt(1-16*self.par['z0']/self.Lo)+1))
        elif (self.zr - self.par['z0'])/self.Lo <=1:
            phiH = -5 * (self.zr - self.par['z0']) / self.Lo
        else:
            phiH = 1 - 5 - (self.zr - self.par['z0']) / self.Lo
        Ra = 0.95 * (log(self.zr / self.par['z0']) - phiH) / (self.kvon * self.ustar)
        D = self.C * self.kbo * self.temp / (3 * pi * self.dvair * self.dim)
        Sc = self.kvair / D
        St = self.Vg * self.ustar ** 2 / (self.g * self.kvair)
        c = 400
        EIM = St ** 2 / (c + St ** 2)
        EIN = 0          # EIN is not used in CMAQ
        convfact = 1 + 0.24 * (self.par['wstar'] / self.ustar) ** 2
        EB = Sc**(-2/3.)
        Rb = 1 / (convfact * self.ustar * (EB + EIM + EIN))
        fustar = (convfact * self.ustar)
        Vd = self.Vg / (1 - exp(-self.Vg * (Ra + Rb)))
        return [Vd, self.Vg, Ra, Rb, fustar, EB, EIM, EIN]


    def CMAQstd(self):
        # phiH is the stability function for heat.
        self.zr = self.par['z'] - self.par['d']

        if (1/self.Lo) < 0:
            phiH = 2 * log((sqrt(1 - 16 * self.zr / self.Lo)+1)/(sqrt(1-16*self.par['z0']/self.Lo)+1))
        elif (self.zr - self.par['z0'])/self.Lo <=1:
            phiH = -5 * (self.zr - self.par['z0']) / self.Lo
        else:
            phiH = 1 - 5 - (self.zr - self.par['z0']) / self.Lo

        Ra = 0.95 * (log(self.zr / self.par['z0']) - phiH) / (self.kvon * self.ustar)
        D = self.C * self.kbo * self.temp / (3 * pi * self.dvair * self.dim)
        Sc = self.kvair / D
        St = self.Vg * self.ustar ** 2 / (self.g * self.kvair)
        c = 400
        EIM = St ** 2 / (c + St ** 2)
        EIN = 0          # EIN is not used in CMAQ
        convfact = 1 + 0.24 * (self.par['wstar'] / self.ustar) ** 2
        EB = Sc ** -(2/3.)
        Rb = 1 / (convfact * self.ustar * (EB + EIM + EIN))
        fustar = (convfact * self.ustar)
        Vd = self.Vg + 1 / (Ra + Rb)
        print(Ra,Rb)
        return [Vd, self.Vg, Ra, Rb, fustar, EB, EIM, EIN, 1/(Ra+Rb)]

    def CMAQmode(self):
        # phiH is the stability function for heat.
        self.zr = self.par['z'] - self.par['d']

        if (1/self.Lo) < 0:
            phiH = 2 * log((sqrt(1 - 16 * self.zr / self.Lo)+1)/(sqrt(1-16*self.par['z0']/self.Lo)+1))
        elif (self.zr - self.par['z0'])/self.Lo <=1:
            phiH = -5 * (self.zr - self.par['z0']) / self.Lo
        else:
            phiH = 1 - 5 - (self.zr - self.par['z0']) / self.Lo
        Ra = 0.95 * (log(self.zr / self.par['z0']) - phiH) / (self.kvon * self.ustar)
        EIN = 0
        convfact = 1 + 0.24 * (self.par['wstar'] / self.ustar) ** 2
        fustar = (convfact * self.ustar)

        #mode options
        logsigmag2 = log(self.sigmag) ** 2

        #moment
        k = 0 # for grass and coniferous forest because they both use number concentration
        #k = 3 # for deciduous forst because they use mass concentration
        if self.modeVg:
            dimg = self.dim
            Vgmean = (dimg) ** 2 * self.dens * self.g / (18 * self.dvair)
            Vghat = Vgmean *  (exp((4*k+4)/2. * logsigmag2) + 1.246 * (2*self.mfpath/dimg) * exp((2*k+1)/2. * logsigmag2))
            myVg = Vghat
        else:
            myVg = self.Vg


        if self.modeEB:
            dimg = self.dim
            Dmean = self.kbo * self.temp / (3 * pi * self.dvair * dimg)
            Dhat = Dmean * (exp((-2*k+1)/2. * logsigmag2) + 1.246 * (2*self.mfpath/dimg) * exp((-4*k+4)/2. * logsigmag2))
            myD = Dhat
            Sc = self.kvair / myD
            EB = Sc**(-2/3.)
        else:
            D = self.kbo * self.temp / (3 * pi * self.dvair * self.dim) * self.C
            Sc = self.kvair / D
            EB = Sc**(-2/3.)

        if self.modeEIM:
            Vgmean = (self.dim) ** 2 * self.dens * self.g / (18 * self.dvair)
            St = Vgmean * self.ustar ** 2 / (self.g * self.kvair)
            if k == 0:
                EIM = min(St **2 /400 * exp(8*logsigmag2), 1.0)
            if k == 3:
                EIM = min(St **2 /400 * exp(20*logsigmag2), 1.0)
        else:
            St = self.ustar ** 2 / (self.g * self.kvair) * self.Vg
            EIM = St ** 2 / (400 + St ** 2)

        Rb = 1 / (convfact * self.ustar * (EB + EIM + EIN))
        Vd = myVg / (1 - exp(-myVg * (Ra + Rb)))
        return [Vd, myVg, Ra, Rb, EB, EIM, EIN]

    def CMAQmodeLAI(self):
        # phiH is the stability function for heat.
        self.zr = self.par['z'] - self.par['d']

        if (1/self.Lo) < 0:
            phiH = 2 * log((sqrt(1 - 16 * self.zr / self.Lo)+1)/(sqrt(1-16*self.par['z0']/self.Lo)+1))
        elif (self.zr - self.par['z0'])/self.Lo <=1:
            phiH = -5 * (self.zr - self.par['z0']) / self.Lo
        else:
            phiH = 1 - 5 - (self.zr - self.par['z0']) / self.Lo
        Ra = 0.95 * (log(self.zr / self.par['z0']) - phiH) / (self.kvon * self.ustar)
        EIN = 0
        fveg = 1 # for single luc
        convfact = (1 + fveg * max(0,self.par['LAI']-1))  + 0.24 * (self.par['wstar'] / self.ustar) ** 2
        fustar = (convfact * self.ustar)

        #mode options
        logsigmag2 = log(self.sigmag) ** 2

        #moment
        k = 0 # for grass and coniferous forest because they both use number concentration
        #k = 3 # for deciduous forst because they use mass concentration

        if self.modeVg:
            dimg = self.dim
            Vgmean = (dimg) ** 2 * self.dens * self.g / (18 * self.dvair)
            Vghat = Vgmean *  (exp((4*k+4)/2. * logsigmag2) + 1.246 * (2*self.mfpath/dimg) * exp((2*k+1)/2. * logsigmag2))
            myVg = Vghat
        else:
            myVg = self.Vg


        if self.modeEB:
            dimg = self.dim
            Dmean = self.kbo * self.temp / (3 * pi * self.dvair * dimg)
            Dhat = Dmean * (exp((-2*k+1)/2. * logsigmag2) + 1.246 * (2*self.mfpath/dimg) * exp((-4*k+4)/2. * logsigmag2))
            myD = Dhat
            Sc = self.kvair / myD
            EB = Sc**(-2/3.)
        else:
            D = self.C * self.kbo * self.temp / (3 * pi * self.dvair * self.dim)
            Sc = self.kvair / D
            EB = Sc**(-2/3.)

        St = myVg * self.ustar / (self.g * self.par['A'])
        EIM = St ** 2 / (1 + St **2)
        Rb = 1 / (convfact * self.ustar * (EB + EIM + EIN))
        Vd = myVg / (1 - exp(-myVg * (Ra + Rb)))
        return [Vd, myVg, Ra, Rb, EB, EIM, EIN]


    def CMAQimp(self):
        # phiH is the stability function for heat.
        self.zr = self.par['z'] - self.par['d']

        if (1/self.Lo) < 0:
            phiH = 2 * log((sqrt(1 - 16 * self.zr / self.Lo)+1)/(sqrt(1-16*self.par['z0']/self.Lo)+1))
        elif (self.zr - self.par['z0'])/self.Lo <=1:
            phiH = -5 * (self.zr - self.par['z0']) / self.Lo
        else:
            phiH = 1 - 5 - (self.zr - self.par['z0']) / self.Lo
        Ra = 0.95 * (log(self.zr / self.par['z0']) - phiH) / (self.kvon * self.ustar)
        D = self.C * self.kbo * self.temp / (3 * pi * self.dvair * self.dim)
        Sc = self.kvair / D
        St = self.Vg * self.ustar ** 2 / (self.g * self.kvair)
        c = self.impc
        EIM = St ** 2 / (c + St ** 2)
        EIN = 0          # EIN is not used in CMAQ
        convfact = 1 + 0.24 * (self.par['wstar'] / self.ustar) ** 2
        EB = Sc**(-2/3.)
        Rb = 1 / (convfact * self.ustar * (EB + EIM + EIN))
        #print(Rb)
        fustar = (convfact * self.ustar)
        Vd = self.Vg / (1 - exp(-self.Vg * (Ra + Rb)))
        #print(Ra, Rb, self.Vg, Vd)
        return [Vd, self.Vg, Ra, Rb, fustar, EB, EIM, EIN]

    def CMAQLAI(self):
        # phiH is the stability function for heat.
        self.zr = self.par['z'] - self.par['d']

        if (1/self.Lo) < 0:
            phiH = 2 * log((sqrt(1 - 16 * self.zr / self.Lo)+1)/(sqrt(1-16*self.par['z0']/self.Lo)+1))
        elif (self.zr - self.par['z0'])/self.Lo <=1:
            phiH = -5 * (self.zr - self.par['z0']) / self.Lo
        else:
            phiH = 1 - 5 - (self.zr - self.par['z0']) / self.Lo
        Ra = 0.95 * (log(self.zr / self.par['z0']) - phiH) / (self.kvon * self.ustar)
        D = self.C * self.kbo * self.temp / (3 * pi * self.dvair * self.dim)
        Sc = self.kvair / D
        St = self.Vg * self.ustar / (self.par['A'] * self.g)
        c = 1
        EIM = St ** 2 / (c + St ** 2)
        EIN = 0          # EIN is not used in CMAQ
        fveg = 1         # for single luc
        convfact = (1 + fveg * max(0,self.par['LAI']-1)) + 0.24 * (self.par['wstar'] / self.ustar) ** 2
        EB = Sc**(-2/3.)
        Rb = 1 / (convfact * self.ustar * (EB + EIM + EIN))
        fustar = (convfact * self.ustar)
        Vd = self.Vg / (1 - exp(-self.Vg * (Ra + Rb)))
        return [Vd, self.Vg, Ra, Rb, fustar, EB, EIM, EIN]


    def CMAQdiag(self):
        # phiH is the stability function for heat.
        self.zr = self.par['z'] - self.par['d']

        if (1/self.Lo) < 0:
            phiH = 2 * log((sqrt(1 - 16 * self.zr / self.Lo)+1)/(sqrt(1-16*self.par['z0']/self.Lo)+1))
        elif (self.zr - self.par['z0'])/self.Lo <=1:
            phiH = -5 * (self.zr - self.par['z0']) / self.Lo
        else:
            phiH = 1 - 5 - (self.zr - self.par['z0']) / self.Lo
        Ra = 0.95 * (log(self.zr / self.par['z0']) - phiH) / (self.kvon * self.ustar)

        D = self.C * self.kbo * self.temp / (3 * pi * self.dvair * self.dim)
        Sc = self.kvair / D

        #CMAQdiag is written for testing different combination of stokes number and impaction
        # taup:  particle relaxation time
        taup = self.Vg / self.g
        if self.soustokes == 'Slinn': #Z01
            St = taup * self.ustar / self.par['A']
        elif self.soustokes == 'Binkowski': # Old CMAQ
            St = self.Vg * self.ustar ** 2 / (self.g * self.kvair)
        elif self.soustokes == 'Giorgi':
            St = taup * self.ustar ** 2 / self.par['dc']
        elif self.soustokes == 'Peters':
            St = self.dens * self.dim ** 2 / (9 * self.dvair * self.par['dc'])
        else:
            print('error, no stokes selected')

        if self.soueim == 'Slinn':
            EIM = St ** 2 / (1 + St ** 2)
        elif self.soueim == 'Wiman':
            EIM = St / (3 + St)
        elif self.soueim == 'Giorgi':
            EIM = (St / (0.6 + St)) ** 3.2
        elif self.soueim == 'Peters':
            EIM = (St / (0.8 + St)) ** 2
        elif self.soueim == 'Pleim':
            EIM = St**2 / (400 + St **2)
        elif self.soueim == 'Z01':
            if self.par['surface'] == 'smooth':
                EIM = 10 ** (-3 / St)
            elif self.par['surface'] == 'vegetation':
                EIM = (St / (self.par['alpha'] + St)) ** (self.par['beta'])
        else:
            print('error, no eim selected')

        EIN = 0
        convfact = 1 + 0.24 * (self.par['wstar'] / self.ustar) ** 2
        EB = Sc ** -(2/3.)
        Rb = 1 / (convfact * self.ustar * (EB + EIM + EIN))
        fustar = (convfact * self.ustar)
        Vd = self.Vg / (1 - exp(-self.Vg * (Ra + Rb)))

        return [Vd, self.Vg, Ra, Rb, fustar, EB, EIM, EIN]


    def Z01(self):
        # phiH is the stability function for heat.
        self.zr = self.par['z'] - self.par['d']
        x = self.par['z'] / self.Lo

        if (x <= 0):
            phiH = 2 * log(0.5 * (1 + (1 - 16 * x) ** (0.5)))
        elif (x > 0):
            phiH = -5 * x
        else:
            print('Error, x is not in the range')

        Ra = (log(self.zr / self.par['z0']) - phiH) / (self.kvon * self.ustar)
        D = self.C * self.kbo * self.temp / (3 * pi * self.dvair * self.dim)
        Sc = self.kvair / D
        EB = Sc ** (-self.par['gamma'])

        if self.par['surface'] == 'smooth':
            St = self.Vg * (self.ustar) ** 2 / self.kvair
            EIM = 10 ** (-3 / St)
        elif self.par['surface'] == 'vegetation':
            St = self.Vg * self.ustar / (self.g * self.par['A'])
            EIM = (St / (self.par['alpha'] + St)) ** (self.par['beta'])
        else:
            print('not valid surface key word')

        EIN = 0.5 * (self.dim / self.par['A']) ** 2
        R1 = exp(-(St ** 0.5))
        e0 = 3  # empirical constant
        Rs = 1 / (e0 * self.ustar * (EB + EIM + EIN) * R1)
        fustar = (e0 * R1 * self.ustar)
        Vd = self.Vg + 1 / (Ra + Rs)
        return [Vd, self.Vg, Ra, Rs, fustar, EB, EIM, EIN, 1/(Ra+Rs)]

    def PZ10(self):
        Vdrift = self.Vg + self.par['Vphor']
        self.zr = self.par['z'] - self.par['d']
        x = self.par['z'] / self.Lo

        if (x <= 0):
            phiH = 2 * log(0.5 * (1 + (1 - 16 * x) ** (0.5)))
        elif (x > 0): #and (x <= 1):
            phiH = -5 * x
        else:
            print('Error, x is not in the range')
        Ra = (log(self.zr / self.par['z0']) - phiH) / (self.kvon * self.ustar)

        # Q, Qg, alpha, eta are dependent on the aerodynamic and surface characteristic features
        # Brownian diffusion
        D = self.C * self.kbo * self.temp / (3 * pi * self.dvair * self.dim)
        Sc = self.kvair / D
        F = (Sc ** (1 / 3.)) / 2.9
        Egb = (Sc ** (-2 / 3.) / 14.5) / (1 / 6. * log((1 + F) ** 2 / (1 - F + F ** 2)) + 1 / sqrt(3) * arctan(
            (2 * F - 1) / sqrt(3)) + pi / (6 * sqrt(3)))

        if self.par['surface'] == 'vegetation':
            # Turbulent impaction
            # phim is the non-dimensional stability function of momentum.
            if (x <= 0):
                phim = (1 - 16 * x) ** (-1 / 4.)
            elif (x > 0): #and (x <= 1):
                phim = 1 + 5 * x
            else:
                print('Error, x is not in the range')

            # faext: aerodynamic extinction coefficient
            faext = (self.par['fcin'] * self.par['LAI'] / (
                    12 * self.kvon ** 2 * (1 - self.par['d'] / self.par['h']) ** 2)) ** (1 / 3.) * \
                    phim ** (2 / 3.) * ((self.par['h'] - self.par['d']) / abs(self.Lo))
            # uf: local friction velocity
            uf = self.ustar * exp(-faext)
            # taup:  particle relaxation time
            taup = self.Vg / self.g
            tauph1 = taup * uf ** 2 / self.kvair
            Egt = 2.5e-3 * self.par['CIT'] * tauph1
            Eg = Egb + Egt

            # phiHlmp is the stability function for heat, calculation is different from that for Ra
            if (x <= 0):
                phiHlmp = (1 - 16 * x) ** (-1 / 2.)
            elif (x > 0): #and (x <= 1):
                phiHlmp = 1 + 5 * x
            else:
                print('Error, x is not in the range')
            # mixing height for particle
            lmph = self.kvon * (self.par['h'] - self.par['d']) / phiHlmp / (
                    (self.par['h'] - self.par['d']) / abs(self.Lo)) #the same as khan code, use abs

            # Reynolds number of the horizontal air flow calculated at top of canopy height h
            Reh = self.par['Uh'] * self.par['L'] / self.kvair
            EB = self.par['fcb'] * Sc ** (-2 / 3.) * Reh ** (-1 / 2.)

            #only for vegetation, if smooth, no need, when needle and leaf mix, how to deal with it
            if self.par['obstacle'] == 'needle':
                EIN = self.par['fcb'] * self.dim / self.par['L']
                #print(EIN)
            elif self.par['obstacle'] == 'leaf':
                EIN = self.par['fcb'] * self.dim / self.par['L'] * (
                        2 + log(4 * self.par['L'] / self.dim))
                #print(EIN)
            else:
                print('Error, cannot identify obstacle')

            # Sth, stokes number on top of the canopy
            Sth = taup * self.par['Uh'] / self.par['L']
            # CIM and betaIM are LUC-dependent coefficients
            EIM = self.par['CIM'] * (Sth / (Sth + self.par['betaIM'])) ** 2

            # tauph2 is different from tauph1
            tauph2 = taup * self.ustar ** 2 / self.kvair
            if (tauph2 <= 20):
                EIT = 2.5e-3 * self.par['CIT'] * tauph2
            else:
                EIT = self.par['CIT']

            # ET, total particle collection efficiency
            ET = self.par['Uh'] / self.ustar * (EB + EIN + EIM) + EIT

            # Q is the ratio the turbulent transport timescale to the vegetation collection timescale
            # Q << 1, homogeneous particle prevails (Aitken and accumulation mode)
            # Q >> 1, inhomogeneous particle prevails (coarse)
            Q = self.par['LAI'] * ET * self.par['h'] / lmph
            Qg = Eg * self.par['h'] / lmph
            eta = sqrt(faext ** 2 / 4. + Q)
            Vds = Eg * self.ustar * ((1 + (Q / Qg - faext / 2.) * tanh(eta) / eta)) / ((1 + (Q / Qg + faext/ 2.) * tanh(eta) / eta))
            Rs = 1 / Vds
            Vd = Vdrift + 1 / (Ra + Rs)
        else:
            Vd = Vdrift + 1 / (Ra + 1 / (Egb * self.ustar))
        return [Vd, self.Vg]

    def ZS14(self):
        Rg = 1 / self.Vg
        self.zr = self.par['z'] - self.par['d']
        # ScT, turbulent Schmidt number
        ScT = (1 + (self.Vg ** 2 / self.ustar ** 2)) ** (0.5)
        if self.par['surface'] == 'vegetation':
            Ra = ScT / (self.kvon * self.ustar) * \
                 log((self.par['z'] - self.par['d']) / (self.par['hrough'] - self.par['d']))
            Vdm = self.ustar / (self.par['Uh'] * self.par['hrough'])
        elif self.par['surface'] == 'smooth':
            Ra = 0.45 * ScT / (self.kvon * self.ustar) * log(self.par['z'] / self.par['z0'])
            Vdm = 3 * self.ustar
        else:
            print('Other surface option is not suitable for ZS14')

        # dc is the diameter of the surface collection element and the value is given by zhang and shao (2014)
        taup = self.Vg / self.g
        St = taup * self.ustar / self.par['dc']
        D = self.C * self.kbo * self.temp / (3 * pi * self.dvair * self.dim)
        Sc = self.kvair / D

        # CB and nB can be got from Zhang and Shao (2014)
        # Reynolds number of the horizontal air flow calculated at top of canopy height h
        # Reh = self.par['Uh'] * self.par['L'] / self.kvair
        Reh = self.par['Uh'] * self.par['dc'] / self.kvair

        if Reh <= 4e3:
            fcb = 0.467
            nB = 1 / 2.
        elif Reh > 4e3 and Reh <= 4e4:
            fcb = 0.203
            nB = 3 / 5.
        elif Reh > 4e4 and Reh <= 4e5:
            fcb = 0.025
            nB = 4 / 5.
        else:
            print('Reh is out of range')
        EB = fcb * Sc ** (-2 / 3.) * Reh ** (nB - 1)
        EIM = (St / (0.6 + St)) ** 2
        # wet diameter after humidity correction
        # rad = self.dim / 2.
        # radw = self.par['C1'] * rad ** self.par['C2'] / (self.par['C3'] * rad ** self.par['C4'] - log10(self.par['RH'])) + rad ** 3
        # dimw =  radw * 2
        # EIN = self.par['Ain'] * self.ustar * 10 ** (-St) * 2 * dimw / self.par['dc']

        EIN = self.par['Ain'] * self.ustar * 10 ** (-St) * 2 * self.dim / self.par['dc']
        E = EB + EIM + EIN

        # efai, effective frontal area index
        # rtau, ratio of the drag on the roof of the roughness element to the total shear stress
        # beta is the ratio of the pressure-drag coefficient to frictio-drag coefficient
        pai = self.par['fai'] / (self.par['hrough'] * self.par['dc']) * ((3.1416 * self.par['dc'] ** 2) / 4)
        efai = self.par['fai'] / (1 - pai) ** 0.1 * exp(-6 * self.par['fai'] / (1 - pai) ** 0.1)
        beta = 200  # Zhang and Shao (2014)
        rtau = beta * efai / (1 + beta * efai)

        # Vgw is Vg after humidity correction, using wet dim to calculate wet Vg( Vgw)
        # Vgw = (dens * dimw ** 2 * self.g * self.C) / (18 * self.dvair)
        # use dry dim for intercomparison, because other models are dry
        Vgw = (self.dens * self.dim ** 2 * self.g * self.C) / (18 * self.dvair)

        taupw = Vgw / self.g
        # near the surface
        taupwsurf = taupw * self.ustar ** 2 / self.kvair

        # Cd is the drag partition coefficient
        Cd = 0.3  # (Shao and Yang, 2008)
        # b is an empirical constant
        #b = 1  # Giorgi(1988) and Zhang et al. (2001)
        b = 0.01 #for plant from ZS14 paper
        # b = 2, #Chamberlain (1967)
        R = exp(-b * sqrt(St))
        Rs = 1 / (R * Vdm * (E / Cd * rtau + (1 + rtau) / Sc + 10 ** (-3 / taupwsurf)) + Vgw)
        Vd = 1 / (Rg + ((Rs - Rg) / exp(Ra / Rg)))
        return [Vd, self.Vg, Ra, Rs]

    """
    def ZH14(self):


     #   for bulk aerosol, not good to compare
     #   ZH14 local params for long grass, only for
     #   LAI=([4],'','change with LAI')
     #   a1=(4.8e-3,'','')
     #   b1=(-7.9e-2,'','')
     #   b2=(1,'','')
     #   b3=(6.6e-1,'','')
     #   c1=(5.1,'','')
     #   c2=(-4.2,'','')
     #   c3=(9.9e-1,'','')
     #   d1=(-2.0,'','')
     #   d2=(6.3e1,'','')
     #   d3=(-1.6e1,'','')
     #   f1=(1.1e1,'','')
     #   f2=(-2.0e1,'','')
     #   f3=(1.1e1,'','')


        # phiH is the stability function for heat.
        x = self.par['z'] / self.Lo
        if (x >= -2) and (x <= 0):
            phiH = 2 * log(0.5 * (1 + (1 - 16 * x) ** (0.5)))
        elif (x > 0) and (x <= 1):
            phiH = -5 * x
        else:
            print('Error, x is not in the range')
        Ra = (log(self.zr / self.par['z0']) - phiH) / (self.kvon * self.ustar)

        # overall structure is similar to Z01 except modifying Rs
        # Rs were modified for three bulk particle sizes (i.e., PM2.5, PM2.5_10, and PM10+)
        if (self.dim <= 2.5e-6):  # unit is meter
            # a1 is an empirical constant derived by regression analysis (from Zhang and He 2014)
            Vds = self.par['a1'] * self.ustar
        elif (self.dim > 2.5e-6) and (self.dim <= 10e-6):
            k1 = self.par['c1'] * self.ustar + \
                 self.par['c2'] * self.ustar ** 2 * self.par['c3'] * self.ustar ** 3
            # LAI need to fix
            Vds = (self.par['b1'] * self.ustar + self.par['b2'] * self.ustar ** 2 * self.par['b3'] * self.ustar ** 3) \
                  * exp(k1 * (self.par['LAI'][0] / max(self.par['LAI']) - 1))
        elif (self.dim > 10e-6):
            k2 = self.par['f1'] * self.ustar + self.par['f2'] * self.ustar ** 2 + self.par['f3'] * self.ustar ** 3
            Vds = (self.par['d1'] * self.ustar + self.par['d2'] * self.ustar ** 2 + self.par['d3'] * self.ustar ** 3) \
                  * exp(k2 * (self.par['LAI'][0] / max(self.par['LAI']) - 1))
        else:
            print('Error, diameter is not in the range')
        Rs = 1 / Vds
        Vd = self.Vg + 1 / (Ra + Rs)
        return [Vd, self.Vg]
    """


class run(__models):
    def __init__(self, dens, dim, temp, press, RH, ustar, Lo, par, model=''):
        """
        model options: Z01 PZ10 ZH14 ZS14
        more models can be added later
        """
        super().__init__(dens=dens, dim=dim, temp=temp, press=press, RH=RH, ustar=ustar, Lo=Lo, par=par)
        self.model = model
        if model == 'Z01':
            Vd = self.Z01()[0]
            Vg = self.Z01()[1]

        elif model == 'PZ10':
            Vd = self.PZ10()[0]
            Vg = self.PZ10()[1]

        elif model == 'ZH14':
            Vd = self.ZH14()[0]
            Vg = self.ZH14()[1]

        elif model == 'ZS14':
            Vd = self.ZS14()[0]
            Vg = self.ZS14()[1]

        elif model == 'CMAQ':
            Vd = self.CMAQ()[0]
            Vg = self.CMAQ()[1]

        elif model == 'CMAQstd':
            Vd = self.CMAQstd()[0]
            Vg = self.CMAQstd()[1]

        elif model == 'CMAQLAI':
            Vd = self.CMAQLAI()[0]
            Vg = self.CMAQLAI()[1]

        else:
            print('No model option has been selected')

class diag(__models):
    def __init__(self, dens, dim, temp, press, RH, ustar, Lo, par, model='', soustokes='', soueim=''):
        super().__init__(dens=dens, dim=dim, temp=temp, press=press, RH=RH, ustar=ustar, Lo=Lo, par=par)
        self.model = model
        self.soustokes = soustokes
        self.soueim = soueim

class mode(__models):
    def __init__(self, dens, dim, temp, press, RH, ustar, Lo, par, sigmag, model='', modeVg=False, modeEB=False, modeEIM=False):
        super().__init__(dens=dens, dim=dim, temp=temp, press=press, RH=RH, ustar=ustar, Lo=Lo, par=par)
        self.model = model
        self.sigmag = sigmag
        self.modeVg = modeVg
        self.modeEB = modeEB
        self.modeEIM = modeEIM

class laimode(__models):
    def __init__(self, dens, dim, temp, press, RH, ustar, Lo, par, sigmag, model='', modeVg=False, modeEB=False):
        super().__init__(dens=dens, dim=dim, temp=temp, press=press, RH=RH, ustar=ustar, Lo=Lo, par=par)
        self.model = model
        self.sigmag = sigmag
        self.modeVg = modeVg
        self.modeEB = modeEB

class imp(__models):
    def __init__(self, dens, dim, temp, press, RH, ustar, Lo, par, impc, model=''):
        super().__init__(dens=dens, dim=dim, temp=temp, press=press, RH=RH, ustar=ustar, Lo=Lo, par=par)
        self.model = model
        self.impc = impc
