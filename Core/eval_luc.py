#land category defination and parameterization


class luc(object):
    def __init__(self, wstar, ustar, Uh, z0, z, d, model=''):
        self.model = model
        self.wstar = wstar
        self.ustar = ustar
        self.Uh = Uh
        self.z0 = z0
        self.z = z
        self.d = d

    def grass(self):
        #combinng all grass landuse type
        if self.model == 'Z01':
            localparams = dict(z0 = self.z0,
                               z = self.z,
                               d = self.d,
                               surface = 'vegetation',
                               gamma = 0.54,
                               alpha = 1.2,
                               beta = 2.,
                               A = 2.0/1000)
            return localparams

        elif self.model == 'ZS14':
            if min(self.d, 230/1000) == 230/1000:
                self.d = 200/1000
            #assume it is the same as plant
            localparams = dict(z0 = self.z0,
                               z= self.z,
                               d=self.d,
                               surface = 'vegetation',
                               hrough = 230/1000,
                               dc = 5 / 1000,
                               Uh = self.Uh,
                               fai = 0.4,
                               Ain = 150)
            return localparams

        elif self.model in ['CMAQ', 'CMAQstd', 'CMAQmode', 'CMAQimp']:
            localparams = dict(z0 = self.z0,
                               d=self.d,
                               z = self.z,
                               wstar = self.wstar)
            return localparams

        elif self.model in ['CMAQLAI', 'CMAQmodeLAI']:
            localparams = dict(z0 = self.z0,
                               d=self.d,
                               z = self.z,
                               wstar = self.wstar,
                               LAI = 4,
                               A=2.0 / 1000)
            return localparams


        elif self.model == 'CMAQdiag':
            localparams = dict(z0 = self.z0,
                               d = self.d,
                               z = self.z,
                               dc = 5/1000,
                               surface = 'vegetation',
                               alpha = 1.2,
                               beta =2.,
                               A = 2.0/1000,
                               wstar = self.wstar)
            return localparams

        else:
            print('no model is defined')

    def cforest(self):
        #evergreen, mix needle and broadleaf
        if self.model == 'Z01':
            localparams = dict(z0 = self.z0,
                               d=self.d,
                               z = self.z,
                               surface = 'vegetation',
                               gamma = (0.56+0.58)/2,
                               alpha = (1.0+0.6)/2,
                               beta = 2.,
                               A = (2+5)/2/1000) #(2.0 + 5.0)/2.
            return localparams


        elif self.model == 'ZS14':
            #assume it is the same as plant
            if min(self.d, 230/1000) == 230/1000:
                self.d = 200/1000
            localparams = dict(z0 = self.z0,
                               d=self.d,
                               z = self.z,
                               surface = 'vegetation',
                               hrough = 230/1000,
                               dc = 5 / 1000,
                               Uh = self.Uh,
                               fai = 0.4,
                               Ain = 150)
            return localparams

        elif self.model in ['CMAQ', 'CMAQstd', 'CMAQmode', 'CMAQimp']:
            localparams = dict(z0 = self.z0,
                               d=self.d,
                               z = self.z,
                               wstar = self.wstar)
            return localparams

        elif self.model in ['CMAQLAI', 'CMAQmodeLAI']:
            localparams = dict(z0 = self.z0,
                               d=self.d,
                               z = self.z,
                               wstar = self.wstar,
                               LAI = 6,
                               A=(2 + 5) / 2 / 1000)
            return localparams

        elif self.model == 'CMAQdiag':
            localparams = dict(z0 = self.z0,
                               d = self.d,
                               z = self.z,
                               dc = 5/1000,
                               surface = 'vegetation',
                               alpha = (1.0+0.6)/2,
                               beta =2.,
                               A =  (2+5)/2/1000,
                               wstar = self.wstar)
            return localparams

        else:
            print('no model is defined')

    def dforest(self):
        #deciduous forest
        #because CMAQ landuse does not classify needleleaf and broadleaf, so params are taken averaged value from
        #both landuse types
        #needle + broad / 2.
        if self.model == 'Z01':
            localparams = dict(z0 = self.z0,
                               d=self.d,
                               z =  self.z,
                               surface = 'vegetation',
                               gamma = 0.56,
                               alpha = (1.1 + .8)/2.,
                               beta = 2.,
                               A = 3.5/1000) #(2.0 + 5.0)/2.
            return localparams

        elif self.model == 'ZS14':
            if min(self.d, 230/1000) == 230/1000:
                self.d = 200/1000
            localparams = dict(z0 = self.z0,
                               d=self.d,
                               z = self.z,
                               surface = 'vegetation',
                               hrough = 230/1000,
                               dc = 5 / 1000,
                               Uh = self.Uh,
                               fai = 0.4,
                               Ain = 150)
            return localparams

        elif self.model in ['CMAQ', 'CMAQstd', 'CMAQmode', 'CMAQimp']:
            localparams = dict(z0 = self.z0,
                               d=self.d,
                               z = self.z,
                               wstar = self.wstar)
            return localparams

        elif self.model in ['CMAQLAI', 'CMAQmodeLAI']:
            localparams = dict(z0 = self.z0,
                               d=self.d,
                               z = self.z,
                               wstar = self.wstar,
                               LAI = 6,
                               A=3.5 / 1000)
            return localparams

        elif self.model == 'CMAQdiag':
            localparams = dict(z0 = self.z0,
                               d = self.d,
                               z = self.z,
                               dc = 5/1000,
                               surface = 'vegetation',
                               alpha = (1.1 + .8)/2.,
                               beta =2.,
                               A =  3.5/1000,
                               wstar = self.wstar)
            return localparams

        else:
            print('no model is defined')

    def water(self):
        #waterface, unavailable parameters will use grass as surrogate
        if self.ustar <= 0.16:
            z0 = 0.021 * self.ustar ** 3.32
        else:
            z0 = 0.00098 * self.ustar ** 1.65

        if self.model == 'Z01':
            localparams = dict(z0 = z0,
                               d=0,
                               z = self.z,
                               surface = 'smooth',
                               gamma = 0.50,
                               alpha = 100.,
                               beta = 2.,
                               A = 2.0/1000)
            return localparams

        elif self.model == 'ZS14':
            localparams = dict(z0 = z0,
                               d = 0,
                               z= self.z,
                               surface = 'smooth',
                               hrough = 30 * z0,
                               dc = 0.1 / 1000,
                               Uh = self.Uh,
                               fai = 0.538,
                               Ain = 100)
            return localparams

        elif self.model in ['CMAQ', 'CMAQstd', 'CMAQmode', 'CMAQimp']:
            localparams = dict(z0 = z0,
                               d = 0,
                               z = self.z,
                               wstar = self.wstar)
            return localparams

        elif self.model in ['CMAQLAI', 'CMAQmodeLAI']:
            localparams = dict(z0 = self.z0,
                               d=self.d,
                               z = self.z,
                               wstar = self.wstar,
                               LAI = 'N/A')
            return localparams

        elif self.model == 'CMAQdiag':
            localparams = dict(z0 = z0,
                               d = 0,
                               z = self.z,
                               dc = 0.1/1000,
                               surface = 'smooth',
                               alpha = 100.,
                               beta =2.,
                               A =  2.0/1000,
                               wstar = self.wstar)
            return localparams

        else:
            print('no model is defined')
