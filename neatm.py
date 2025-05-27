import astropy.units as u
import astropy.constants as const
import numpy as np
from scipy.integrate import dblquad
from uncertainties import ufloat
from uncertainties.umath import *
from astroquery.jplhorizons import Horizons
from datetime import datetime, timedelta

class neatm(object):
    """
    Calculate properties associated with a small body using the NEATM (Delbo & Harris 2002)
    Eqns 19 - 21
    """

    def __init__(self,objectName,obsDate,obsTime,nu,pv,phase_int,bolo_emissivity,rad_emissivity,eta,uncertainties=True):

        self.objectName = objectName
        self.obsDate = obsDate
        self.obsTime = obsTime
        self.nu = nu
        self.pv = pv
        self.phase_int = phase_int
        self.bolo_emissivity = bolo_emissivity
        self.rad_emissivity = rad_emissivity
        self.eta = eta
        self.uncertainties = uncertainties
        self.lamda = self.nu.to(u.m,equivalencies=u.spectral())

        self.get_eph()
        self.subsolar_temp()

    def get_eph(self):
        #Find the time one minute ahead so we can generate a JPL/HORIZONS ephemeris

        dObs = datetime.strptime(self.obsDate+' '+self.obsTime, '%Y-%m-%d %H:%M:%S')
        deltaObs = dObs + timedelta(minutes=1)
        startTime = dObs.strftime('%Y-%m-%d %H:%M')
        endTime = deltaObs.strftime('%Y-%m-%d %H:%M')

        #Run the query to JPL/HORIZONS and store the results in a dataframe
        obj = Horizons(id=self.objectName, id_type='designation', location='ALMA, center',
                    epochs = {'start': startTime, 'stop': endTime, 'step': '1m'})
        eph = obj.ephemerides(quantities='19,20,24', no_fragments=True, closest_apparition=True)

        self.df_eph = eph.to_pandas()
        self.rh = self.df_eph['r'][0].item()
        self.delta = (self.df_eph['delta'][0].item()*u.au).to(u.m)
        self.alpha = self.df_eph['alpha'][0].item()*u.deg  

    def subsolar_temp(self):
        S_sun = const.L_sun / (4 * np.pi * ((1*u.au).to(u.m)**2))
        t_ss_u =  ((1-self.pv*self.phase_int) * S_sun / (self.bolo_emissivity * const.sigma_sb * self.eta * self.rh**2))**(1/4) 
        print(f'Subsolar temperature = {t_ss_u}')
        if self.uncertainties:
            self.t_ss = ((1-self.pv.nominal_value*self.phase_int) * S_sun / (self.bolo_emissivity.nominal_value * const.sigma_sb * self.eta.nominal_value * self.rh**2))**(1/4)
        else:
            self.t_ss = ((1-self.pv*self.phase_int) * S_sun / (self.bolo_emissivity * const.sigma_sb * self.eta * self.rh**2))**(1/4)


    def neatm_integral(self,theta,phi):
        alpha = (self.alpha.to(u.rad)).value
        numerator = (np.cos(phi)**2) * np.cos(theta-alpha)
        if (theta >= alpha - np.pi/2) & (theta <= alpha + np.pi/2):
            temp = self.t_ss * (np.cos(phi)**(1/4)) * (np.cos(theta)**(1/4))
        else:
            temp = 1e-60*u.K
        denominator = np.exp ( (const.h * const.c) / (self.lamda * const.k_B * temp) ) - 1

        return numerator/denominator
    

    def calcDiameter(self,flux):
        #Input flux in Jy - first convert from W/m^2/Hz to W/m^2/m (in place of W/m^2/um so that astropy.units can handle it)
        lamda_flux = flux.to(u.J/u.s/u.m**3, equivalencies=u.spectral_density(self.nu))
        
        #Calculate the integral
        neatm_int = dblquad(self.neatm_integral,0,np.pi/2,-np.pi/2,np.pi/2)
        #Calculate the diameter
        diameter = (lamda_flux * (self.lamda)**5 * (self.delta)**2 / (self.rad_emissivity * const.h * (const.c**2) * neatm_int[0]))**(1/2)
        print(f'Diameter = {diameter.to(u.km)}')

    def tempGas(self,vexp,gamma,mass):
        v = vexp.to(u.m/u.s)
        m = (mass / const.N_A).to(u.kg) 

        self.tgas = v**2 * (gamma - 1) * m / ( (gamma + 1) * gamma * const.k_B.to(u.kg*u.m**2/u.s**2/u.K))
        print(f'T = {self.tgas}')

    def gas_integral(self,theta,phi):
        alpha = (self.alpha.to(u.rad)).value
        numerator = (np.cos(phi)**2) * np.cos(theta-alpha)
        if (theta >= alpha - np.pi/2) & (theta <= alpha + np.pi/2):
            temp = self.tgas * (np.cos(phi)**(1/4)) * (np.cos(theta)**(1/4))
        else:
            temp = 1e-60*u.K
        denominator = np.exp ( (const.h * const.c) / (self.lamda * const.k_B * temp) ) - 1

        return numerator/denominator

    def calcDiameterGas(self,flux,vexp,gamma,mass):

        #Input flux in Jy - first convert from W/m^2/Hz to W/m^2/m (in place of W/m^2/um so that astropy.units can handle it)
        lamda_flux = flux.to(u.J/u.s/u.m**3, equivalencies=u.spectral_density(self.nu))

        #Calculate the subsolar temperature
        self.tempGas(vexp,gamma,mass)
        print(self.tgas)

        #Calculate the integral
        neatm_int = dblquad(self.gas_integral,0,np.pi/2,-np.pi/2,np.pi/2)


        #Calculate the diameter
        diameter = (lamda_flux * (self.lamda)**5 * (self.delta)**2 / (self.rad_emissivity * const.h * (const.c**2) * neatm_int[0]))**(1/2)
        print(f'Diameter = {diameter.to(u.km)}')


    def calcNucleusFlux(self,diameter):
        
        #Calculate the neatm integral
        neatm_int = dblquad(self.neatm_integral,0,np.pi/2,-np.pi/2,np.pi/2)
        flux = ( self.rad_emissivity * (diameter**2) * const.h * (const.c**2) * neatm_int[0]) / ( (self.delta**2) * (self.lamda**5) )
        flux_jy = flux.to(u.mJy,equivalencies=u.spectral_density(self.nu))
        print(f'Flux = {flux_jy:.3f}')

    def calcMassKappa(self,K,kappa,phi):

        
        K_j = K.to(u.J/u.m**2)

        if self.uncertainties:
            temp = 277 * (1 - self.pv.nominal_value * self.phase_int)**(1/4) / (self.rh)**0.5 * u.K
        else:
            temp = 277 * (1 - self.pv * self.phase_int)**(1/4) / (self.rh)**0.5 * u.K

        pbRad = phi.to(u.rad).value

        m = K_j * np.sqrt(np.pi/(4*np.log(2)))*np.pi*pbRad * (self.lamda**2) * (self.delta**2) / (2 * const.k_B * temp * kappa)
        print(f'Mass = {m:.3e}')

    def calcDustFlux(self,mass,kappa):

        if self.uncertainties:
            temp = 277 * (1 - self.pv.nominal_value * self.phase_int)**(1/4) / (self.rh)**0.5 * u.K
        else:
            temp = 277 * (1 - self.pv * self.phase_int)**(1/4) / (self.rh)**0.5 * u.K

        flux = 2 * const.k_B * temp * mass * kappa / (self.lamda**2 * self.delta**2)

        print(f'Flux = {flux.to(u.Jy):.3e}')
        
        return (flux.to(u.Jy))
    
    def calcMassFlux(self,flux,kappa):

        if self.uncertainties:
            temp = 277 * (1 - self.pv.nominal_value * self.phase_int)**(1/4) / (self.rh)**0.5 * u.K
        else:
            temp = 277 * (1 - self.pv * self.phase_int)**(1/4) / (self.rh)**0.5 * u.K

        mass = flux * (self.lamda**2 * self.delta**2) / (2*const.k_B * temp * kappa)
        print(f'Mass = {mass.to(u.kg):.3e}')
        return (mass.to(u.kg))


        