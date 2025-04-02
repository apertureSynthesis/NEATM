import astropy.units as u
import astropy.constants as const
import numpy as np
from scipy.integrate import dblquad
from uncertainties import ufloat
from uncertainties.umath import *

phase_integral = 0.4
emis_rad = 0.70
demis_rad = 0.13
emis = 0.90
demis = 0.06
eta = 1.175
deta = 0.42
pv = 0.053
dpv = 0.012
A = phase_integral * pv
diam = 137*u.km
d_diam = 17*u.km
diam_m = diam.to(u.m)
d_diam_m = d_diam.to(u.m)

#Putting Hz as 1/s for astropy
# nu1 = 216e9*1/u.s
# nu2 = 232e9*1/u.s
nu3 = 233e9*1/u.s

#Lellouch
r1 = 20.0046*u.au
d1 = 19.6776*u.au
delta1 = d1.to(u.m)
phi1 = 2.77*u.deg
alpha1 = phi1.to(u.rad)
print('delta1 = {:.3e}'.format(delta1))

#Calculate T_ss
S_sun = const.L_sun / (4*np.pi*((1*u.au).to(u.m)**2))

T_ss1 = ( (1-A)*S_sun / (emis*const.sigma_sb*eta*(r1.value**2)) )**(1/4)

print('T_ss1 = {}'.format(T_ss1))

pv_u = ufloat(pv,dpv)
emiss_u = ufloat(emis,demis)
eta_u = ufloat(eta,deta)
diam_m_u = ufloat(diam_m.value,d_diam_m.value)
emiss_rad_u = ufloat(emis_rad,demis_rad)

T_ss1_u = ( (1 - (pv_u*phase_integral)) * S_sun.value / (emiss_u*const.sigma_sb.value*eta_u*(r1.value**2)))**(1/4)

print('T_ss1_u = {}'.format(T_ss1_u))

def neatm_int(theta,phi,alpha,t_ss,nu):

    numerator = (np.cos(phi)**2) * np.cos(theta-alpha.value)
    if (theta >= alpha.value - np.pi/2) & (theta <= alpha.value + np.pi/2):
        temp = t_ss * (np.cos(phi)**(1/4)) * (np.cos(theta)**(1/4))
    else:
        temp = 0*u.K
    denominator = np.exp( (const.h*nu) / (const.k_B * temp) ) - 1

    return numerator / denominator

def flux(emis,diam,delta,nu,neatm):

    return ( emis * (diam**2) * const.h * (nu**3) * neatm ) / ( (delta**2) * (const.c**2) )

neatm_int1_nu1 = dblquad(neatm_int,0,np.pi/2,-np.pi/2,np.pi/2,args=(alpha1,T_ss1,nu3))



flux1_nu1 = flux(emiss_rad_u, diam_m_u, delta1, nu3, neatm_int1_nu1[0])


print('March 8 flux at 216 GHz = {:.3f} mJy'.format(1e29*flux1_nu1.value))
