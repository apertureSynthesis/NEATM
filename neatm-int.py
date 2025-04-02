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

# nu1 = 216*u.GHz
# nu2 = 232*u.GHz
nu3 = 233*u.GHz
# lamda1 = (nu1.to(u.m,equivalencies=u.spectral()))
# lamda2 = (nu2.to(u.m,equivalencies=u.spectral()))
lamda3 = (nu3.to(u.m,equivalencies=u.spectral()))

#March 8
r1 = 20.0046*u.au
d1 = 19.6776*u.au
delta1 = d1.to(u.m)
phi1 = 2.77*u.deg
alpha1 = phi1.to(u.rad)
print('delta1 = {:.3e}'.format(delta1))

# #March 17
# r2 = 16.58*u.au
# d2 = 16.88*u.au
# delta2 = d2.to(u.m)
# phi2 = 3.24*u.deg
# alpha2 = phi2.to(u.rad)


#Calculate T_ss
S_sun = const.L_sun / (4*np.pi*((1*u.au).to(u.m)**2))

#Assuming NEATM formalism
T_ss1 = ( (1-A)*S_sun / (emis*const.sigma_sb*eta*(r1.value**2)) )**(1/4)
#T_ss2 = ( (1-A)*S_sun / (emis*const.sigma_sb*eta*(r2.value**2)) )**(1/4)

pv_u = ufloat(pv,dpv)
emiss_u = ufloat(emis,demis)
eta_u = ufloat(eta,deta)
diam_m_u = ufloat(diam_m.value,d_diam_m.value)
emiss_rad_u = ufloat(emis_rad,demis_rad)

T_ss1_u = ( (1 - (pv_u*phase_integral)) * S_sun.value / (emiss_u*const.sigma_sb.value*eta_u*(r1.value**2)))**(1/4)
#T_ss2_u = ( (1 - (pv_u*phase_integral)) * S_sun.value / (emiss_u*const.sigma_sb.value*eta_u*(r2.value**2)))**(1/4)

print('T_ss1 = {}'.format(T_ss1))
#print('T_ss2 = {}'.format(T_ss2))

print('T_ss1_u = {}'.format(T_ss1_u))
#print('T_ss2_u = {}'.format(T_ss2_u))

#Using vexp measured from CO
# v1_u = ufloat(389,7)
# v2_u = ufloat(384,14)

# v1 = 386
# v2 = 394

# gamma = 1.4
# R = const.R.value
# m_co = 0.028010 #kg/mol

# T_ss1 = ( v1**2 * (gamma-1)/2 * m_co / (gamma*R)) * u.K
# T_ss2 = ( v2**2 * (gamma-1)/2 * m_co / (gamma*R)) * u.K

# T_ss1_u = ( v1_u**2 * (gamma-1)/2 * m_co / (gamma*R))
# T_ss2_u = ( v2_u**2 * (gamma-1)/2 * m_co / (gamma*R))

# print('T_ss1 = {}'.format(T_ss1_u))
# print('T_ss2 = {}'.format(T_ss2_u))

#NEATM integral numerator

def neatm_numerator(theta,phi,alpha):
    return (np.cos(phi)**2) * np.cos(theta-alpha)

def neatm_denominator(theta,phi,t_ss,lamda):
    temp = t_ss * (np.cos(phi))**(1/4) * (np.cos(theta))**(1/4)
    return np.exp( (const.h*const.c) / (lamda * const.k_B * temp) ) - 1

def neatm_int(theta,phi,alpha,t_ss,lamda):

    numerator = (np.cos(phi)**2) * np.cos(theta-alpha.value)
    if (theta >= alpha.value - np.pi/2) & (theta <= alpha.value + np.pi/2):
        temp = t_ss * (np.cos(phi)**(1/4)) * (np.cos(theta)**(1/4))
    else:
        temp = 0*u.K
    denominator = np.exp( (const.h*const.c) / (lamda * const.k_B * temp) ) - 1

    return numerator / denominator




# neatm_int1_lamda1 = dblquad(neatm_int,0,np.pi/2,-np.pi/2,np.pi/2,args=(alpha1,T_ss1,lamda1))
# neatm_int1_lamda2 = dblquad(neatm_int,0,np.pi/2,-np.pi/2,np.pi/2,args=(alpha1,T_ss1,lamda2))
neatm_int1_lamda3 = dblquad(neatm_int,0,np.pi/2,-np.pi/2,np.pi/2,args=(alpha1,T_ss1,lamda3))

# neatm_int2_lamda1 = dblquad(neatm_int,0,np.pi/2,-np.pi/2,np.pi/2,args=(alpha2,T_ss2,lamda1))
# neatm_int2_lamda2 = dblquad(neatm_int,0,np.pi/2,-np.pi/2,np.pi/2,args=(alpha2,T_ss2,lamda2))
#neatm_int2_lamda3 = dblquad(neatm_int,0,np.pi/2,-np.pi/2,np.pi/2,args=(alpha2,T_ss2,lamda3))



#print(neatm_int1_lamda1[0])
#print(emis * (diam_m**2) * const.h * (const.c**2) / ((delta1**2) * (lamda1**5)))


def flux(emis,diam,delta,lamda,neatm):

    return ( emis * (diam**2) * const.h * (const.c**2) * neatm ) / ( (delta**2) * (lamda**5) )

# flux1_lamda1 = flux(emiss_rad_u, diam_m_u, delta1, lamda1, neatm_int1_lamda1[0])
# flux1_lamda2 = flux(emiss_rad_u, diam_m_u, delta1, lamda2, neatm_int1_lamda2[0])
flux1_lamda3 = flux(emiss_rad_u, diam_m_u, delta1, lamda3, neatm_int1_lamda3[0])

# print('March 8 flux at 216 GHz = {:.3f} mJy'.format(1e29 * flux1_lamda1.value * (lamda1.value) / (nu1.to(u.Hz).value)))
# print('March 8 flux at 232 GHz = {:.3f} mJy'.format(1e29 * flux1_lamda2.value * (lamda2.value) / (nu2.to(u.Hz).value)))
print('March 8 flux at 233 GHz = {:.3f} mJy'.format(1e29 * flux1_lamda3.value * (lamda3.value) / (nu3.to(u.Hz).value)))

# flux2_lamda1 = flux(emiss_rad_u, diam_m_u, delta2, lamda1, neatm_int2_lamda1[0])
# flux2_lamda2 = flux(emiss_rad_u, diam_m_u, delta2, lamda2, neatm_int2_lamda2[0])
# flux2_lamda3 = flux(emiss_rad_u, diam_m_u, delta2, lamda3, neatm_int2_lamda3[0])

# print('Mar 17 flux at 216 GHz = {:.3f} mJy'.format(1e29 * flux2_lamda1.value * (lamda1.value) / (nu1.to(u.Hz).value)))
# print('Mar 17 flux at 232 GHz = {:.3f} mJy'.format(1e29 * flux2_lamda2.value * (lamda2.value) / (nu2.to(u.Hz).value)))
# print('Mar 17 flux at 233 GHz = {:.3f} mJy'.format(1e29 * flux2_lamda3.value * (lamda3.value) / (nu3.to(u.Hz).value)))
