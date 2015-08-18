#!/usr/bin/env python

import numpy as np
from scipy import integrate
from scipy import interpolate 
from scipy import optimize
import mpmath
from mpmath import mp
import matplotlib.pyplot as plt 

mp.dps = 30

# Physical constants in natural units 
me = 0.51
Alphae = 1/137.0
Ee = mp.sqrt(4 *mp.pi)*Alphae
Bc = me**2/Ee
gauss = Bc/(4.41* 10**13)
fm = 1/197.327
cm = 10**13 * fm
m = 100 *cm
hbarc = 197.327
km = 10**18 *fm
GeV = 10**3
KeV = 10**(-3.0)
rho_0 = 0.16
erg = 624150.9
gram = 5.6*10**(26)
s = (1/6.58 )*10**(22)
minute = 60.0 *s
hour = 60.0 *minute
day = 24*hour
year = 365* day
kg = 10**3 *gram
Msun = 1.9891* 10**30 *kg
mPi = 134.98
fPi = 130.0/np.sqrt(2.0)
z = 0.568
wv = 0.029
M = 938.918;
MN = 939.565
MP = 938.272
Mapprox = 939.0
mP = 1.22*10**(19.0) *GeV
G = mP**(-2.0)
Mz = 91.1876* GeV
Mw = 80.4 *GeV
gw = 0.653
ge = Ee/np.sqrt(4.0* np.pi)
Thetaw = 28.75/180.0
gz = ge/(np.sin(Thetaw)* np.cos(Thetaw))
GF = 1.16636*10**(-11.0)
kelvin = 0.862*10**(-9.0)

#The refitted Skyrme parametrizations
ska25s2009=np.array([-137.191, 153.932, -2181.28, 0.340552, 14605.8, 0.46161, 0.25],float)
ska35s2009=np.array([-172.485, 172.087, -1767.71, 0.282732, 12899.2, 0.413266, 0.35],float)
skt109=np.array([-112.324, 142.467, -1810.72, 0.28223, 12863., 0.392585, 1/3],float)
skt209=np.array([-113.857, 143.999, -1807.87, 0.267778, 12802.4, 0.366144, 1/3],float)
skt309=np.array([-124.432, 148.492, -1812.16, 0.288584, 12906.6, 0.416129, 1/3],float)
eos1=[ska25s2009,ska35s2009,skt109,skt209,skt309]
eosname1=("Ska25s20-0.9","Ska35s20-0.9","SKT1-0.9","SKT2-0.9","SKT3-0.9")
ska25s2010=np.array([7.88383, -0.431231, -2182.96, 0.293491, 14651.7, 0.304538, 0.25],float)
ska35s2010=np.array([-2.41114, -0.507978, -1767.92, 0.247025, 12910.2, 0.220377, 0.35],float)
skt110=np.array([25.175, -12.2601, -1815.65, 0.24924, 12984.8, 0.217553, 1/3],float)
skt210=np.array([-28.3961, 16.7196, -1822.27, 0.258133, 13155.8, 0.24668, 1/3],float)
skt310=np.array([4.51143, 0.0220502, -1816.13, 0.256993, 13025.9, 0.247513, 1/3],float)
eos2=[ska25s2010,ska35s2010,skt110,skt210,skt310]
eosname2=("Ska25s20-1.0","Ska35s20-1.0","SKT1-1.0","SKT2-1.0","SKT3-1.0")
eosname=[eosname1,eosname2]
eos=[eos1,eos2]

# eos1 -> Skyrme Parametirzations with effective neutron mass of 0.9 M
# eos2 -> Skyrme Parametirzations with effective neutron mass of 1.0 M

# Routine to select an interaction for calculations
a=b=t0=x0=t3=x3=alpha=0
def init_skyrme(eff_mass,inter):
	global a,b,t0,x0,t3,x3,alpha 
	a,b,t0,x0,t3,x3,alpha=eos[eff_mass][inter][0],eos[eff_mass][inter][1],eos[eff_mass][inter][2],eos[eff_mass][inter][3]\
							,eos[eff_mass][inter][4],eos[eff_mass][inter][5],eos[eff_mass][inter][6]
	return 


# Effective Masses, Mean Field Shifts based on the parametrizations
def mp(rho,y):
	return M /(1.0 + M /(4.0 * hbarc**2) *rho* (a + 2.0 * y * b))
	
def mn(rho,y):
	return M /(1 + M /(4 * hbarc**2) *rho* (a + 2 * (1-y) * b))

def dmpdrho(rho,y):
	return 0.0-(a+2.0*b*y)/(4.0*hbarc**2)*mp(rho,y)**2

def dmndrho(rho,y):
	return 0.0-(a+2.0*b*(1.0-y))/(4.0*hbarc**2)*mn(rho,y)**2

def dmpdy(rho,y):
	return 0.0-b*rho/(2.0*hbarc**2)*mp(rho,y)**2
	
def dmndy(rho,y):
	return 0.0-b*rho/(2.0*hbarc**2)*mn(rho,y)**2

def up(rho,y):
	return (1/2.0* t0* (2.0 + x0 - y - 2.0* x0* y)* rho - 
  1/24.0* t3* (-4.0 - alpha + 2.0* y* (1.0 + (-1.0 + y) *alpha) + \
     x3* (-1.0 + 2.0* y)* (2.0 + (-1.0 + 2.0* y)* alpha))* rho**(1.0 + alpha))
	
def un(rho,y):
	return (0.5* t0* (1.0+y+x0 *(-1.0+2.0* y))* rho-1/24.0* t3 *(-2.0 *(1+y+x3 *(-1.0+2.0 *y))\
	+(-1.0+x3 *(1.0-2.0 *y)**2.0+2.0 *(-1.0+y) *y)*alpha) *rho**(1.0+alpha))

def dupdrho(rho,y):
	return up(rho,y)/rho-1.0/24.0*t3*alpha*(-4.0-alpha+2.0*y*(1.0-(1.0-y)*alpha)\
	+x3*(2.0*y-1.0)*(2.0+(2.0*y-1.0)*alpha))*rho**alpha

def dundrho(rho,y):
	return un(rho,y)/rho-1.0/24.0*t3*alpha*(-2.0*(1.0+y+x3*(-1.0+2.0*y))\
	+(-1.0+x3*(1.0-2.0*y)**2.0+2.0*(-1.0+y)*y)*alpha)*rho**alpha
	
def dupdy(rho,y):
	return -1.0/12.0*rho*(6.0*t0*(1.0+2.0*x0)+t3*(1.0+2.0*x3)*(1.0+(-1.0+2.0*y)*alpha)*rho**alpha)

def dundy(rho,y):
	return 1.0/12.0*rho*(6.0*t0*(1.0+2.0*x0)-t3*(1.0+2.0*x3)*(-1.0+(-1.0+2.0*y)*alpha)*rho**alpha)

#Integrands for desnity, kinetic energy, entropy calculations

def fermi(k,mass,ushift,Tt,mut):
	arg=(k**2/(2.0*mass)-mut+ushift)/Tt
	f=1.0/(1.0+mpmath.exp(arg))
	return f

def fermi_der(k,mass,ushift,Tt,mut):
	arg=(k**2/(2.0*mass)-mut+ushift)/Tt
	return -mpmath.exp(arg)*fermi(k,mass,ushift,Tt,mut)**2.0

def rho_int(k,mass,ushift,Tt,mut):
	return (k/np.pi)**2*fermi(k,mass,ushift,Tt,mut)*fm**3.0

def drhodrho_int(k,mass,dmdrho,ushift,dudrho,Tt,mut):
	return (k/np.pi)**2*fermi_der(k,mass,ushift,Tt,mut)*\
	(dudrho/T-dmdrho*k**2/2.0/mass**2/T)*fm**3.0

def drhody_int(k,mass,dmdy,ushift,dudy,Tt,mut):
	return (k/np.pi)**2*fermi_der(k,mass,ushift,Tt,mut)*\
	(dudy/T-dmdy*k**2/2.0/mass**2/T)*fm**3.0
	
def drhodrho(mass,dmdrho,ushift,dudrho,Tt,mut):
	return integrate.quad(drhodrho_int,0,np.inf,args=(mass,dmdrho,ushift,dudrho,Tt,mut),\
	epsabs=1.49e-06, epsrel=1.49e-06,limit=50)[0]

def drhody(mass,dmdrho,ushift,dudrho,Tt,mut):
	return integrate.quad(drhody_int,0,np.inf,args=(mass,dmdrho,ushift,dudrho,Tt,mut),\
	epsabs=1.49e-06, epsrel=1.49e-06,limit=50)[0]
	
def rho_f(mass,ushift,T,mu):
	return integrate.quad(rho_int,0,np.inf,args=(mass,ushift,T,mu),\
	epsabs=1.49e-06, epsrel=1.49e-06,limit=50)[0]

def rho_n(mun,rho,y,T):
	return rho_f(mn(rho,y),un(rho,y),T,mun)

def rho_p(mup,rho,y,T):
	return rho_f(mn(rho,y),up(rho,y),T,mup)

def tau_int(k,mass,ushift,T,mu):
	return k**4/np.pi**2*fermi(k,mass,ushift,T,mu)*fm**5.0

def tau_f(mass,ushift,T,mu):
	return integrate.quad(tau_int,0,np.inf,args=(mass,ushift,T,mu),\
	epsabs=1.49e-06, epsrel=1.49e-06,limit=50)[0]

def tau_n_rho(rho,y,mun,T):
	return tau_f(mn(rho,y),un(rho,y),T,mun)

def tau_p_rho(rho,y,mup,T):	
	return tau_f(mp(rho,y),up(rho,y),T,mup)
	
def s_int(k,mass,ushift,T,mu):
	f=fermi(k,mass,ushift,T,mu)
	integrand=-k**2/np.pi**2*(mpmath.log(f**f)+mpmath.log((1.0-f)**(1.0-f)))*fm**3.0
	return integrand
	
def s_f(mass,ushift,T,mu):
	sf=integrate.quad(s_int,0.0,np.inf,args=(mass,ushift,T,mu),\
	epsabs=1.49e-06, epsrel=1.49e-06,limit=50)[0]
	return sf

# This section treats rho, y as independent parameters and gives 
# the rest of the thermodynamic potentials as functions of (rho,y,T)

# Routines for the user:
# mu_N(rho,y,T), mu_N(rho,y,T):
#		the chemical potentials from total baryon density, 
#          proton fraction, and temperature
# Ek_rho(rho,y,T),Ei_rho(rho,y,T), Et_rho(rho,y,T):
#		the energy per baryon -> kinetic, interaction, total
# s_N_rho(rho,y,T), s_P_rho(rho,y,T), St_rho(rho,y,T):
#		the entropy per baryon -> neutron, proton, total
# Ft_rho(rho,y,T),Pt_rho(rho,y,T) -> free energy per bryon , pressure	
# skyrme_data_rho(rho,y,T) -> all the data:
# [rho, yp, T, mun, mup, Et, Ft, St, Pt]
MU_MIN = MU_MAX=0
MU_MIN_C=-1000.0
MU_MAX_C = 1000.0
RHO_MIN=0

def mu_nt0(rho,y,T):
	return un(rho,y)+(3.0**(2.0/3.0)*np.pi**(4.0/3.0)*rho**(2.0/3.0)*fm**(-2.0))/(2.0*mn(rho,y))

def mu_pt0(rho,y,T):
	return up(rho,y)+(3.0**(2.0/3.0)*np.pi**(4.0/3.0)*rho**(2.0/3.0)*fm**(-2.0))/(2.0*mp(rho,y))


def mu_eq(mu,rho_type,rho,y,T):
	if rho_type == rho_n:
		rho_target = (1.0-y)*rho
	else:
		rho_target = y*rho
	return rho_target-rho_type(mu,rho,y,T)

def mu_rho(rho_type,rho,y,T):
	if rho_type == rho_n:
		MU_MIN = MU_MIN_C*mu_nt0(rho,y,T)
		MU_MAX = MU_MAX_C*mu_nt0(rho,y,T)
	else:
		MU_MIN = MU_MIN_C*mu_pt0(rho,y,T)
		MU_MAX = MU_MAX_C*mu_pt0(rho,y,T)
	mumin,mumax=MU_MIN,MU_MAX
	dmu=0.001
	while(mu_eq(mumin,rho_type,rho,y,T)*mu_eq(mumax,rho_type,rho,y,T)>0 and mumax>mumin):
		mumax=mumax-dmu
	return optimize.brentq(mu_eq,mumin,mumax,args=(rho_type,rho,y,T))


def mu_N(rho,y,T):
	return mu_rho(rho_n,rho,y,T)

def mu_N_list(rho_array,y,T):
	mu_array = np.zeros(len(rho_array))
	for i in range(0,len(rho_array)):
		mu_array[i]=mu_N(rho_array[i],y,T)
	return mu_array
	
def mu_P(rho,y,T):
	return mu_rho(rho_p,rho,y,T)

def mu_P_list(rho_array,y,T):
	mu_array = np.zeros(len(rho_array))
	for i in range(0,len(rho_array)):
		mu_array[i]=mu_P(rho_array[i],y,T)
	return mu_array
	
def rho_N_y(rho,y,T):
	return rho_n(mu_N(rho,y,T),rho,y,T)

def rho_P_y(rho,y,T):
	return rho_p(mu_P(rho,y,T),rho,y,T)
	
def tau_N_rho_y(rho,y,T):
	return tau_f(mn(rho,y),un(rho,y),T,mu_N(rho,y,T))

def tau_P_rho_y(rho,y,T):	
	return tau_f(mp(rho,y),up(rho,y),T,mu_P(rho,y,T))
	
def Ek_rho_y(rho,y,T):
	return (tau_N_rho_y(rho,y,T)+tau_P_rho_y(rho,y,T))*hbarc**2/(2.0*M*rho)
		
def Ei_rho_y(rho,y,T):
	taun=tau_N_rho(rho,y,T)
	taup=tau_P_rho(rho,y,T)
	return (-0.25 *t0 *(-1.0+x0 *(1.0-2.0 *y)**2 + 2.0*(-1.0 + y)*y)*\
	 rho-1.0/24.0*t3*(-1.0+x3*(1.0-2.0*y)**2 \
     +2.0*(-1.0 + y)*y)* rho**(1.0 +  alpha) + \
     1/8.0*rho*(a*(taun+taup) + \
     2.0*b*((1.0-y)*taun+y*taup)))

def Et_rho_y(rho,y,T):
	taun=tau_N_rho(rho,y,T)
	taup=tau_P_rho(rho,y,T)
	ek=(taun+taup)*hbarc**2/(2.0*M*rho)
	ei=(-0.25 *t0 *(-1.0+x0 *(1.0-2.0 *y)**2 + 2.0*(-1.0 + y)*y)*\
	 rho-1.0/24.0*t3*(-1.0+x3*(1.0-2.0*y)**2 \
     +2.0*(-1.0 + y)*y)* rho**(1.0 +  alpha) + \
     1/8.0*rho*(a*(taun+taup) + \
     2.0*b*((1.0-y)*taun+y*taup)))
	return ek+ei	

def s_N_rho_y(rho,y,T):
	return s_f(mn(rho,y),un(rho,y),T,mu_N(rho,y,T))/rho

def s_P_rho_y(rho,y,T):
	return s_f(mp(rho,y),up(rho,y),T,mu_P(rho,y,T))/rho	

def St_rho_y(rho,y,T):
	return (s_f(mn(rho,y),un(rho,y),T,mu_N(rho,y,T))\
	+s_f(mp(rho,y),up(rho,y),T,mu_P(rho,y,T)))/rho

def St_rho_list(rho_array,y,T):
	s_array = np.zeros(len(rho_array))
	for i in range(0,len(rho_array)):
		s_array[i]=St_rho_y(rho_array[i],y,T)
	return s_array

def Ft_rho_y(rho,y,T):
	taun=tau_N_rho(rho,y,T)
	taup=tau_P_rho(rho,y,T)
	ek=(taun+taup)*hbarc**2/(2.0*M*rho)
	ei=(-0.25 *t0 *(-1.0+x0 *(1.0-2.0 *y)**2 + 2.0*(-1.0 + y)*y)*\
	rho-1.0/24.0*t3*(-1.0+x3*(1.0-2.0*y)**2 \
	+2.0*(-1.0 + y)*y)* rho**(1.0 +  alpha) + \
	1/8.0*rho*(a*(taun+taup) + \
	2.0*b*((1.0-y)*taun+y*taup)))
	st=(s_f(mn(rho,y),un(rho,y),T,mu_N(rho,y,T))\
	+s_f(mp(rho,y),up(rho,y),T,mu_P(rho,y,T)))/rho
	return (ek_ei-T*st)

def Pt_rho_y(rho,y,T):
	MN,UN,MUN=mn(rho,y),un(rho,y),mu_N(rho,y,T)
	MP,UP,MUP=mp(rho,y),up(rho,y),mu_P(rho,y,T)
	taun=tau_f(MN,UN,T,MUN)
	taup=tau_f(MP,UP,T,MUP)
	sn=s_f(MN,UN,T,MUN)
	sp=s_f(MP,UP,T,MUP)
	ek=(taun+taup)*hbarc**2/(2.0*M)
	ei=rho*(-0.25 *t0 *(-1.0+x0 *(1.0-2.0 *y)**2 + 2.0*(-1.0 + y)*y)*\
	 rho-1.0/24.0*t3*(-1.0+x3*(1.0-2.0*y)**2 \
     +2.0*(-1.0 + y)*y)* rho**(1.0 +  alpha) + \
     1/8.0*rho*(a*(taun+taup)) + \
     2.0*b*((1.0-y)*taun+y*taup))
	return rho*(MUN*(1.0-y)+y*MUP)+T*(sn+sp)-ek-ei

def Pt_rho_mun_mup(rho,y,MUN,MUP,T):
	MN,UN=mn(rho,y),un(rho,y)
	MP,UP=mp(rho,y),up(rho,y)
	taun=tau_f(MN,UN,T,MUN)
	taup=tau_f(MP,UP,T,MUP)
	sn=s_f(MN,UN,T,MUN)
	sp=s_f(MP,UP,T,MUP)
	ek=(taun+taup)*hbarc**2/(2.0*M)
	ei=rho*(-0.25 *t0 *(-1.0+x0 *(1.0-2.0 *y)**2 + 2.0*(-1.0 + y)*y)*\
	 rho-1.0/24.0*t3*(-1.0+x3*(1.0-2.0*y)**2 \
     +2.0*(-1.0 + y)*y)* rho**(1.0 +  alpha) + \
     1/8.0*rho*(a*(taun+taup)) + \
     2.0*b*((1.0-y)*taun+y*taup))
	return rho*(MUN*(1.0-y)+y*MUP)+T*(sn+sp)-ek-ei

def rhoy_ptmu_eq(p,Pt,mun,mup,T):
	rho,y= p
	return (100.0*abs((Pt-Pt_rho_mun_mup(rho,y,mun,mup,T))/(Pt+10.0**(-6.0))),(abs(mu_eq(mun,rho_n,rho,y,T))/(rho*(1.0-y)+10**(-6.0))\
	+abs(mu_eq(mup,rho_p,rho,y,T))/(y*rho+10.0**(-6.0)))*100.0)
		
def rho_y_PTMU(rho,y,Pt,mun,mup,T):
	rho_sol,y_sol=optimize.fsolve(rhoy_ptmu_eq,(rho,y),args=(Pt,mun,mup,T))
	return (rho_sol,y_sol)

def Pt_rho_y_list(rho_array,y,T):
	p_array = np.zeros(len(rho_array))
	for i in range(0,len(rho_array)):
		p_array[i]=Pt_rho_y(rho_array[i],y,T)
	return p_array
	

def skyrme_data_rho_y(rho,y,T):
	return np.array([rho,y,T,mu_N(rho,y,T),mu_P(rho,y,T),Et_rho_y(rho,y,T),\
	Ft_rho_y(rho,y,T),St_rho_y(rho,y,T),Pt_rho_y(rho,y,T)])

def skyrme_data_table_rho_y(rho_array,y_array,T_array):
	 data=np.zeros( (len(rho_array),len(y_array),len(T_array)) )
	 for i in len(rho_array):
		 for j in len(y_array):
			 for k in len(T_array):
				 data[i][j][k]=skyrme_data_rho_y(rho_array[i],y_array[j],T_array[k])
	 return data
	

eff_mass=1
inter=3
init_skyrme(eff_mass,inter)
YPT= 0.3
T=2.0
rho_b=0.005
pt=Pt_rho_y(rho_b,YPT,T)
MUN=mu_N(rho_b,YPT,T)
MUP=mu_P(rho_b,YPT,T)
#rhoy=rho_y_PTMU(0.2,0.5,pt,MUN,MUP,T)
print("Parametrization chosen: ",eosname[eff_mass][inter])
print("The input parameters for the gas:")
print("rho = ",rho_b," fm^-3, yp = ",YPT,", T = ",T," MeV")
print("Output values from input:")
print("mu_n = ",MUN," MeV, mu_p = ",MUP," MeV, P = ",pt," MeV fm^-3")
#print("Using P,mu_n,mu_p find (rho,y) for the liquid by solving the system of equations:")
#print( "(ABS( (P-P(rho,y))/P )*100 , ABS( (rho_n-rho_n(rho,y,mu_n))/rho_n ) + ABS( (rho_p-rho_p(rho,y,mu_n))/rho_p )*100")
#print("where, rho_n = (1-y) * rho, rho_p = y * rho")
#print("Initial guess for (rho,y): (", 0.2," fm^-3, ",0.5,")")
#print("The result of the equations: ",rhoy_ptmu_eq(rhoy,pt,MUN,MUP,T))
#pt2=Pt_rho_mun_mup(rhoy[0],rhoy[1],MUN,MUP,T)
#print("Parameters for liquid phase: (rho, y, mu_n, mu_p, P)")
#print((rhoy[0],rhoy[1],MUN,MUP,pt2))
