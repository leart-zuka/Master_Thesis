# -*- coding: utf-8 -*-
"""
Created on Wed May 15 18:12:56 2024

@author: tfrank
"""




from qutip import *
import pandas as pd
from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np


from CavityOperators import operators

###################
# Model parameters
###################

N = 2  # number of cavity fock states (0,1)
d = 0 # atom-cavity detuning


# Clebsch-Gorden coefficient D2 line (F2 -> F'3)
mu0 = np.sqrt(1/6)  # pi (mf2 -> mf2)
mu1 = np.sqrt(1/3)  # simga+ (mf1 -> mf2)
mu2 = np.sqrt(1/10) # simga- (mf1 -> mf0)
mu3 = np.sqrt(1/10) # sigma+ (mf-1 -> mf0)

g0_KC = 2 * np.pi * 0.0438 # coupling strength (F2mf2 -> F'3mf2)
g1_KC = g0_KC*(mu1/mu0) # coupling strength (F1mf1 -> F'3mf2)
g2_KC = g0_KC*(mu2/mu0) # coupling strength (F1mf1 -> F'3mf0)
g3_KC = g0_KC*(mu3/mu0) # coupling strength (F1mf-1 -> F'3mf0)

G = [["g0_KC/2pi", g0_KC*1000/(2*np.pi), "pi","(F2mf2 -> F'3mf2)"],
        ["g1_KC/2pi", g1_KC*1000/(2*np.pi), "sigma+","(F1mf1 -> F'3mf2)"],
        ["g2_KC/2pi", g2_KC*1000/(2*np.pi), "sigma-","(F1mf1 -> F'3mf0)"],
        ["g3_KC/2pi", g3_KC*1000/(2*np.pi), "sigma+","(F1mf-1 -> F'3mf0)"]]

headers = ["coupling","(MHz)", "pol","transition"]

df = pd.DataFrame(G, columns=headers)
print(df)

kappa = 2 * np.pi * 0.063 # cavity dissipation rate
kappa_oc = 2 * np.pi * 0.054
gamma_5P32_5S = 2 * np.pi* 0.006065 / 2  # atom dissipation rate

MM_rf = 1
MM_fc = 0.9

atomdim = [5] #5 atomic levels
photondim = [2, 2] # 2modes: a[0] sigma +/ sigma- a[1], photons [0,1]

#%%


tlist = np.linspace(0, 2000, 100)


#Setup the operators, the Hamiltonian and initial state

# creating operators
a, S, astate = operators(photondim, atomdim)

# JC Hamiltonian in dipole interaction form and rotating wave approximation
args = {'sigma': 100, 't0': 1000}

def Field_in(t,args):
    sigma = args['sigma']
    t0 = args['t0']
    return(np.sqrt(1/(sigma * np.sqrt(2*np.pi)) * np.exp(-(t-t0)**2 / sigma**2 / 2)))
           
def r_theory(d,g,kappa, kappa_oc, MM_rf, MM_fc,gamma):
    C_d = g**2/(2*(gamma+1j*d)*(kappa+1j*d))
    return(MM_rf-(MM_fc**2)*2*kappa_oc/((kappa+1j*d)*(2*C_d+1)))



#Create a list of collapse operators that describe the dissipation
c_ops = [np.sqrt(2 * kappa) * a[0], 
         np.sqrt(2 * kappa) * a[1], #outcoupling +/- mode     
         np.sqrt(2 * gamma_5P32_5S * 1/6) * S[0][4][1],
         np.sqrt(2 * gamma_5P32_5S * 1/3) * S[0][0][1],
         np.sqrt(2 * gamma_5P32_5S * 1/10) * S[0][0][2],
         np.sqrt(2 * gamma_5P32_5S * 1/10) * S[0][3][2],
         np.sqrt(2 * gamma_5P32_5S * 3/10) * S[0][4][2],
         ]


# Observables
obs = [
       S[0][0][0], 
       S[0][1][1], 
       S[0][2][2],     
       S[0][3][3], 
       S[0][4][4],
       a[0].dag() * a[0], 
       a[1].dag() * a[1], 
       a[0], 
       a[1]
       ]

# Initialm state
psi0 = tensor(astate[0][0],
              fock(photondim[0], 0),
              fock(photondim[0], 0))

r_sigmaM =[]
r_sigmaP =[]
r_th_sigmaP = []
r_th_sigmaM = []
R_sigmaM =[]
R_sigmaP =[]
R_th_sigmaP = []
R_th_sigmaM = []

D = np.linspace(-.2,.2,160)
for d in D:
    H_D = 2* np.pi * d * (a[0].dag()*a[0] - S[0][1][1]) + 2* np.pi * d * (a[1].dag()*a[1] - S[0][0][0]) + 2* np.pi * d * (a[0].dag()*a[0] - S[0][3][3])
    H_g = g1_KC * S[0][1][0] * a[0] + g2_KC * S[0][2][0] * a[1] + g3_KC * S[0][3][0] * a[0]
    H_drive = 0.8*np.sqrt(2 * kappa_oc)*(a[0] + a[0].dag()+a[1] + a[1].dag())
    H = [H_D + H_g + H_g.dag(),[H_drive, Field_in]]

    #Evolve the system
    output = mesolve(H, psi0, tlist, c_ops, obs, args)
    R_sigmaP.append(sum((np.nan_to_num(abs(Field_in(tlist,args)-0.8*1j*np.sqrt(2*kappa_oc)*output.expect[7]))**2))*(tlist[1]-tlist[0]))
    r_th_sigmaP.append(r_theory(d*2* np.pi,g1_KC,kappa, kappa_oc, MM_rf, MM_fc,gamma_5P32_5S/2))
    R_th_sigmaP.append(abs(r_theory(d*2* np.pi,g1_KC,kappa, kappa_oc, MM_rf, MM_fc,gamma_5P32_5S/2))**2)
    
    R_sigmaM.append(sum((np.nan_to_num(abs(Field_in(tlist,args)-0.8*1j*np.sqrt(2*kappa_oc)*output.expect[8]))**2))*(tlist[1]-tlist[0]))
    r_th_sigmaM.append(r_theory(d*2* np.pi,g2_KC,kappa, kappa_oc, MM_rf, MM_fc,gamma_5P32_5S/2))
    R_th_sigmaM.append(abs(r_theory(d*2* np.pi,g2_KC,kappa, kappa_oc, MM_rf, MM_fc,gamma_5P32_5S/2))**2)
    print(R_sigmaP[-1])



d = 0.5
H_D = 2* np.pi * d * (a[0].dag()*a[0] - S[0][0][0])
H_g = g1_KC * S[0][1][0] * a[0] #+ g2_KC * S[0][2][0] * a[1] + g3_KC * S[0][3][2] * a[0]
H_drive = 0.8*np.sqrt(kappa_oc)*(a[0] + a[0].dag())
H = [H_D + H_g + H_g.dag(),[H_drive, Field_in]]


#Evolve the system
output = mesolve(H, psi0, tlist, c_ops, obs, args)
#%%
fig1,ax = plt.subplots(2,2,figsize=(16,5))
ax[0][0].plot(tlist,Field_in(tlist,args)**2, label="field_in", color = "r")
ax[0][0].plot(tlist,output.expect[2], label="outcoupled field")
ax[0][0].plot(tlist,(Field_in(tlist,args)-0.01*np.sqrt(2*kappa_oc)*output.expect[4])**2, label="R", color = 'g')
ax[0][0].legend()

R = np.nan_to_num(Field_in(tlist,args)-0.8*np.sqrt(2*kappa_oc)*output.expect[2])**2
print(sum(R), sum(abs(Field_in(tlist,args))**2))
# ax[0][1].plot(tlist,Field_in(tlist,args), label="field_in", color = "r")
# ax[0][1].plot(tlist,Field_in(tlist,args)-output.expect[2], label="r", color = "g")
# ax[0][1].legend()
# ax[1][0].plot(tlist,abs(Field_in(tlist,args)**2), label="drive", color = "r")
# ax[1][0].plot(tlist,abs(output.expect[2])**2, label="P_cavity")
# ax[1][0].plot(tlist,output.expect[0], label="P_cavity")
# ax[1][0].legend()
# ax[1][1].plot(tlist,abs(Field_in(tlist,args)**2), label="drive", color = "r")
# ax[1][1].plot(tlist,abs(Field_in(tlist,args)-output.expect[2])**2, label="R", color = "g")
# ax[1][1].legend()
#ax[0].tight_layout()
#Visualize the results

fig, ax = plt.subplots(figsize=(8, 5))
# ax.plot(D, r, label="simulation_steady")
ax.plot(D, R_th_sigmaP, label="theory+")
ax.plot(D, R_sigmaP, label="simulation_pulse+")
ax.legend()

fig2, ax = plt.subplots(figsize=(8, 5))
# ax.plot(D, r, label="simulation_steady")
ax.plot(D, R_th_sigmaM, label="theory-")
ax.plot(D, R_sigmaM, label="simulation_pulse-")
ax.legend()

