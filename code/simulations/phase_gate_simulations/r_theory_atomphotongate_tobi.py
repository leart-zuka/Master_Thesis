# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 15:04:32 2024

@author: tfrank
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.offsetbox import OffsetImage, AnnotationBbox

#%% Parameters


# involved Clebsch-Gorden coefficient D2 line (F2 -> F'3):
mu0 = -np.sqrt(1/6)  # pi (mf2 -> mf2) #reference transition
mu1 = -np.sqrt(3/10) # pi (mf0 -> mf0)
mu2 = np.sqrt(1/5)   # sigma+/- (mf0 -> mf+/-1)

g0_KC = 2 * np.pi * 0.0438 # coupling strength (F2mf2 -> F'3mf2) # reference transition measured value from M. Brekenfeld thesis
g_pi_KC = g0_KC*(mu1/mu0) # coupling strength (F2mf0 -> F'3mf0) sigma- good
g_sp_KC = g0_KC*(mu2/mu0) # coupling strength (F2mf0 -> F'3mf+1) sigma+ good
g_sm_KC = g0_KC*(mu2/mu0) # coupling strength (F2mf0 -> F'3mf-1) sigma- good



G = [["g0_KC/2pi", g0_KC*1000/(2*np.pi), "pi","(F2mf2 -> F'3mf2)"],
     ["g_pi_KC/2pi", g_pi_KC*1000/(2*np.pi), "pi","(F2mf0 -> F'3mf0)"],
     ["g_sp_KC/2pi", g_sp_KC*1000/(2*np.pi), "sigma+","(F2mf0 -> F'3mf+1)"],
     ["g_sm_KC/2pi", g_sm_KC*1000/(2*np.pi), "sigma-","(F2mf0 -> F'3mf-1)"]
     ]


headers = ["coupling","(MHz)", "pol","transition"]
df = pd.DataFrame(G, columns=headers)
print(df)



kappa = 2 * np.pi * 0.063 # cavity dissipation rate
kappa_oc = 2 * np.pi * 0.054 # outcoupling rate
gamma_5P32_5S = 2 * np.pi* 0.006065 / 2  # atom dissipation rate

MM_rf = 0.99
MM_fc = 0.9

# approximation for small detuning?
def r_atom_in_cavity(d_c,d_a,g,kappa, kappa_oc, MM_rf, MM_fc,gamma):
    C_d = g**2/(2*(gamma+1j*d_a)*(kappa+1j*d_c))
    return(MM_rf-(MM_fc**2)*2*kappa_oc/((kappa+1j*d_c)*(2*C_d+1)))

def r_empty_cavity(d_c,kappa, kappa_oc, MM_rf, MM_fc):
    return(MM_rf-(MM_fc**2)*2*kappa_oc/(kappa+1j*d_c))

#%% Atom in F = 2 mf = 0
d_list = np.linspace(-0.5,0.5,1000)

d_c_V = 2 * np.pi * 0.5
d_c_pi = 2 * np.pi * 0
d_a_V = 2 * np.pi * 0
d_a_pi = 2 * np.pi * 0

r1_sp = r_atom_in_cavity(d_c_V,d_a_V,g_sp_KC,kappa, kappa_oc, MM_rf, MM_fc,gamma_5P32_5S*(1/5+4/15+1/30)) # equal to sm
r2 = r_atom_in_cavity(d_c_pi,d_a_pi,g_sm_KC,kappa, kappa_oc, MM_rf, MM_fc,gamma_5P32_5S*(3/10+1/10+1/10))

# Reflection curve
r_up_sp_list = []
r_up_pi_list = []

for d in d_list:
    r_up_sp_list.append(r_atom_in_cavity(d_c_V+2*np.pi*d,d_a_V+2*np.pi*d,g_sp_KC,kappa, kappa_oc, MM_rf, MM_fc,gamma_5P32_5S*(1/3+1/6))) # equal to sm
    r_up_pi_list.append(r_atom_in_cavity(d_c_pi+2*np.pi*d,d_a_pi+2*np.pi*d,g_pi_KC,kappa, kappa_oc, MM_rf, MM_fc,gamma_5P32_5S*(1/10+1/10+3/10)))



#%% Atom in mf = -1
d_c_V = 2 * np.pi * 0.5
d_c_pi = 2 * np.pi * 0
d_a_V = 2 * np.pi * 6.8
d_a_pi = 2 * np.pi * 6.3

r3 = r_atom_in_cavity(d_c_V,d_a_V,g_sp_KC,kappa, kappa_oc, MM_rf, MM_fc,gamma_5P32_5S)
r4 = r_atom_in_cavity(d_c_pi,d_a_pi,g_sm_KC,kappa, kappa_oc, MM_rf, MM_fc,gamma_5P32_5S)

# Reflection curve
r_down_sp_list = [] #atom sits in down state, photon is sigma+ (similar for sigma-)
r_down_pi_list = [] #atom sits in down state, photon is pi

for d in d_list:
    r_down_sp_list.append(r_atom_in_cavity(d_c_V+2*np.pi*d,d_a_V+2*np.pi*d,g_sp_KC,kappa, kappa_oc, MM_rf, MM_fc,gamma_5P32_5S*(1/3+1/6)))
    r_down_pi_list.append(r_atom_in_cavity(d_c_pi+2*np.pi*d,d_a_pi+2*np.pi*d,g_pi_KC,kappa, kappa_oc, MM_rf, MM_fc,gamma_5P32_5S*(1/2)))


# plt.plot(d_list,r3_list)
# plt.plot(d_list,r4_sp_list)
# plt.axvline(x=0.0)

#%% store data
r_list = [r1_sp,r2,r3,r4]
r_data = [["r1_sp",r1_sp],
          ["r2",r2],
          ["r3",r3],
          ["r4",r4]
          ]

headers2 = ["coefficient","theory"]
df2 = pd.DataFrame(r_data, columns=headers2)

print(df2)
#%% Plotting
image_pi_up = plt.imread('r_pi_up.png')
image_sp_up = plt.imread('r_sp_up.png')
image_pi_down = plt.imread('r_pi_down.png')
image_sp_down = plt.imread('r_sp_down.png')
#
fig, ax = plt.subplots(2, 2,figsize= (18,9))

ax[0, 0].plot(d_list,r_up_sp_list, 'r') #row=0, col=0
ax[0, 0].axvline(x=0.0)
ax[0, 0].add_artist(AnnotationBbox(OffsetImage(image_sp_up,zoom=0.34), (0.42,0.35),frameon=True))
ax[0, 0].set_xlabel('detuning in (GHz)')
ax[0, 0].set_ylabel('reflectivity')
# ax[0, 0].text(r_list[0],fontsize=40)

ax[1, 0].plot(d_list,r_up_pi_list, 'b') #row=1, col=0
ax[1, 0].axvline(x=0.0)
ax[1, 0].add_artist(AnnotationBbox(OffsetImage(image_pi_up,zoom=0.34), (0.42,0.35),frameon=True))
ax[1, 0].set_xlabel('detuning in (GHz)')
ax[1, 0].set_ylabel('reflectivity')

ax[0, 1].plot(d_list,r_down_sp_list, 'g') #row=0, col=1
ax[0, 1].axvline(x=0.0)
ax[0, 1].add_artist(AnnotationBbox(OffsetImage(image_sp_down,zoom=0.34), (0.42,0.35),frameon=True))
ax[0, 1].set_xlabel('detuning in (GHz)')
ax[0, 1].set_ylabel('reflectivity')

ax[1, 1].plot(d_list,r_down_pi_list, 'k') #row=1, col=1
ax[1, 1].axvline(x=0.0)
ax[1, 1].add_artist(AnnotationBbox(OffsetImage(image_pi_down,zoom=0.34), (0.42,0.35),frameon=True))
ax[1, 1].set_xlabel('detuning in (GHz)')
ax[1, 1].set_ylabel('reflectivity')

plt.show()

