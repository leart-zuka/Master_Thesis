import matplotlib.pyplot as plt
import numpy as np
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import cmath

# -----------------------
# ------ Constants ------ 
# -----------------------

mu0 = -np.sqrt(1/6)  # pi (mf2 -> mf2) #reference transition
mu1 = -np.sqrt(3/10) # pi (mf0 -> mf0)
kappa = 2 * np.pi * 0.063 # cavity dissipation rate
kappa_oc = 2 * np.pi * 0.054 # outcoupling rate
gamma_5P32_5S = 2 * np.pi* 0.006065 / 2  # atom dissipation rate
g0_KC = 2 * np.pi * 0.0438 # coupling strength (F2mf2 -> F'3mf2) # reference transition measured value from M. Brekenfeld thesis
g_pi_KC = g0_KC*(mu1/mu0) # coupling strength (F2mf0 -> F'3mf0) sigma- good
mu_rf = 0.99
mu_fc = 0.9
d_list = np.linspace(-0.5,0.5,1000)

# -----------------------
# ------ Functions ------ 
# -----------------------

def r_empty_cavity(mu_rf:float, mu_fc:float, kappa_m1:float, kappa:float, d_c:float): 
    return mu_rf - mu_fc**2 * (2*kappa_m1/kappa + 1j*d_c)

def r_non_empty_cavity(g:float, gamma:float, d_a: float, d_c: float, kappa: float, kappa_m1: float, mu_rf: float, mu_fc: float):
    C_c = g**2/(2*(gamma+1j*d_a)*(kappa+1j*d_c))
    return mu_rf - mu_fc**2 * (2*kappa_m1)/(kappa+1j*d_c) * 1/(2*C_c + 1)

# ----------------------------
# ------ Atom in Cavity ------
# ----------------------------

d_c_pi = 2 * np.pi * 0
d_a_pi = 2 * np.pi * 0

in_cavity = []
in_cavity_phase = []
for d in d_list:
    val = r_non_empty_cavity(g=g_pi_KC,gamma=gamma_5P32_5S*(1/10+1/10+3/10),d_a=d_a_pi+2*np.pi*d, d_c=d_c_pi+2*np.pi*d,kappa=kappa,kappa_m1=kappa_oc,mu_rf=mu_rf,mu_fc=mu_fc)
    in_cavity.append(np.abs(val))
    in_cavity_phase.append(np.abs(cmath.phase(val)/np.pi))

# --------------------------------
# ------ Atom not in Cavity ------
# --------------------------------

d_c_pi = 2 * np.pi * 0
d_a_pi = 2 * np.pi * 6.3

not_in_cavity = []
not_in_cavity_phase = []
for d in d_list:
    val = r_non_empty_cavity(g=g_pi_KC,gamma=gamma_5P32_5S*(1/2),d_a=d_a_pi+2*np.pi*d,d_c=d_c_pi+2*np.pi*d,kappa=kappa,kappa_m1=kappa_oc,mu_rf=mu_rf,mu_fc=mu_fc)
    not_in_cavity.append(np.abs(val))
    not_in_cavity_phase.append(np.abs(cmath.phase(val)/np.pi))


# ------------------------ 
# ------- Plotting ------- 
# ------------------------ 

image_pi_up = plt.imread('r_pi_up.png')
image_pi_down = plt.imread('r_pi_down.png')

fig, ax = plt.subplots(2, 2,figsize= (18,9))

ax[0,0].plot(d_list,in_cavity, 'r') 
ax[0,0].axvline(x=0.0)
ax[0,0].add_artist(AnnotationBbox(OffsetImage(image_pi_up,zoom=0.34), (0.42,0.35),frameon=True))
ax[0,0].set_xlabel('detuning in (GHz)')
ax[0,0].set_ylabel('reflectivity')

ax[0,1].plot(d_list,not_in_cavity, 'b') 
ax[0,1].axvline(x=0.0)
ax[0,1].add_artist(AnnotationBbox(OffsetImage(image_pi_down,zoom=0.34), (0.42,0.35),frameon=True,zorder=10))
ax[0,1].set_xlabel('detuning in (GHz)')
ax[0,1].set_ylabel('reflectivity')

ax[1,0].plot(d_list, in_cavity_phase, "g")
ax[1,0].axvline(x=0.0)
ax[1,0].set_xlabel("detuning in (GHz)")
ax[1,0].set_ylabel("phase shift in (rad)")

ax[1,1].plot(d_list, not_in_cavity_phase, "k")
ax[1,1].axvline(x=0.0)
ax[1,1].set_xlabel("detuning in (GHz)")
ax[1,1].set_ylabel("phase shift in (rad)")

plt.show()


