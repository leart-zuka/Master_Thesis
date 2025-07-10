from typing import Dict
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import cmath

# -----------------------
# ------ Constants ------ 
# -----------------------

# Clebsch-Gorden coefficient D2 line (F2 -> F'3)
Mu0 = np.sqrt(1/6)  # pi (mf2 -> mf2)
Mu1 = np.sqrt(1/3)  # simga+ (mf1 -> mf2)
Mu2 = np.sqrt(1/10) # simga- (mf1 -> mf0)
Mu3 = np.sqrt(1/10) # sigma+ (mf-1 -> mf0)

G0_kc = 2 * np.pi * 0.0438 # coupling strength (F2mf2 -> F'3mf2)
G1_kc = G0_kc*(Mu1/Mu0) # coupling strength (F1mf1 -> F'3mf2)
G2_kc = G0_kc*(Mu2/Mu0) # coupling strength (F1mf1 -> F'3mf0)
G3_kc = G0_kc*(Mu3/Mu0) # coupling strength (F1mf-1 -> F'3mf0)

Kappa = 2 * np.pi * 0.063 # cavity dissipation rate
Kappa_oc = 2 * np.pi * 0.054
Gamma_5P32_5S = 2 * np.pi* 0.006065 / 2  # atom dissipation rate

Mu_rf = 1
Mu_fc = 0.9

Atom_dimensions = 3 # |0>, |1>, |e>
Photon_dimensions = [2] # only π-pol. light is able to enter our cavity

def r_non_empty_cavity(g:float, gamma:float, d_a: float, d_c: float, kappa: float, kappa_m1: float, mu_rf: float, mu_fc: float):
    C_c = g**2/(2*(gamma+1j*d_a)*(kappa+1j*d_c))
    return mu_rf - mu_fc**2 * (2*kappa_m1)/(kappa+1j*d_c) * 1/(2*C_c + 1)

def input_shape(t:float, args: Dict[str,float]) -> float:
    t0 = args['t0']
    return np.exp(-(t-t0/2)**2/(t0/5)**2)
