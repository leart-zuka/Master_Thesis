# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 17:02:47 2024

@author: tfrank
"""

'''reflection coefficients'''


import matplotlib.pyplot as plt
import numpy as np

R_V = (21.8-2.7)/(97.2-2.7)
print(R_V)
print(np.sqrt(R_V))

#memory based photon-photon gate
r1_sm,r1_sp,r2,r3,r4=[(0.8259169681505598+0j),
                    (0.9353272870477435+0j),
                    (0.9240692099477851+0.08044411008325653j),
                    (-0.16258865581210658-0.01032288939993842j),
                    (0.9745117519410239+0.13271463768855818j)]
r1 = (r1_sp+r1_sm)/2
theory = [(r1_sp+r1_sm)/2,r2,r3*2.75,r4]

#Atom-photon gate
simulation = [(0.939707+0.015902j),
(0.934679+0.000000j),
(0.968271+0.172339j),
(-0.398521-0.008054j)]

theory = [(0.939707+0.015902j),
(0.934679+0.000000j),
(0.468271+0.172339j),
(-0.398521-0.008054j)]

ideal = [1,1,1,-1]

r_coefficient_list = [ideal,theory,simulation]

def normalize(v):
    return(v / np.sqrt(np.sum(abs(v)**2)))

list_tt = []
for r in r_coefficient_list:
    
    
    """
    sp1_v2 = np.array([1,0,0,0])
    sp1_pi = np.array([0,1,0,0])
    sm1_v2 = np.array([0,0,1,0])
    sm1_pi = np.array([0,0,0,1])
    
    basis1 = [sp1_v2,sp1_pi,sm1_v2,sm1_pi]"""
    
    M1_r_exp = np.array([
         [r[0], 0, 0, 0],
         [0, r[1], 0, 0],
         [0, 0, r[2], 0],
         [0, 0, 0, r[3]]
         ])
    
    
    F2_v2 = np.array([1,0,0,0])
    F2_pi = np.array([0,1,0,0])
    F1_v2 = np.array([0,0,1,0])
    F1_pi = np.array([0,0,0,1])
    
    basis2 = [F2_v2,F2_pi,F1_v2,F1_pi]
    
    T = 1/np.sqrt(2)*np.array([
         [1j, 0, 1, 0],
         [0, 1j, 0, 1 ],
         [-1j, 0, 1, 0],
         [0, -1j, 0, 1]
         ])
    
    T = 1/np.sqrt(2)*np.array([
         [1j, -1j, 0, 0],
         [1, 1, 0, 0 ],
         [0, 0, 1j, -1j],
         [0, 0, 1, 1]
         ])
    
    T_inv = T.conjugate().transpose()
    
    np.set_printoptions(precision =3)
    
    #Transformation to basis2:
    M2_r_exp = np.matmul(T_inv,np.matmul(M1_r_exp,T))
    
    trueth_table = []
    for x in basis2:
        for y in basis2:
            x_r = np.matmul(M2_r_exp,x)
            x_r_norm = normalize(x_r)
            X=round(abs(np.matmul(y,x_r_norm))**2,3)
            trueth_table.append(X)
    
    tt = np.reshape(np.array(trueth_table), (4,4))
    list_tt.append(tt)
    print(tt)


#%% plotting
figure = plt.figure(figsize=(16,6))
figure.suptitle('Cnot gate', fontsize=25)
pos = 131
title = ['ideal','theory', 'simulation']
s = 0
for tt in list_tt[0:3]:
    
    axes = figure.add_subplot(pos)
    caxes = axes.matshow(tt, interpolation ='nearest', cmap = 'coolwarm', vmin = 0, vmax = 1)
    figure.colorbar(caxes)
    
    axes.tick_params(axis="x", bottom=True, top=False, labelbottom=True, labeltop=False)
    axes.title.set_text(title[s])
    ticklabels=[r'$\left| F2 \right\rangle_1 \left| L \right\rangle_2$',
              r'$\left| F2 \right\rangle_1 \left| R \right\rangle_2$',
              r'$\left| F1 \right\rangle_1 \left| L \right\rangle_2$',
              r'$\left| F1 \right\rangle_1 \left| R \right\rangle_2$']
    
    axes.set_xticklabels(['']+ticklabels)
    axes.set_yticklabels(['']+ticklabels)
    

    
    axes.set_xlabel('Output', fontsize=16, color='red')
    axes.set_ylabel('Input', fontsize=16, color='blue')
    
    for i in range(4):
        for j in range(4):
            c = tt[j,i]
            axes.text(i, j, str(c)[0:4], va='center', ha='center', fontsize = 14)
    pos = pos + 1
    s = s + 1
    #plt.xticks(rotation=-65)

plt.tight_layout()
plt.show()


