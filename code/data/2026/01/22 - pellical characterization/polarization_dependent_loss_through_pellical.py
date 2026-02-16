from rich import print

print("-----------------------------------------")
print("-------Values for angle of ~4 Deg-----")
print("-----------------------------------------")
baseline_reflection = 1.736e-6
baseline_transmission = 145.5e-9

H_trans = 6.608e-3 - baseline_transmission
V_trans = 6.609e-3 - baseline_transmission
A_trans = 6.615e-3 - baseline_transmission
D_trans = 6.583e-3 - baseline_transmission
R_trans = 6.658e-3 - baseline_transmission
L_trans = 6.586e-3 - baseline_transmission

H_refl = 114.4e-6 - baseline_reflection
V_refl = 113.3e-6 - baseline_reflection
A_refl = 115.2e-6 - baseline_reflection
D_refl = 114.0e-6 - baseline_reflection
R_refl = 114.9e-6 - baseline_reflection
L_refl = 113.4e-6 - baseline_reflection


H_in = H_trans + H_refl
V_in = V_trans + V_refl
A_in = A_trans + A_refl
D_in = D_trans + D_refl
R_in = R_trans + R_refl
L_in = L_trans + L_refl

print(f"Transmission for H: {H_trans / H_in * 100:.2f}%")
print(f"Transmission for V: {V_trans / V_in * 100:.2f}%")
print(f"Transmission for A: {A_trans / A_in * 100:.2f}%")
print(f"Transmission for D: {D_trans / D_in * 100:.2f}%")
print(f"Transmission for R: {R_trans / R_in * 100:.2f}%")
print(f"Transmission for L: {L_trans / L_in * 100:.2f}%")

print("-------------------------------------")

print(f"Reflection for H: {H_refl / H_in * 100:.2f}%")
print(f"Reflection for V: {V_refl / V_in * 100:.2f}%")
print(f"Reflection for A: {A_refl / A_in * 100:.2f}%")
print(f"Reflection for D: {D_refl / D_in * 100:.2f}%")
print(f"Reflection for R: {R_refl / R_in * 100:.2f}%")
print(f"Reflection for L: {L_refl / L_in * 100:.2f}%")
