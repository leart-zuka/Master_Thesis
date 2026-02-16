from rich import print


baseline = 1.4e-6

H_in = 6.727e-3 - baseline
V_in = 6.610e-3 - baseline
A_in = 6.562e-3 - baseline
D_in = 6.687e-3 - baseline
R_in = 6.655e-3 - baseline
L_in = 6.742e-3 - baseline

H_trans = 6.633e-3 - baseline
V_trans = 6.584e-3 - baseline
A_trans = 6.522e-3 - baseline
D_trans = 6.580e-3 - baseline
R_trans = 6.549e-3 - baseline
L_trans = 6.618e-3 - baseline

H_refl = 12.54e-6 - baseline
V_refl = 14.48e-6 - baseline
A_refl = 13.62e-6 - baseline
D_refl = 13.60e-6 - baseline
R_refl = 13.24e-6 - baseline
L_refl = 13.15e-6 - baseline


print("-------------------------------------")
print("-------Values for angle of 4 Deg-----")
print("-------------------------------------")

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


print("-----------------------------------------")
print("-------Values for angle of 6.703 Deg-----")
print("-----------------------------------------")
baseline = 1.4e-6

H_trans = 9.754e-3 - baseline
V_trans = 9.839e-3 - baseline
A_trans = 9.726e-3 - baseline
D_trans = 9.799e-3 - baseline
R_trans = 9.985e-3 - baseline
L_trans = 9.829e-3 - baseline

H_refl = 48.92e-6 - baseline
V_refl = 52.97e-6 - baseline
A_refl = 51.51e-6 - baseline
D_refl = 50.83e-6 - baseline
R_refl = 53.10e-6 - baseline
L_refl = 50.78e-6 - baseline


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


print("-----------------------------------------")
print("-------Values for angle of ~4 Deg-----")
print("-----------------------------------------")
baseline_reflection = 1.736e-6
baseline_transmission = 145.5e-9

H_trans = 6.651e-3 - baseline_transmission
V_trans = 6.699e-3 - baseline_transmission
A_trans = 6.689e-3 - baseline_transmission
D_trans = 6.649e-3 - baseline_transmission
R_trans = 6.658e-3 - baseline_transmission
L_trans = 6.705e-3 - baseline_transmission

H_refl = 115.0e-6 - baseline_reflection
V_refl = 114.7e-6 - baseline_reflection
A_refl = 115.7e-6 - baseline_reflection
D_refl = 114.2e-6 - baseline_reflection
R_refl = 114.9e-6 - baseline_reflection
L_refl = 115.5e-6 - baseline_reflection


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
