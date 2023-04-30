import matplotlib.pyplot as plt
import numpy as np

pressure_lex = "D:\YandexDisk-zhukov.ia@phystech.edu\MIPT\VI Semestr\Count_Maths\Count_Maths_ex_1\Empty\CIR_scheme_p.txt"
dense_lex = "D:\YandexDisk-zhukov.ia@phystech.edu\MIPT\VI Semestr\Count_Maths\Count_Maths_ex_1\Empty\CIR_scheme_rho.txt"
energy_lex = "D:\YandexDisk-zhukov.ia@phystech.edu\MIPT\VI Semestr\Count_Maths\Count_Maths_ex_1\Empty\CIR_scheme_e.txt"
velocity_lex = "D:\YandexDisk-zhukov.ia@phystech.edu\MIPT\VI Semestr\Count_Maths\Count_Maths_ex_1\Empty\CIR_scheme_u.txt"

p = []
pr = []
u = []
ur = []
ρ = []
ρr = []
e = []
er = []

with open(pressure_lex, "r") as f:
    for line in f:
        p.append(line.split())
for el in p[0]:
    el = float(el)
    pr.append(el)

with open(dense_lex, "r") as f:
    for line in f:
        ρ.append(line.split())
for el in ρ[0]:
    el = float(el)
    ρr.append(el)

with open(energy_lex, "r") as f:
    for line in f:
        e.append(line.split())
for el in e[0]:
    el = float(el)
    er.append(el)

with open(velocity_lex, "r") as f:
    for line in f:
        u.append(line.split())
for el in u[0]:
    el = float(el)
    ur.append(el)

L = np.linspace(-10, 10, 100)

plt.subplot(2, 2, 1)
plt.scatter(L, pr)
plt.grid()
plt.title("Pressure")

plt.subplot(2, 2, 2)
plt.scatter(L, ρr)
plt.grid()
plt.title("Density")

plt.subplot(2, 2, 3)
plt.scatter(L, ur)
plt.grid()
plt.title("Velocity")

plt.subplot(2, 2, 4)
plt.scatter(L, er)
plt.grid()
plt.title("Energy")

plt.show()