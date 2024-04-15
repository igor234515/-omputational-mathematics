import matplotlib.pyplot as plt
import numpy as np

pressurefile = "D:\YandexDisk-zhukov.ia@phystech.edu\MIPT\VI Semestr\Count_Maths\Count_Maths_ex_1\Riman\pressure.txt"
densefile = "D:\YandexDisk-zhukov.ia@phystech.edu\MIPT\VI Semestr\Count_Maths\Count_Maths_ex_1\Riman\dense.txt"
energyfile = "D:\YandexDisk-zhukov.ia@phystech.edu\MIPT\VI Semestr\Count_Maths\Count_Maths_ex_1\Riman\energy.txt"
velocityfile = "D:\YandexDisk-zhukov.ia@phystech.edu\MIPT\VI Semestr\Count_Maths\Count_Maths_ex_1\Riman\dense2.txt"

pressure_lex = "D:\YandexDisk-zhukov.ia@phystech.edu\MIPT\VI Semestr\Count_Maths\Count_Maths_ex_1\Empty\CIR_scheme_p.txt"
dense_lex = "D:\YandexDisk-zhukov.ia@phystech.edu\MIPT\VI Semestr\Count_Maths\Count_Maths_ex_1\Empty\CIR_scheme_rho.txt"
energy_lex = "D:\YandexDisk-zhukov.ia@phystech.edu\MIPT\VI Semestr\Count_Maths\Count_Maths_ex_1\Empty\CIR_scheme_e.txt"
velocity_lex = "D:\YandexDisk-zhukov.ia@phystech.edu\MIPT\VI Semestr\Count_Maths\Count_Maths_ex_1\Empty\CIR_scheme_u.txt"

u = []
P = []
ρ = []
e = []

with open(pressurefile, "r") as f:
    for line in f:
        P.append(float(line))

with open(densefile, "r") as f:
    for line in f:
        ρ.append(float(line))

with open(velocityfile, "r") as f:
    for line in f:
        u.append(float(line))

with open(energyfile, "r") as f:
    for line in f:
        e.append(float(line))
L = np.linspace(-10, 10, 100)

plt.subplot(2, 2, 1)
plt.scatter(L, P)
plt.grid()
plt.title("Pressure")

plt.subplot(2, 2, 2)
plt.scatter(L, ρ)
plt.grid()
plt.title("Density")

plt.subplot(2, 2, 3)
plt.scatter(L, u)
plt.grid()
plt.title("Velocity")

plt.subplot(2, 2, 4)
plt.scatter(L, e)
plt.grid()
plt.title("Energy")

plt.show()
