import numpy as np
import matplotlib.pyplot as plt
pressure = "D:\YandexDisk-zhukov.ia@phystech.edu\MIPT\VI Semestr\Count_Maths\Count_Maths_ex_2\Viscous layer\Pressure.txt"
p = []
internal = []
with open(pressure, "r") as f:
    for i, line in enumerate(f):
        internal.append(float(line.rstrip()))
        if (i+1)%100 == 0:
            p.append(internal)
            internal = []

print(len(p))
x = np.linspace(0, 99, 100)

plt.subplot(2, 2, 1)
plt.scatter(x, p[0])
plt.grid()

plt.subplot(2, 2, 2)
plt.scatter(x, p[5])
plt.grid()

plt.subplot(2, 2, 3)
plt.scatter(x, p[11])
plt.grid()

plt.subplot(2, 2, 4)
plt.scatter(x, p[20])
plt.grid()


plt.show()
