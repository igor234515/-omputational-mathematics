import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D
from matplotlib import cm

Sub = []
data = []

angle = "D:\\YandexDisk-zhukov.ia@phystech.edu\\MIPT\\VI Semestr\\Count_Maths\\Count_Maths_ex_1\\XXVI.10.2\\Solution.txt"
laks = "D:\YandexDisk-zhukov.ia@phystech.edu\MIPT\VI Semestr\Count_Maths\Count_Maths_ex_1\XXVI.10.2\Laks.txt"
with open(laks,'r') as f:
    data = f.readlines()
    for i in range(len(data)):
        data[i] = np.array([float(item) for item in data[i].split()])
spacescale = np.arange(len(data[0]))
data = np.array(data)

S = np.linspace(0, 20, len(data[0]))
T = np.linspace(0, 18, len(data[:, 0]))
space, time = np.meshgrid(S, T)

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
ax.plot_surface(time, space, data)
ax.set_xlabel('time')
ax.set_ylabel('space')
plt.show()