import numpy as np
import matplotlib.pyplot as plt

x, y, u = np.loadtxt("zad_5.dat", unpack=True)

plt.tricontourf(x, y, u, 30)
plt.colorbar()
plt.show()