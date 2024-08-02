## :setlocal makeprg=cd\ script\ &&\ python3.12\ radianza.py
import matplotlib.pyplot as plt
import numpy as np

h = 4.135668e-21 # MeV s
c = 2.9979245800e10 # cm/s

def funB(nu, h, c):
    return (2 * h * nu^3) / c^2 / (np.exp(h * nu) - 1)

nu = range(0, 1000, 1)
B = []
for ii in nu:
    B.append(funB(ii, h, c))

plt.figure()
plt.plot(nu, B)
plt.show()
