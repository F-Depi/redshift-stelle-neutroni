import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

''' Andamento P(r), m(r), rho(r) di una stella
data = pd.read_csv('../data/data.csv')

print("M = " + str(data['m'].iloc[-1]))
print("R = " + str(data['r'].iloc[-1]))

fig, axs = plt.subplots(1, 3, figsize=(20, 5))
axs[0].plot(data['r'], data['P'], linestyle='', marker='.', markersize=0.5)
axs[0].axhline(y=0, color='k')
#axs[0].set_ylim(-5, 200)
axs[0].set_xlabel('Raggio')
axs[0].set_ylabel('Pressione')
axs[0].set_title('P(r)')

axs[1].plot(data['r'], data['m'], linestyle='', marker='.', markersize=0.5)
#axs[1].set_ylim(-0.1, 1)
axs[1].set_xlabel('Raggio')
axs[1].set_ylabel('Massa')
axs[1].set_title('m(r)')

axs[2].plot(data['r'], data['rho'], linestyle='', marker='.', markersize=0.5)
#axs[2].set_ylim(-1, 15)
axs[2].set_xlabel('Raggio')
axs[2].set_ylabel('Densità')
axs[2].set_title('rho(r)')

plt.show()
'''

''' h diversi '''
data_cvg = pd.read_csv('../data/data_cvg.csv')

fig, axs = plt.subplots(1, 2, figsize=(15, 5))
axs[0].plot(data_cvg['h'], data_cvg['R'], linestyle='', marker='o')
axs[0].set_xscale('log')
axs[0].set_xlabel('Incremento')
axs[0].set_ylabel('Raggio')

axs[1].plot(data_cvg['h'], data_cvg['M'], linestyle='', marker='o')
axs[1].set_xscale('log')
axs[1].set_xlabel('Incremento')
axs[1].set_ylabel('Massa')

plt.show()


exit()
''' Grafico M-R '''
data3 = pd.read_csv('../data/RM2.csv')

plt.figure()
plt.plot(data3['R'], data3['M'], linestyle='', marker='o')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Raggio')
plt.ylabel('Massa')
plt.title('M-R')
plt.show()



