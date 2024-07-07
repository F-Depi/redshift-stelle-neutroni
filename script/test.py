import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

data = pd.read_csv('../data/rho-P.csv')
data2 = pd.read_csv('../data/P-rho.csv')

plt.plot(data['rho'], data['P'], linestyle='', marker='.')
plt.plot(data2['rho'], data2['P'], linestyle='', marker='.')
plt.axhline(y=0, color='k')
plt.xlabel('Density')
plt.ylabel('Pressure')
plt.show()
