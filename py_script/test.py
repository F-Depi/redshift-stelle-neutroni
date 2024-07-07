import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

data = pd.read_csv('../data/rho-P.csv')
data2 = pd.read_csv('../data/P-rho.csv')

plt.plot(data['rho'], data['P'])
plt.plot(data2['rho'], data2['P'])
plt.xlabel('Density')
plt.ylabel('Pressure')
plt.show()
