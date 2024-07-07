import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv('../data/data.csv')

plt.figure()
plt.plot(data['r'], data['P'])
plt.axhline(y=0, color='k')
plt.xlabel('Raggio')
plt.ylabel('Pressione')
plt.title('P(r)')
plt.savefig('../images/P(r).png')

plt.figure()
plt.plot(data['r'], data['m'])
plt.xlabel('Raggio')
plt.ylabel('Massa')
plt.title('m(r)')
plt.savefig('../images/m(r).png')

plt.figure()
plt.plot(data['r'], data['rho'])
plt.xlabel('Raggio')
plt.ylabel('Densità')
plt.title('rho(r)')
plt.savefig('../images/rho(r).png')

plt.show()
