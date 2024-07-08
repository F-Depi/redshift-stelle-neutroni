import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def test_P_rhi():
    data = pd.read_csv('../data/rho-P.csv')
    data2 = pd.read_csv('../data/P-rho.csv')

    plt.plot(data['rho'], data['P'], linestyle='', marker='.')
    plt.plot(data2['rho'], data2['P'], linestyle='', marker='.')
    plt.axhline(y=0, color='k')
    plt.xlabel('Density')
    plt.ylabel('Pressure')
    plt.show()


def test_stella():
    data = pd.read_csv('../data/data0.csv')

    print("M = " + str(data['m'].iloc[-1]))
    print("R = " + str(data['r'].iloc[-1]))

    fig, axs = plt.subplots(1, 3, figsize=(20, 5))
    axs[0].plot(data['r'], data['P'], linestyle='', marker='.')
    axs[0].axhline(y=0, color='k')
    #axs[0].set_ylim(-5, 200)
    axs[0].set_xlabel('Raggio')
    axs[0].set_ylabel('Pressione')
    axs[0].set_title('P(r)')

    axs[1].plot(data['r'], data['m'], linestyle='', marker='.')
    #axs[1].set_ylim(-0.1, 1)
    axs[1].set_xlabel('Raggio')
    axs[1].set_ylabel('Massa')
    axs[1].set_title('m(r)')

    axs[2].plot(data['r'], data['rho'], linestyle='', marker='.')
    #axs[2].set_ylim(-1, 15)
    axs[2].set_xlabel('Raggio')
    axs[2].set_ylabel('Densità')
    axs[2].set_title('rho(r)')

    data = pd.read_csv('../data/data1.csv')

    print("M = " + str(data['m'].iloc[-1]))
    print("R = " + str(data['r'].iloc[-1]))

    axs[0].plot(data['r'], data['P'], linestyle='', marker='.')
    axs[1].plot(data['r'], data['m'], linestyle='', marker='.')
    axs[2].plot(data['r'], data['rho'], linestyle='', marker='.')

    data = pd.read_csv('../data/data2.csv')

    print("M = " + str(data['m'].iloc[-1]))
    print("R = " + str(data['r'].iloc[-1]))

    axs[0].plot(data['r'], data['P'], linestyle='', marker='.')
    axs[1].plot(data['r'], data['m'], linestyle='', marker='.')
    axs[2].plot(data['r'], data['rho'], linestyle='', marker='.')

    plt.show()


def test_cvg():
    data_cvg = pd.read_csv('../data/data_cvg.csv')

    fig, axs = plt.subplots(1, 2, figsize=(15, 5))
    axs[0].plot(data_cvg['h'], data_cvg['R'], linestyle='', marker='.')
    axs[0].set_xscale('log')
    axs[0].set_xlabel('Incremento')
    axs[0].set_ylabel('Raggio')

    axs[1].plot(data_cvg['h'], data_cvg['M'], linestyle='', marker='.')
    axs[1].set_xscale('log')
    axs[1].set_xlabel('Incremento')
    axs[1].set_ylabel('Massa')

    plt.show()


def plot_MR():
    data = np.genfromtxt('../data/MR_0.csv', delimiter=',', skip_header=1, dtype=float)
    P0 = data[:, 1] * 1.733e15
    R = data[:, 2]
    M = data[:, 3]

    fig, axs = plt.subplots(3, 1, figsize=(18, 10))

    axs[0].plot(R, M, linestyle='', marker='.')
    axs[0].plot(R, R/2, linestyle='-', color='red')
    axs[1].plot(R, P0, linestyle='', marker='.')
    axs[2].plot(P0, M, linestyle='', marker='.')

    data = np.genfromtxt('../data/MR_1.csv', delimiter=',', skip_header=1, dtype=float)
    P0 = data[:, 1] * 1.733e15
    R = data[:, 2]
    M = data[:, 3]

    axs[0].plot(R, M, linestyle='', marker='.')
    axs[1].plot(R, P0, linestyle='', marker='.')
    axs[2].plot(P0, M, linestyle='', marker='.')

    data = np.genfromtxt('../data/MR_2.csv', delimiter=',', skip_header=1, dtype=float)
    P0 = data[:, 1] * 1.733e15
    R = data[:, 2]
    M = data[:, 3]

    axs[0].plot(R, M, linestyle='', marker='.')
    axs[1].plot(R, P0, linestyle='', marker='.')
    axs[2].plot(P0, M, linestyle='', marker='.')

    ## Make it nicer

    axs[0].set_xscale('log')
    axs[0].set_yscale('log')
    axs[0].grid('on')
    axs[0].set_xlabel('Raggio [km]')
    axs[0].set_ylabel('Massa [M_sun]')
    axs[0].set_ylabel(r'Massa $[M_{\odot}]$')
    axs[0].legend([r'Modello $a \alpha b \beta$', 'Shwarzschild Radius', r'$\lambda = 5/3, K = 0.05$', r'$\lambda = 2.54, K = 0.01$'])

    axs[1].set_xscale('log')
    axs[1].set_yscale('log')
    axs[1].grid('on')
    axs[1].set_xlabel('Raggio [km]')
    # axs[1].set_ylabel(r'Pressione iniziale $\left[\frac{Mev}{c^2 fm^3}\right]$')
    axs[1].set_ylabel(r'Pressione iniziale $\left[\frac{Kg}{m^3}\right]$')

    axs[2].set_xscale('log')
    axs[2].set_yscale('log')
    #axs[2].set_xlim(-2, 1e18)
    axs[2].grid('on')
    #axs[2].set_xlabel(r'Pressione iniziale $\left[\frac{Mev}{c^2 fm^3}\right]$')
    axs[2].set_xlabel(r'Pressione iniziale $\left[\frac{Kg}{m^3}\right]$')
    axs[2].set_ylabel(r'Massa $[M_{\odot}]$')

    plt.show()


''' P(rho) vs rho(P) '''
#test_P_rhi()

''' Andamento P(r), m(r), rho(r) di una stella '''
#test_stella()

''' h diversi '''
#test_cvg()

''' Grafico M-R '''
plot_MR()

''' Grafico M-R generale '''

