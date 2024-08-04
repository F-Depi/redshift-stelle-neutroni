## ::setlocal makeprg=cd\ script\ &&\ python3.12\ plot.py
from matplotlib.lines import lineStyles
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

SMALL_SIZE = 13
MEDIUM_SIZE = 13
BIGGER_SIZE = 13

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

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


def plot_MR_rela(save=['yes', 'no']):
    plt.figure()

    data = np.genfromtxt('../data/MR_0.csv', delimiter=',', skip_header=1, dtype=float)
    R = data[:, 2]
    M = data[:, 3]
    plt.plot(R, R/2, linestyle='-', color='black', label='Raggio di Shwarzschild')
    plt.plot(R, M, linestyle='', marker='.', color='r', label = r'$\epsilon_1$')
    data = np.genfromtxt('../data/MR_1.csv', delimiter=',', skip_header=1, dtype=float)
    R = data[:, 2]
    M = data[:, 3]
    plt.plot(R, M, linestyle='', marker='.', color='b', label = r'$\epsilon_2$')
    data = np.genfromtxt('../data/MR_2.csv', delimiter=',', skip_header=1, dtype=float)
    R = data[:, 2]
    M = data[:, 3]
    plt.plot(R, M, linestyle='', marker='.', color='g', label = r'$\epsilon_3$')

    plt.title('Curva massa raggio per stelle con 3 equazioni di stato diverse')
    plt.xlabel('Raggio [km]')
    plt.ylabel(r'Massa [$M_\odot$]')
    plt.xscale('log')
    plt.yscale('log')
    plt.grid()
    plt.legend()
    plt.legend()
    if save == 'yes': plt.savefig('../report/Figures/MR.eps', format='eps')
    plt.show()


def plot_Phi(save=['yes', 'no']):
    colors = [['purple', '#32CD32'], ['orange', 'b'], ['r', 'cyan']]
    k = 0
    plt.figure()
    for i in ['0', '1', '2']:
        data = np.genfromtxt('../data/Phi_int_' + i + '.csv', delimiter=',', skip_header=1, dtype=float)
        r_int = data[:, 0]
        Phi_int = data[:, 1]

        data = np.genfromtxt('../data/Phi_ext_' + i + '.csv', delimiter=',', skip_header=1, dtype=float)
        r_ext = data[:, 0]
        Phi_ext = data[:, 1]

        plt.plot(r_int, Phi_int, color=colors[k][0], label=r'$\Phi_\text{int} (r)$ per $\epsilon_{'+str(k+1)+'}$')
        plt.plot(r_ext, Phi_ext, color=colors[k][1], label=r'$\Phi_\text{ext} (r)$ per $\epsilon_{'+str(k+1)+'}$')

        k += 1
    plt.title('Potenziale gravitazionale delle stelle con massa maggiore')
    plt.xlabel('Raggio [km]')
    plt.ylabel(r'$\Phi$')
    plt.grid()
    plt.legend()
    if save == 'yes': plt.savefig('../report/Figures/Phi.eps', format='eps')
    plt.show()


def plot_B(save=['yes', 'no']):
    M = [r'$14\text{M}_\odot$', r'$0.93\text{M}_\odot$', r'$1.5\text{M}_\odot$']
    R = [r'$59$km', r'$11$km', r'$8.6$km']
    kk = 0

    for tipo in ['1', '2', '3']:
        data = np.genfromtxt('../data/radianza/B.csv', delimiter=',', skip_header=1, dtype=float)
        data_corr1 = np.genfromtxt('../data/radianza/Bcorr_'+tipo+'_1.5.csv', delimiter=',', skip_header=1, dtype=float)
        data_corr2 = np.genfromtxt('../data/radianza/Bcorr_'+tipo+'_8.0.csv', delimiter=',', skip_header=1, dtype=float)
        data_corr3 = np.genfromtxt('../data/radianza/Bcorr_'+tipo+'_-1.0.csv', delimiter=',', skip_header=1, dtype=float)
        plt.figure()
        plt.plot(data[:,0], data[:,1], label=r'B$(\nu)$')
        plt.plot(data_corr1[:,0], data_corr1[:,1], label=rf'$B_{tipo}(\nu)$, $r = 1.5R$')
        plt.plot(data_corr2[:,0], data_corr2[:,1], color='purple', label=rf'$B_{tipo}(\nu)$, $r = 8R$')
        plt.plot(data_corr3[:,0], data_corr3[:,1], 'g--', label=rf'$B_{tipo}(\nu)$, $r = \infty$')
        #plt.ticklabel_format(axis='y', scilimits=(-3, -3))
        plt.title(rf'Radianza in funzione della frequenza, $M_{tipo} = ${M[kk]} $R_{tipo} = ${R[kk]}')
        plt.xlabel('f [Hz]')
        plt.ylabel(r'Radianza $[\frac{\text{MeV}}{\text{fm}^2}]$')
        plt.legend()

        kk += 1

        if save == 'yes': plt.savefig('../report/Figures/radianza'+tipo+'.eps', format='eps')
    plt.show()


def plt_P_test_cvgN(save=['yes', 'no']):
    dataT = np.genfromtxt('../data/potenza/test_cvg_N_trap.csv', delimiter=',', skip_header=1, dtype=float)
    dataS = np.genfromtxt('../data/potenza/test_cvg_N_simp.csv', delimiter=',', skip_header=1, dtype=float)

    plt.figure()

    metodo = ['T:   ', 'S:   ']
    kk = 0

    for data in [dataT, dataS]:
        P = abs(data[-1,0] - data[:-1,0]) / data[-1,0]
        N = data[:-1,1]

        plt.plot(N, P, linestyle='', marker='o', label=metodo[kk] + r'$\frac{| \mathcal{P} - \mathcal{P}_{N = 1e8} |}{\mathcal{P}_{N = 1e8}}$')
        kk += 1
    plt.title(r"Convergenza dell'integrale per $r = 1.5R$, $R = 59km$, $M = 14 M_\odot$")
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('log(N)')
    plt.ylabel(r'$\log_{10} \left( \frac{Potenza}{\text{MeVs}^{-1}\text{fm}^{-2}} \right)$')
    plt.legend(fontsize=20)
    plt.tight_layout()
    if save == 'yes': plt.savefig('../report/Figures/Pot_cvgN.eps', format='eps')
    plt.show()


def plt_P_test_cvgA(save=['yes', 'no']):
    dataT = np.genfromtxt('../data/potenza/test_cvg_A_trap.csv', delimiter=',', skip_header=1, dtype=float)
    dataS = np.genfromtxt('../data/potenza/test_cvg_A_simp.csv', delimiter=',', skip_header=1, dtype=float)

    plt.figure()

    metodo = ['T:   ', 'S:   ']
    kk = 0

    for data in [dataT, dataS]:
        P = abs(data[-1,0] - data[:-1,0]) / data[-1,0]
        nu_max = data[:-1,2]

        plt.plot(nu_max, P, linestyle='', marker='o', label=metodo[kk]+r'$\frac{| \mathcal{P} - \mathcal{P} _{\hat \nu_\text{max} = 200} |}{\mathcal{P} _{\hat \nu_\text{max} = 200}}$')
        kk += 1
    plt.axvline(20, color='r')
    plt.title(r"Convergenza dell'integrale per $r = 1.5R$, $R = 59km$, $M = 14 M_\odot$")
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'$\log(\hat \nu_\text{max}$)')
    plt.ylabel(r'$\log_{10} \left( \frac{Potenza}{\text{MeVs}^{-1}\text{fm}^{-2}} \right)$')
    plt.legend(fontsize=20)
    plt.tight_layout()
    if save == 'yes': plt.savefig('../report/Figures/Pot_cvgA.eps', format='eps')
    plt.show()


def plt_Pot(save=['yes','no']):

    plt.figure(figsize=(12,6))
    colors = ['orange', 'purple', 'black', 'g', 'r', 'b']

    kk = 0
    for metodo in ['simp', 'trap']:
        for i in ['1', '2', '3']:
            data = np.genfromtxt('../data/potenza/Pot_'+metodo+'_'+i+'.csv', delimiter=',', skip_header=1, dtype=float)
            r = data[:-1,0]
            Pot = data[:-1,1]
            # The last row contains R_star and Pot ad infinity
            R = data[-1,0]
            Pot_inf = data[-1,1]

            plt.plot(r, Pot, color=colors[kk], marker='.', linestyle='', label=metodo[0].upper()+': 'r'$\mathcal{P}_'+i+r'(r)$')
            plt.axvline(R, color=colors[kk], linestyle='-')
            plt.axhline(Pot_inf, color=colors[kk], linestyle='--')
            kk += 1

    plt.title('Potenza irradiata dalle 3 stelle più massive')
    plt.xlabel('raggio [km]')
    plt.ylabel(r'potenza $\left[ \frac{\text{MeV}}{\text{sfm}^3} \right]$')
    plt.legend(loc='upper right')
    plt.tight_layout()
    if save == 'yes': plt.savefig('../report/Figures/Pot.eps', format='eps')
    plt.show()


''' P(rho) vs rho(P) '''
#test_P_rhi()

''' Andamento P(r), m(r), rho(r) di una stella '''
#test_stella()

''' h diversi '''
#test_cvg()

''' Grafico M-R e altro'''
#plot_MR()

''' Grafico M-R generale '''
#plot_MR_rela()

''' Grafico potenziale gravitazionale '''
#plot_Phi()

''' Grafico della radianza '''
#plot_B()

''' Convergenza dell'integrale della potenza'''
#plt_P_test_cvgN('yes')
#plt_P_test_cvgA('yes')

''' Potenza in funzione di r '''
#plt_Pot('yes')


