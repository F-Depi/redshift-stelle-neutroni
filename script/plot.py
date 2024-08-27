## ::setlocal makeprg=cd\ script\ &&\ python3.12\ plot.py
from matplotlib.lines import lineStyles
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

SMALL_SIZE = 13
MEDIUM_SIZE = 18
BIGGER_SIZE = 18

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)    # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the x tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the y tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize

def test_rhovsP():
    data = pd.read_csv('../data/test/rho_of_P.csv')
    data2 = pd.read_csv('../data/test/P_of_rho.csv')

    plt.plot(data['rho'], data['P'], linestyle='', marker='.', label=r'($\rho$, P_of_rho()')
    plt.plot(data2['rho'], data2['P'], linestyle='', marker='.', label='(findRho(P), P)')
    plt.axhline(y=0, color='k')
    plt.xlabel('Density')
    plt.ylabel('Pressure')
    plt.legend()
    plt.show()


def compare_eneries():
    plt.figure()

    labels = [r'Modello $a \alpha b \beta$',
              r'$\lambda = 5/3, K = 0.05$',
              r'$\lambda = 2.54, K = 0.01$']
    for i in [0, 1, 2]:
        data = np.genfromtxt('../data/test/E_of_P_'+str(i+1)+'.csv',
                             delimiter=',', skip_header=1, dtype=float)
        P = data[:,0]
        E = data[:,1]
        plt.plot(P, E, label=labels[i])

    plt.title('Confronto tra le densità di energia E(P)')
    plt.xlabel(r'P $\left[ \frac{Mev}{fm^3} \right]$')
    plt.ylabel(r'E $\left[ \frac{Mev}{fm^3} \right]$')
    plt.legend()
    plt.tight_layout()
    plt.show()


def test_cvg_stella():

    labels = [r'Modello $a \alpha b \beta$',
              r'$\lambda = 5/3, K = 0.05$',
              r'$\lambda = 2.54, K = 0.01$']

    fig, axs = plt.subplots(1, 3, figsize=(20, 5))

    for kk in [1, 2, 3]:

        data = pd.read_csv(f'../data/test/data{kk:d}.csv')

        print("M = " + str(data['m'].iloc[-1]))
        print("R = " + str(data['r'].iloc[-1]))

        axs[0].plot(data['r'], data['P'],
                    linestyle='', marker='.', label=labels[kk-1])
        axs[1].plot(data['r'], data['m'],
                    linestyle='', marker='.', label=labels[kk-1])
        axs[2].plot(data['r'], data['rho'],
                    linestyle='', marker='.', label=labels[kk-1])

    axs[0].axhline(y=0, color='k')
    axs[0].set_xlabel('Raggio [km]')
    axs[0].set_ylabel(r'Pressione $\left[ \frac{Mev}{fm^3} \right]$')
    axs[0].set_title('P(r)')

    axs[1].set_xlabel('Raggio [km]')
    axs[1].set_ylabel(r'Massa $[M_\odot]$')
    axs[1].set_title('m(r)')

    axs[2].set_xlabel('Raggio [km]')
    axs[2].set_ylabel(r'Densità $\left[ \frac{Mev}{fm^3} \right]$')
    axs[2].set_title('rho(r)')

    plt.tight_layout()
    axs[0].legend()
    axs[1].legend()
    axs[2].legend()

    plt.show()


def test_cvg_h():
    data_cvg = pd.read_csv('../data/test/cvg_1.csv')

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


def plot_MR_PR_MP():

    labels = [r'Modello $a \alpha b \beta$',
              r'$\lambda = 5/3, K = 0.05$',
              r'$\lambda = 2.54, K = 0.01$']

    fig, axs = plt.subplots(3, 1, figsize=(18, 10))

    for kk in [1, 2, 3]:
        data = np.genfromtxt(f'../data/MR_{kk:d}.csv', delimiter=',',
                             skip_header=1, dtype=float)
        P0 = data[:, 1]
        R = data[:, 2]
        M = data[:, 3]

        axs[0].plot(R, M, linestyle='', marker='.', label=labels[kk - 1])
        axs[1].plot(R, P0, linestyle='', marker='.', label=labels[kk - 1])
        axs[2].plot(P0, M, linestyle='', marker='.', label=labels[kk - 1])

    ## Shwarzschild radius
    R = np.arange(3, 40, 0.01)
    axs[0].plot(R, R/2, linestyle='-', color='red', label='Shwarzschild Radius')

    ## Make it nicer

    axs[0].set_xscale('log')
    axs[0].set_yscale('log')
    axs[0].grid('on')
    axs[0].set_xlabel('Raggio [km]')
    axs[0].set_ylabel('Massa [M_sun]')
    axs[0].set_ylabel(r'Massa $[M_{\odot}]$')
    axs[0].legend(loc='upper left')

    axs[1].set_xscale('log')
    axs[1].set_yscale('log')
    axs[1].grid('on')
    axs[1].set_xlabel('Raggio [km]')
    axs[1].set_ylabel(r'Pressione iniziale $\left[ \frac{Mev}{fm^3} \right]$')
    axs[1].legend(loc='upper left')

    axs[2].set_xscale('log')
    axs[2].set_yscale('log')
    #axs[2].set_xlim(-2, 1e18)
    axs[2].grid('on')
    axs[2].set_xlabel(r'Pressione iniziale $\left[ \frac{Mev}{fm^3} \right]$')
    axs[2].set_ylabel(r'Massa $[M_{\odot}]$')
    axs[2].legend(loc='upper left')

    plt.tight_layout()
    plt.show()


def plot_MR(save=['yes', 'no']):
    labels = [r'$\epsilon_1$: Modello $a \alpha b \beta$',
              r'$\epsilon_2:~\lambda = 5/3, K = 0.05$',
              r'$\epsilon_3:~\lambda = 2.54, K = 0.01$']

    plt.figure(figsize=(9,5))

    for kk in [1, 2, 3]:
        data = np.genfromtxt(f'../data/MR_{kk:d}.csv', delimiter=',',
                             skip_header=1, dtype=float)
        R = data[:, 2]
        M = data[:, 3]
        plt.plot(R, M, linestyle='', marker='.', label=labels[kk - 1])

    R = np.arange(3, 40, 0.01)
    plt.plot(R, R/2, linestyle='-', color='black', label='Raggio di Shwarzschild')

    plt.title('Curva massa raggio per stelle con 3 equazioni di stato diverse')
    plt.xlabel('Raggio [km]')
    plt.ylabel(r'Massa [$M_\odot$]')
    plt.xscale('log')
    plt.yscale('log')
    plt.grid()
    plt.tight_layout()
    plt.legend()
    if save == 'yes': plt.savefig('../report/Figures/MR.eps', format='eps')
    plt.show()


def plot_Phi(save=['yes', 'no']):

    colors = [['purple', '#32CD32'], ['orange', 'b'], ['r', 'cyan']]

    plt.figure(figsize=(9,5))

    for kk in [1, 2, 3]:
        data = np.genfromtxt(f'../data/Phi_int_{kk:d}.csv', delimiter=',',
                             skip_header=1, dtype=float)
        r_int = data[:, 0]
        Phi_int = data[:, 1]

        data = np.genfromtxt(f'../data/Phi_ext_{kk:d}.csv', delimiter=',',
                             skip_header=1, dtype=float)
        r_ext = data[:, 0]
        Phi_ext = data[:, 1]

        plt.plot(r_int, Phi_int, color=colors[kk-1][0],
                 label=r'$\Phi_\text{int} (r)$ per $\epsilon_{'+str(kk)+'}$')
        plt.plot(r_ext, Phi_ext, color=colors[kk-1][1],
                 label=r'$\Phi_\text{ext} (r)$ per $\epsilon_{'+str(kk)+'}$')

    plt.title('Potenziale gravitazionale delle stelle con massa maggiore')
    plt.xlabel('Raggio [km]')
    plt.ylabel(r'$\Phi$')
    plt.grid()
    plt.legend()
    if save == 'yes': plt.savefig('../report/Figures/Phi.eps', format='eps')
    plt.show()


def plot_B(save=['yes', 'no']):
    M = [r'$2.5\text{M}_\odot$', r'$0.99\text{M}_\odot$', r'$1.6\text{M}_\odot$']
    R = [r'$11.0$km', r'$10.8$km', r'$8.5$km']
    kk = 0

    for tipo in ['1', '2', '3']:
        data = np.genfromtxt('../data/radianza/B.csv',
                                   delimiter=',', skip_header=1, dtype=float)
        data_corr1 = np.genfromtxt('../data/radianza/Bcorr_'+tipo+'_1.5.csv',
                                   delimiter=',', skip_header=1, dtype=float)
        data_corr2 = np.genfromtxt('../data/radianza/Bcorr_'+tipo+'_8.0.csv',
                                   delimiter=',', skip_header=1, dtype=float)
        data_corr3 = np.genfromtxt('../data/radianza/Bcorr_'+tipo+'_-1.0.csv',
                                   delimiter=',', skip_header=1, dtype=float)

        plt.figure()
        plt.plot(data[:,0], data[:,1], label=r'B$(\nu)$')
        plt.plot(data_corr1[:,0], data_corr1[:,1],
                 label=rf'$B_{tipo}(\nu)$, $r = 1.5R$')
        plt.plot(data_corr2[:,0], data_corr2[:,1],
                 color='purple', label=rf'$B_{tipo}(\nu)$, $r = 8R$')
        plt.plot(data_corr3[:,0], data_corr3[:,1], 'g--',
                 label=rf'$B_{tipo}(\nu)$, $r = \infty$')
        plt.title(rf'Radianza stella $M_{tipo} = ${M[kk]} $R_{tipo} = ${R[kk]}')
        plt.xlabel('f [Hz]')
        plt.ylabel(r'Radianza $[\frac{\text{MeV}}{\text{fm}^2}]$')
        plt.tight_layout()
        plt.legend()

        kk += 1

        if save == 'yes':
            plt.savefig('../report/Figures/radianza'+tipo+'.eps', format='eps')
    plt.show()


def plt_P_test_cvgN(save=['yes', 'no']):
    dataT = np.genfromtxt('../data/potenza/test_cvg_N_trap.csv',
                          delimiter=',', skip_header=1, dtype=float)
    dataS = np.genfromtxt('../data/potenza/test_cvg_N_simp.csv',
                          delimiter=',', skip_header=1, dtype=float)

    plt.figure()

    metodo = ['T:   ', 'S:   ']
    lab = r'$\frac{|\mathcal{I} - \mathcal{I}_{N =1e8}|}{\mathcal{I}_{N =1e8}}$'
    kk = 0

    for data in [dataT, dataS]:
        P = abs(data[-1,0] - data[:-1,0]) / data[-1,0]
        N = data[:-1,1]

        plt.plot(N, P, linestyle='', marker='o', label=metodo[kk] + lab)
        kk += 1
    plt.title('Errore sul metodo di integrazione')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('log(N)')
    plt.ylabel(r'$\log_{10} \left( \text{Errore relativo} \right)$')
    plt.legend(fontsize=20)
    plt.tight_layout()
    if save == 'yes':
        plt.savefig('../report/Figures/Pot_cvgN.eps', format='eps')
    plt.show()


def plt_P_test_cvgA(save=['yes', 'no']):
    dataT = np.genfromtxt('../data/potenza/test_cvg_A_trap.csv',
                          delimiter=',', skip_header=1, dtype=float)
    dataS = np.genfromtxt('../data/potenza/test_cvg_A_simp.csv',
                          delimiter=',', skip_header=1, dtype=float)

    plt.figure()

    lab_num = r'| \mathcal{I} - \mathcal{I} _{\hat \nu_\text{max} = 200} |'
    lab_den = r'\mathcal{I} _{\hat \nu_\text{max} = 200}'
    lab = r'$\frac{' + lab_num + r'}{' + lab_den + r'}$'
    metodo = ['T:   ' + lab, 'S:   ' + lab]
    kk = 0

    for data in [dataT, dataS]:
        P = abs(data[-1,0] - data[:-1,0]) / data[-1,0]
        nu_max = data[:-1,2]

        plt.plot(nu_max, P, linestyle='', marker='o', label = metodo[kk])
        kk += 1
    plt.axvline(20, color='r')
    plt.title('Errore sul metodo di integrazione')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'$\log(\hat \nu_\text{max}$)')
    plt.ylabel(r'$\log_{10} \left( \text{Errore relativo} \right)$')
    plt.legend(fontsize=20)
    plt.tight_layout()
    if save == 'yes':
        plt.savefig('../report/Figures/Pot_cvgA.eps', format='eps')
    plt.show()


def plt_Pot(save=['yes','no']):

    plt.figure(figsize=(12,6))
    colors = ['g', 'r', 'b']

    kk = 0
    for i in ['1', '2', '3']:
        data = np.genfromtxt('../data/potenza/Pot_trap_'+i+'.csv', delimiter=',', skip_header=1, dtype=float)
        r = data[:-1,0]
        Pot = data[:-1,1]# * 1.602176e-13 * 1e30     # MeV to J and fm^-2 to m^-2
        # The last row contains R_star and Pot ad infinity
        R = data[-1,0]
        Pot_inf = data[-1,1]# * 1.602176e-13 * 1e30  # MeV to J and fm^-2 to m^-2

        plt.plot(r, Pot, color=colors[kk], marker='.', linestyle='', label=r'$\mathcal{P}_'+i+r'(r)$')
        plt.axvline(R, color=colors[kk], linestyle='-')
        plt.axhline(Pot_inf, color=colors[kk], linestyle='--')
        kk += 1

    plt.title('Potenza irradiata dalle 3 stelle più massive')
    plt.xlabel('raggio [km]')
    plt.ylabel(r'potenza $\left[ \frac{\text{MeV}}{\text{s fm}^2} \right]$')
    plt.legend(loc='upper right')
    plt.tight_layout()
    if save == 'yes': plt.savefig('../report/Figures/Pot.eps', format='eps')
    plt.show()


def plt_Teff(save=['yes','no']):

    plt.figure(figsize=(12,6))
    colors = ['g', 'r', 'b']

    kk = 0
    for i in ['1', '2', '3']:
        data = np.genfromtxt('../data/potenza/Teff_'+i+'.csv', delimiter=',', skip_header=1, dtype=float)
        T = data[:,0]
        Teff = data[:,1]

        plt.plot(T, Teff, color=colors[kk], marker='.', linestyle='', label=r'$T_{eff} (r = \infty, T)$, stella '+i)
        kk += 1
    plt.plot(T, T, color = 'black', linestyle='--', label='T = T')

    plt.title('Temperatura efficace dalle 3 stelle più massive')
    plt.xlabel('Temperatura Propria [MeV]')
    plt.ylabel('Temperatura all\'infinito [MeV]')
    plt.legend(loc='upper left', fontsize=16)
    plt.tight_layout()
    if save == 'yes': plt.savefig('../report/Figures/Teff.eps', format='eps')
    plt.show()


def plt_Teff_su_Pressione(T0, save=['yes','no']):
    plt.figure()
    #colors = ['g', 'r', 'b']
    lab = r'$T_{eff} (r = \infty, T = 1 \text{MeV}, P), ~ \epsilon_'

    for kk in [1, 2, 3]:
        data = np.genfromtxt(f'../data/punto7/Teff{T0:.2f}_su_P_{kk:d}.csv',
                             delimiter=',', skip_header=1, dtype=float)
        P = data[:,0]
        Teff = data[:,1]

        plt.plot(P, Teff, marker='.', linestyle='', label = lab+str(kk)+r'$')

    plt.title('Temperatura efficacie - pressione centrale\n'
              rf'stelle con $k_b T^* = {T0:.2f}$MeV')
    plt.xlabel(r'Pressione $\left[ \frac{\text{MeV}}{\text{fm}^3} \right]$')
    plt.ylabel(r'Temperatura $\left[ \text{MeV} \right]$')
    plt.legend(loc='upper right')
    plt.tight_layout()
    if save == 'yes':
        plt.savefig(f'../report/Figures/Teff{T0:.2f}_su_P.eps', format='eps')
    plt.show()




## I dati sono generati con gli script in C indicati

#################################### test.c ####################################
''' P(rho) vs rho(P) '''
#test_rhovsP()
#compare_eneries()

''' Andamento P(r), m(r), rho(r) per 3 politropiche a 1 pressione '''
#test_cvg_stella()

''' h diversi '''
#test_cvg_h()


#################################### main.c ####################################
''' Grafico M-R P-R M-P '''
#plot_MR_PR_MP()

''' Grafico M-R generale '''
#plot_MR('yes')

''' Grafico potenziale gravitazionale '''
#plot_Phi('yes')


################################## radianza.c ##################################
''' Grafico della radianza '''
#plot_B('no')

''' Convergenza dell'integrale della potenza'''
#plt_P_test_cvgN('yes')
#plt_P_test_cvgA('yes')

''' Potenza in funzione di r '''
#plt_Pot('yes')

''' Temperatura efficace '''
#plt_Teff('yes')


################################### punto7.c ###################################
''' Temperatura efficace rispetto a Pressione centrale delle stelle '''
#plt_Teff_su_Pressione(1, 'yes')
#plt_Teff_su_Pressione(0.01, 'yes')


