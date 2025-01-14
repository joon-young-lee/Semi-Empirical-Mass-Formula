import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import ScalarFormatter
import numpy as np

BPS = np.loadtxt('BPS.txt', skiprows=1)
EoS_data = np.loadtxt('EoS.txt')
BPS_MassDensity = BPS[:44, 0]
BPS_Pressure = BPS[:44, 1]
MassDensity = EoS_data[:, 0]
Pressure = EoS_data[:, 1]
BPS_Z_array = BPS[:44, 3]
BPS_A_array = BPS[:44, 4]
Z_array = EoS_data[:, 3]
A_array = EoS_data[:, 4]

if __name__ == '__main__':
    mpl.rcParams['figure.dpi'] = 200
    # fig, axs = plt.subplots(2, 2, figsize=(12, 7))
    # fig.suptitle('EoS', fontsize=20)
    # axs[0, 0].scatter(BPS_Pressure, BPS_MassDensity, color = 'b')
    # axs[0, 1].scatter(Pressure, MassDensity, color='r')
    # axs[0, 0].set_xlabel(r'Pressure $\mathrm{erg/cm^3}$', fontsize = 10)
    # axs[0, 0].set_ylabel(r'Mass Density $\mathrm{erg/cm^3}$', fontsize = 10)
    # axs[0, 1].set_xlabel(r'Pressure $\mathrm{erg/cm^3}$', fontsize = 10)
    # axs[0, 1].set_ylabel(r'Mass Density $\mathrm{erg/cm^3}$', fontsize = 10)
    
    # axs[1, 0].scatter(BPS_Z_array, BPS_A_array, color='b')
    # axs[1, 0].set_xlabel(r'Proton Number', fontsize = 10)
    # axs[1, 0].set_ylabel(r'Mass Number', fontsize = 10)
    # axs[1, 1].scatter(Z_array, A_array, color='r')
    # axs[1, 1].set_xlabel(r'Proton Number', fontsize = 10)
    # axs[1, 1].set_ylabel(r'Mass Number', fontsize = 10)
    # plt.show()
    
    fig, ax = plt.subplots(figsize=(10, 5.5))
    fontsize = 15
    # Plot the data
    ax.plot(BPS_Pressure, BPS_MassDensity, '-', color='b', label='BPS')
    ax.plot(Pressure, MassDensity, '-', color='r', label='SEMF')

    # Add labels
    ax.set_xlabel(r'Pressure $\mathrm{[erg/cm^3]}$', fontsize=fontsize)
    ax.set_ylabel(r'Energy Density $\mathrm{[erg/cm^3]}$', fontsize=fontsize)

    # Add scientific notation to ticks
    ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(axis='both', style='sci', scilimits=(0, 0))

    # Customize ticks
    ax.tick_params(axis='both', which='major', labelsize=fontsize)

    # Add legend
    ax.legend(loc='upper left', fontsize=20)

    # Add grid
    ax.grid()

    # Show plot
    # plt.savefig('../../../Desktop/research/presentation_1226/beamer/eos_crust.png')
    plt.savefig('./eos_crust.png')
    plt.show()