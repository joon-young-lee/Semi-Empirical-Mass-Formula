import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from modules import binding_energy
from modules import standard_deviation
from modules import binding_energy_only_5
import pandas as pd
import cgs

# Data
df = pd.read_fwf('BE_data.txt', usecols=(2,3,11),
              widths=(1,3,5,5,5,1,3,4,1,13,11,11,9,1,2,11,9,1,3,1,12,11,1),
              skiprows=39, header=None,
              index_col=False)
df = np.array(df)
N_data0 = np.array(df[:, 0])
Z_data0 = np.array(df[:, 1])

BE0 = np.array(df[:, 2])

for i in range(len(BE0)):
    if BE0[i][-1] == '#':
        BE0[i] = -1.0
N_data = np.array([])
Z_data = np.array([])
BE = np.array([])
for i in range(len(BE0)):
    if BE0[i] != -1.0:
        N_data = np.append(N_data, N_data0[i])
        Z_data = np.append(Z_data, Z_data0[i])
        BE = np.append(BE, BE0[i])
BE = BE.astype('float64')
Bexp_data = BE
# print(BE[200])


N_data = N_data[7:]
Z_data = Z_data[7:]
Bexp_data = Bexp_data[7:] * (N_data + Z_data) / 1000
# Combine Z and N into a 2D array for curve_fit
ZN_data = np.vstack((Z_data, N_data)).T

# Perform the least-squares fit
params_initial = [15, 17, 0.7, 11, 1.5, 1.0, 30, 20.0, -20.0, -1.0]  # Initial guesses for av, as_, ac, ap, kv, ks, W
# params_initial = [15, 17, 0.001, 1.5, 1.0]  # Initial guesses for av, as_, ac, kv, ks
params_opt, params_cov = curve_fit(lambda ZN, av, as_, ac, ap, kv, ks, W, ak, a0, fp
                                   : binding_energy(ZN[:, 0], ZN[:, 1], av, as_, ac, ap, kv, ks, W, ak, a0, fp), 
                                   ZN_data, Bexp_data, p0=params_initial, maxfev = int(1e5), method = 'lm')



# for Z = 40
N_test = 200
Z_val = 40 # Decide the Z value to see when binding energy decreases
Z = np.zeros(N_test)
N = np.zeros(N_test)
for j in range(N_test):
    Z[j] = Z_val
    N[j] = j + 1


BE_cal = binding_energy(Z, N, *params_opt)
total_mass = cgs.m_p * Z + cgs.m_n * N - binding_energy(Z, N, *params_opt)
delta_mass = np.zeros(N_test)
for k in range(N_test-1):
    delta_mass[k] = total_mass[k+1] - total_mass[k]
# print(delta_mass[0:-1])

print(len(BE_cal))
print(len(N))
plt.plot(N, -BE_cal)
# plt.plot(N[0:-1], delta_mass[0:-1])
plt.title(f'Z = {Z_val}')
# plt.axhline(cgs.m_n, color='k')
plt.xlabel('N varied')
plt.ylabel(r'$\Delta$ Mass (MeV)')
plt.show()

all_Z = np.linspace(8, 118, num=111, endpoint=True)
varied_N = np.linspace(1, 201, num=200, endpoint=False)
print(varied_N)
# find the allowed added neutrons
