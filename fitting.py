import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from modules import binding_energy
from modules import standard_deviation
from modules import binding_energy_only_5
import pandas as pd

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
print(BE[200])


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
# params_opt, params_cov = curve_fit(lambda ZN, av, as_, ac, kv, ks
#                                    : binding_energy_only_5(ZN[:, 0], ZN[:, 1], av, as_, ac, kv, ks), 
#                                    ZN_data, Bexp_data, p0=params_initial, maxfev = int(1e5))
Bfit = binding_energy(Z_data, N_data, *params_opt)
# # standard deviation
sigma = standard_deviation(Bfit, Bexp_data)

# # Output the fitted coefficients
print("Fitted coefficients:")
print(f"av: {params_opt[0]}")
print(f"as: {params_opt[1]}")
print(f"ac: {params_opt[2]}")
print(f"ap: {params_opt[3]}")
print(f"kv: {params_opt[4]}")
print(f"ks: {params_opt[5]}")
print(f"W: {params_opt[6]}")
print(f"ak: {params_opt[7]}")
print(f"a0: {params_opt[8]}")
print(f"fp: {params_opt[9]}")
print(f"Standard Deviation: {sigma}")


# # Plot experimental vs. fitted values

fig, axs = plt.subplots(1, 2, figsize=(5, 5))
fig.suptitle('Empirical Data and SEMF', fontsize=20)
axs[0].scatter(Z_data + N_data, Bexp_data, label='Empirical', color='red', s=1)
axs[1].scatter(Z_data + N_data, Bfit, label='SEMF', color='blue', s=1)
axs[0].set_xlabel('Mass Number A (Z + N)', fontsize = 14)
axs[0].set_ylabel('Binding Energy (MeV)', fontsize = 14)
axs[1].set_xlabel('Mass Number A (Z + N)', fontsize = 14)
axs[1].set_ylabel('Binding Energy (MeV)', fontsize = 14)
axs[0].text(0, 1700, f'Number of data: {len(Z_data)}', fontsize = 15)
axs[1].text(-50, 1000, f"\n \
            $a_v$: {params_opt[0]:.2f}\n \
            $a_s$: {params_opt[1]:.2f}\n \
            $a_c$: {params_opt[2]:.2f}\n \
            $a_p$: {params_opt[3]:.2f}\n \
            $k_v$: {params_opt[4]:.2f}\n \
            $k_s$: {params_opt[5]:.2f}\n \
            $W$: {params_opt[6]:.2f}\n \
            $a_k$: {params_opt[7]:.2f}\n \
            $a_0$: {params_opt[8]:.2f}\n \
            $f_p$: {params_opt[9]:.2f}\n \
            $\sigma$: {sigma:.2f}\n", fontsize=13)
plt.savefig('./results/BE_data_1027.png')
axs[0].legend(fontsize = 15, loc='lower right')
axs[1].legend(fontsize = 15, loc='lower right')
plt.show()

#########################################################################
print(len(Z_data))
# for Z = 40
N_test = 80
Z_val = 40 # Decide the Z value to see when binding energy decreases
Z = np.zeros(N_test)
N = np.zeros(N_test)
for j in range(N_test):
    Z[j] = Z_val
    N[j] = j + 1
BE_cal = binding_energy(Z, N, *params_opt)


print(BE_cal)

# plt.plot(N, BE_cal)
# plt.xlabel('N varied')
# plt.ylabel('BE')
# plt.show()