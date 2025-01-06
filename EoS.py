import numpy as np
import cgs
import modules
import matplotlib.pyplot as plt
import multiprocessing as mp
from fitting import a_v, a_s, a_c, a_p, k_v, k_s, W, a_k, a_0, f_p, b_1, b_2

def E_tot_pressure(Z, N, n_b):
    A = N + Z
    n_N = n_b / A
    n_e = n_b * Z / A
    # Energy of Nuclei
    # Binding Energy Already Negative
    E_N = cgs.m_p * Z * cgs.c2 + cgs.m_n * N * cgs.c2 - \
        cgs.MeV * modules.binding_energy(Z, N, a_v, a_s, a_c, a_p, k_v, k_s, W, a_k, a_0, f_p, b_1, b_2)
    
    
    # Lattice Energy
    a = (2/n_N) ** (1/3)
    E_L = -1.819620 * Z ** 2 * cgs.e ** 2 / a
    
    
    # Electron Energy
    k_e = (3 * np.pi ** 2 * n_e) ** (1/3)
    t = cgs.hbar * k_e / (cgs.m_e * cgs.c)
    # E_e = cgs.m_e ** 4 * cgs.c ** 5 / \
    # (8 * np.pi ** 2 * cgs.hbar ** 3) * ((2 * t ** 3 + t) * np.sqrt(t ** 2 + 1) - np.log(t + np.sqrt(t**2+1)))
    E_e = 3/4 * n_e * cgs.m_e * cgs.c2 * t

    E_tot = n_N * (E_N + E_L) + E_e
    pressure = 1 / 3 * E_L * n_N + n_e * cgs.c * cgs.hbar * (3 * np.pi ** 2 * n_e) ** (1/3) - E_e
    
    return E_tot, pressure


def minimum_tot(lower_Z, upper_Z, lower_N, upper_N, n_b): # Bounds should be int
    E, P = E_tot_pressure(np.array([lower_Z]), np.array([lower_N]), n_b)
    for i in range(lower_Z, upper_Z + 1):
        for j in range(lower_N, upper_N + 1):
            if E_tot_pressure(np.array([i]), np.array([j]), n_b)[0] < E:
                E, P = E_tot_pressure(np.array([i]), np.array([j]), n_b)
                Z = i
                N = j
            else:
                continue
    A = Z + N
    return E/cgs.c2, P, Z, A

def find_minimum_for_nb(n_b, lower_Z, upper_Z, lower_N, upper_N):
    E_min, P_min = float('inf'), None
    Z_min, N_min = None, None
    for Z in range(lower_Z, upper_Z + 1):
        for N in range(lower_N, upper_N + 1):
            E, P = E_tot_pressure(np.array([Z]), np.array([N]), n_b)
            if E < E_min:
                E_min, P_min = E, P
                Z_min, N_min = Z, N
    return E_min, P_min, Z_min, N_min


def process_nb(args):
    n_b, lower_Z, upper_Z, lower_N, upper_N = args
    return find_minimum_for_nb(n_b, lower_Z, upper_Z, lower_N, upper_N)



if __name__ == "__main__":
    data = np.loadtxt('BPS.txt')
    data_MassDensity = data[:, 0]
    data_Pressure = data[:, 1]
    data_nb = data[:, 2]

    lower_Z, upper_Z = 24, 80
    lower_N, upper_N = 30, 100
    E_41, P = E_tot_pressure(np.array([41]), np.array([82]), 1.105E+35)
    E_42, P = E_tot_pressure(np.array([42]), np.array([82]), 1.105E+35)
    print(E_41)
    print(E_42)

    # Prepare arguments for multiprocessing
    args = [(n_b, lower_Z, upper_Z, lower_N, upper_N) for n_b in data_nb]

    with mp.Pool(mp.cpu_count()) as pool:
        results = pool.map(process_nb, args)

    # Process results
    MassDensity = np.zeros_like(data_nb)
    Pressure = np.zeros_like(data_nb)
    Z_array = np.zeros_like(data_nb, dtype=int)
    A_array = np.zeros_like(data_nb, dtype=int)

    for k, (E, P, Z, N) in enumerate(results):
        MassDensity[k] = E / cgs.c2
        Pressure[k] = P
        Z_array[k] = Z
        A_array[k] = Z + N

        if MassDensity[k] <= 4.3 * 1.e11:
            print('---------------------------------------')
            print(f'Baryon Density: {data_nb[k]}')
            print(f'Mass Density: {MassDensity[k]:e}')
            print(f'Pressure: {Pressure[k]:e}')
            print(f'Proton Number: {Z_array[k]}')
            print(f'Total Mass Number: {A_array[k]}')
            print('---------------------------------------')
            l = k
        else:
            break

    # Save results
    calculations = np.column_stack((MassDensity[:k], Pressure[:k], data_nb[:k], Z_array[:k], A_array[:k]))
    np.savetxt('EoS.txt', calculations, fmt=['%.3E', '%.3E', '%.3E', '%d', '%d'], delimiter=' ',
               header='# rho (g/cm^3)  P (dynes/cm^2)  n_B (cm^-3) Z A', comments='')