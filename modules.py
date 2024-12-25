import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# Define magic numbers
MAGIC_NUMBERS = [2, 8, 20, 28, 50, 82, 126, 184]

def count_valence_nucleons(N, Z):
    """Calculate the valence nucleons for neutrons and protons."""
    def valence(n, magic):
        distances = np.array([min(np.abs(n - magic[i]), np.abs(n - magic[i + 1]))
                              if i + 1 < len(magic) else np.inf for i in range(len(magic))])
        return distances.min(axis=0)

    nv = np.array([valence(n, MAGIC_NUMBERS) for n in np.atleast_1d(N)])
    zv = np.array([valence(z, MAGIC_NUMBERS) for z in np.atleast_1d(Z)])
    return nv, zv

# Define the binding energy model
def binding_energy(Z, N, av, as_, ac, ap, kv, ks, W, ak, a0, fp, b1, b2):
    A = Z + N  # Mass number
    I = (N - Z) / A  # Neutron-proton asymmetry

    # Volume term
    volume = av * (1 - kv * I**2) * A
    # Surface term
    surface = -as_ * (1 - ks * I**2) * A**(2 / 3)
    # Coulomb term
    coulomb = -ac * Z**2 / A**(1 / 3)
    # Pairing term
    pairing = np.zeros_like(A)
    pairing[(Z % 2 == 0) & (N % 2 == 0)] = ap / A[(Z % 2 == 0) & (N % 2 == 0)]**0.5  # Even Z, even N
    pairing[(Z % 2 == 1) & (N % 2 == 1)] = -ap / A[(Z % 2 == 1) & (N % 2 == 1)]**0.5  # Odd Z, odd N
    # Wigner term
    wigner = -W * np.abs(N - Z) / A
    # Curvature and higher-order terms
    curvature = -ak * A**(1 / 3)
    higher_order = -a0 * A**0
    f_p = -fp * Z**2 / A

    # Shell corrections
    nv, zv = count_valence_nucleons(N, Z)
    shell = b1 * (nv + zv) + b2 * (nv + zv)**2

    return volume + surface + coulomb + pairing + wigner + curvature + higher_order + f_p + shell

def standard_deviation(fit, data):
    sigma = 0
    n = len(fit)
    for i in range(n):
        sigma += np.abs(fit[i] - data[i]) ** 2
    sigma /= n
    sigma = np.sqrt(sigma)
    
    
    return sigma