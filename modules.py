import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


# Define the binding energy model
def binding_energy(Z, N, av, as_, ac, ap, kv, ks, W, ak, a0, fp):
    A = Z + N  # Mass number
    I = (N - Z) / A  # Neutron-proton asymmetry

    # Volume term
    volume = av * (1 - kv * I**2) * A
    # Surface term
    surface = -as_ * (1 - ks * I**2) * A**(2/3)
    # Coulomb term
    coulomb = -ac * Z ** 2 / A ** (1/3)
    # Pairing energy term (vectorized)
    pairing = np.zeros_like(A)
    pairing[(Z % 2 == 0) & (N % 2 == 0)] = ap / A[(Z % 2 == 0) & (N % 2 == 0)]**0.5  # Even Z, even N
    pairing[(Z % 2 == 1) & (N % 2 == 1)] = -ap / A[(Z % 2 == 1) & (N % 2 == 1)]**0.5  # Odd Z, odd N
    # Wigner energy term
    wigner = -W * np.abs(N - Z) / A
    
    curvature_energy = - ak * A ** (1/3)

    higher_order = - a0 * A ** 0

    f_p = - fp * Z**2 / A

    return volume + surface + coulomb + pairing + wigner + curvature_energy + higher_order + f_p

def binding_energy_only_5(Z, N, av, as_, ac, kv, ks):
    A = Z + N  # Mass number
    I = (N - Z) / A  # Neutron-proton asymmetry

    # Volume term
    volume = av * (1 - kv * I**2) * A
    # Surface term
    surface = -as_ * (1 - ks * I**2) * A**(2/3)
    # Coulomb term
    coulomb = -ac * Z ** 2 / A ** (1/3)
    
    return volume + surface + coulomb

def standard_deviation(fit, data):
    sigma = 0
    n = len(fit)
    for i in range(n):
        sigma += np.abs(fit[i] - data[i]) ** 2
    sigma /= n
    sigma = np.sqrt(sigma)
    
    
    return sigma