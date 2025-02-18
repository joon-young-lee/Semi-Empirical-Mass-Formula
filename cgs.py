import numpy as np

c = 2.99792458 * 10 ** 10
c2 = c**2
e = 4.8032 * 10 ** (-10) # cgs 
# m_n = 939.57/c2 # MeV/c^2
m_n = 1.67492747 * 1.e-24
m_p = 1.67262192 * 1.e-24
#m_p = 938.27/c2 # MeV/c^2
m_e = 9.1093837 * 1.e-28
# m_e = 0.5110/c2 # MeV/c^2
hbar = 1.054572 * 1.e-27
# hbar = 6.582 * 1.e-10 # MeV * s
MeV = 1.60218e-6 # erg
fm = 1.e-13 # cm
def main():
    print(4.11 * 10 ** 11 * c2/ MeV * fm**3)
    return 0

if __name__ == '__main__':
    main()

