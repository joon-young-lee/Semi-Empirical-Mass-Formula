import numpy as np

c = 2.99792458 * 10 ** 10
c2 = c**2
e = 4.8032 * 10 ** (-10)
m_n = 939.57
m_p = 938.27
m_e = 9.1094 * 10 ** (-28)
MeV = 1.60218e-6 # erg
fm = 1.e-13 # cm
def main():
    print(4.11 * 10 ** 11 * c2/ MeV * fm**3)
    return 0

if __name__ == '__main__':
    main()

