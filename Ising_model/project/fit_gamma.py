import numpy as np
import pandas as pd

from find_gamma_nu import find_gamma_nu


def main():
    df = pd.read_csv("chi_at_betac.csv")
    popt, pcov = find_gamma_nu(df)
    print(popt[1], np.sqrt(pcov[1][1]))
    gamma_nu = popt[1]
    nu = 1
    d_nu = 0.003
    d_gamma_nu = np.sqrt(pcov[1][1])

    gamma = nu * gamma_nu
    d_gamma = gamma * np.abs((d_gamma_nu/gamma_nu)+(d_nu/nu))
    print(gamma, d_gamma)


if __name__ == "__main__":
    main()
    # result 1.74 pm 0.01