import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def M_at_betac():
    sizes = [20, 30, 40, 50]
    beta_c = 0.4407
    m_at_betac = []
    for size in sizes:
        df = pd.read_csv(f"magnetizzazione_data_L{size}.csv")
        f = interp1d(df["beta"], df["magnetizzazione"])
        df["sort_axis"] = np.abs(df["beta"] - beta_c)
        delta_M = df.sort_values(by="sort_axis")["delta_magnetizzazione"].iloc[0]
        m_at_betac.append(
            {
                "M": float(f(beta_c)),
                "delta_M": delta_M
            }
        )

    return pd.DataFrame(m_at_betac)


def fit_beta_nu_function(L, a, alpha):
    return a * L ** alpha


def fit_beta_nu():
    df_m = M_at_betac()
    sizes = [20, 30, 40, 50]
    popt, pcov = curve_fit(fit_beta_nu_function, sizes, df_m["M"], p0=[1, -0.125], sigma=df_m["delta_M"])
    plt.errorbar(sizes, df_m["M"], df_m["delta_M"], fmt=".", capsize=4)
    size_array = np.linspace(19, 51, 1000)
    plt.plot(size_array, fit_beta_nu_function(size_array, *popt))
    plt.title(r"Fit della magnetizzazione al punto critico in funzione di $L$")
    plt.xlabel(r"$L$")
    plt.ylabel(r"$\langle |M| \rangle_C$")
    plt.savefig("/Users/marcoparrinello/Desktop/ising_images/magnetizzazione_critica.png")

    chi_square = np.sum((df_m["M"] - fit_beta_nu_function(sizes, *popt)) ** 2 / df_m["delta_M"] ** 2)
    print(chi_square / 2)

    print(popt[1], np.sqrt(pcov[1][1]))
# il risultato è beta/nu = 0.1232 pm 0.0009, chisq = 0.11
fit_beta_nu()