import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


def fit_function(L, a, b):
    # https://arxiv.org/pdf/1706.02541.pdf
    return a + b * np.log(L)


def C_at_betac():
    sizes_f = [20, 30, 40, 50]
    beta_c = 0.4407
    c_at_betac = []
    for size_f in sizes_f:
        df_f = pd.read_csv(f"capacita_data_L{size_f}.csv")
        f = interp1d(df_f["beta"], df_f["capacita"])
        df_f["sort_axis"] = np.abs(df_f["beta"] - beta_c)
        delta_C = df_f.sort_values(by="sort_axis")["delta_capacita"].iloc[0]
        c_at_betac.append(
            {
                "C": float(f(beta_c)),
                "delta_C": delta_C
            }
        )

    return pd.DataFrame(c_at_betac)


df_c = C_at_betac()
sizes = [20, 30, 40, 50]

plt.errorbar(sizes, df_c["C"], df_c["delta_C"], fmt=".", capsize=4)

popt, pcov = curve_fit(fit_function, sizes, df_c["C"], p0=[0.5, 2.5])
size_list = np.linspace(19, 51, 1000)
plt.plot(size_list, fit_function(size_list, *popt))
plt.title(r"Grafico della capacità al punto critico in funzione di $L$")
plt.xlabel(r"$L$")
plt.ylabel(r"$C_c$")
print(popt, np.sqrt(pcov.diagonal()))
# il risultato è 0.54 pm 0.25, 2.59 pm 0.07

chi_square = np.sum((df_c["C"] - fit_function(sizes, *popt)) ** 2 / df_c["delta_C"] ** 2)
print(chi_square / 2)
# il risultato è 0.74

# plt.savefig("/Users/marcoparrinello/Desktop/ising_images/Cc_L.png")