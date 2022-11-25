import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d


def find_gamma_nu(data):
    # using the relation chi_betac = a L^gamma/nu
    chi_max = data["chi"]
    delta_chi_max = data["delta_chi"]
    sizes_f = [20, 30, 40, 50]
    fit_params = curve_fit(fit_exponents, sizes_f, chi_max,
                           p0=[0.12, 1.75], sigma=delta_chi_max)

    return fit_params


def fit_exponents(L, a, alpha):
    return a * L ** alpha


if __name__ == "__main__":
    sizes = [20, 30, 40, 50]
    beta_c = 0.4407
    new_data = []
    for size in sizes:
        df = pd.read_csv(f"suscettivita_data_L{size}.csv")
        f = interp1d(df["beta"], df["suscettivita"])
        df["sort_axis"] = np.abs(df["beta"] - beta_c)
        delta_chi = df.sort_values(by="sort_axis")["delta_suscettivita"].iloc[0]
        new_data.append(
            {
                "chi": float(f(beta_c)),
                "delta_chi": delta_chi
            }
        )

    new_df = pd.DataFrame(new_data)
    popt, pcov = find_gamma_nu(new_df)

    size_list = np.linspace(19, 51, 1000)

    plt.errorbar(sizes, new_df["chi"], new_df["delta_chi"], fmt=".", capsize=4)
    plt.plot(size_list, fit_exponents(size_list, *popt))
    plt.title(r"Grafico di $\chi$ nel punto critico in funzione di $L$")
    plt.xlabel(r"$L$")
    plt.ylabel(r"$\chi(\beta_c)$")
    # plt.savefig("/Users/marcoparrinello/Desktop/ising_images/gamma_nu.png")

    # print(popt[1], np.sqrt(pcov[1][1]))
    chi_square = np.sum((new_df["chi"] - fit_exponents(sizes, *popt)) ** 2 / new_df["delta_chi"] ** 2)
    print(chi_square / 2)
    # new_df.to_csv("chi_at_betac.csv", index=False)


