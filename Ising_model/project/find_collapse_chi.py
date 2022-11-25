import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import quad
from scipy.optimize import minimize, curve_fit
import matplotlib.pyplot as plt


def fit_function(x, a, b, c):
    return a + c * (x - b) ** 2


def transform_beta(L, beta, beta_c, inverse_nu):
    return (beta - beta_c) * L ** inverse_nu


def distance(size_list, beta_c, inverse_nu):
    final_distance = 0
    gamma_nu_ = 1.75
    for size_1 in size_list:
        for size_2 in size_list:
            if size_1 != size_2:
                df_1 = pd.read_csv(f"suscettivita_data_L{size_1}.csv")
                beta = transform_beta(size_1, df_1["beta"], beta_c, inverse_nu)
                f_1 = interp1d(beta, df_1["suscettivita"] / size_1 ** gamma_nu_, bounds_error=False, fill_value=0)

                df_2 = pd.read_csv(f"suscettivita_data_L{size_2}.csv")
                beta = transform_beta(size_2, df_2["beta"], beta_c, inverse_nu)
                f_2 = interp1d(beta, df_2["suscettivita"] / size_2 ** gamma_nu_, bounds_error=False, fill_value=0)

                final_distance += quad(lambda _x: (f_1(_x) - f_2(_x)) ** 2, -1.5, 1)[0]

    return final_distance


def func_to_minimize(params):
    size_list = [20, 30, 40, 50]
    return distance(size_list, params[0], params[1])


def find_collapsing_params(func):
    bounds = ((0.42, 0.46), (0.9, 1.1))  # bounds of params
    x0 = np.array([0.44, 1])  # initial guess of beta_c and inverse_nu

    output = minimize(fun=func, x0=x0, bounds=bounds, tol=10e-5)
    return output.x


def fit_beta_nu(ell, betaC, nu):
    return betaC - 0.412 * ell ** (- 1 / nu)


def main():
    sizes = [20, 30, 40, 50]
    markers = ["d", "v", "^", "s"]
    colours = ["r", "g", "y", "b"]
    gamma_nu = 1.75
    output = find_collapsing_params(func_to_minimize)
    beta_c = output[0]
    inverse_nu = output[1]

    df = pd.read_csv("beta_pc_vs_chi_max.csv")

    df["x_max"] = (df["beta_pc"] - beta_c) * df["size"] ** inverse_nu

    x_max = df["x_max"].mean()

    # we find x_max = -0.412

    popt, pcov = curve_fit(fit_beta_nu, df["size"], df["beta_pc"], p0=[0.44, 1], sigma=df["delta_beta_pc"])
    plt.errorbar(df["size"], df["beta_pc"], df["delta_beta_pc"],
                 fmt=markers[0], c=colours[0], markersize=0, capsize=4)

    x_list = np.linspace(19, 51, 1000)
    plt.plot(x_list, fit_beta_nu(x_list, *popt))
    plt.title(r"Fit di $\beta_{pc}$ in funzione di $L$")
    plt.xlabel(r"$L$")
    plt.ylabel(r"$\beta_{pc}$")

    x = np.array(sizes)
    y = df["beta_pc"].to_numpy()
    dy = df["delta_beta_pc"].to_numpy()

    chisq = np.sum ((y - fit_beta_nu(x, *popt))**2 / dy**2)

    ndof = 2

    print(chisq / ndof)

    # plt.savefig("/Users/marcoparrinello/Desktop/ising_images/beta_pc_vs_L.png")


if __name__ == "__main__":
    # main()
    sizes = [20, 30, 40, 50]
    markers = ["d", "v", "^", "s"]
    colours = ["r", "g", "y", "b"]
    gamma_nu = 1.725
    beta_c = 0.44
    nu = 1

    for i, size in enumerate(sizes):
        df = pd.read_csv(f"suscettivita_data_L{size}.csv")

        y = df["suscettivita"] / size ** gamma_nu
        dy = df["delta_suscettivita"] / size ** gamma_nu
        x = (df["beta"] - beta_c) * size

        plt.errorbar(x, y, dy, fmt=markers[i], c=colours[i], markersize=0, capsize=4)

    plt.title("Grafici della suscettibilità collassati")
    plt.xlabel(r"$x=(\beta-\beta_c) L^{1/\nu}$")
    plt.ylabel(r"$\chi/L^{\gamma/\nu}$")
    plt.legend(["L20", "L30", "L40", "L50"])
    plt.savefig("/Users/marcoparrinello/Desktop/ising_images/susceptiblity_collapsed.png")


