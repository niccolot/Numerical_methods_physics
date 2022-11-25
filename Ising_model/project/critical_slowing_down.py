import os
import re

from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from threading import Thread

import numpy as np
import pandas as pd


def correlation_magnetization(k_corr: int, m_corr):
    corr = 0
    N = len(m_corr) - k_corr
    m_avg = np.mean(m_corr)
    for index in range(N):
        corr += (m_corr[index] - m_avg) * (m_corr[index + k_corr] - m_avg)

    return corr / N


def save_correlation(beta_save, f_save, size_save):
    print("inizio")
    beta_float = float("0." + beta_save)
    data_func = np.loadtxt(f_save, delimiter="\t")
    df_func = pd.DataFrame(data_func)
    df_func.columns = ["index", "m", "m2", "m4", "e", "e2"]
    df_func["m"].apply(np.abs)
    k_list = [i_ + 1 for i_ in range(50)]
    corr_list = []
    for k in k_list:
        corr_list.append(correlation_magnetization(k, df_func["m"]))

    data_func = {
        "k": k_list,
        "correlation": corr_list,
        "beta": [beta_save for _ in range(len(k_list))]
    }
    beta_str = "{:.4f}".format(beta_float)
    pd.DataFrame(data_func).to_csv(f"correlation_data/beta={beta_str}_L={size_save}.csv", index=False)
    print(f"fine {beta_str} e {size_save}")


def main():
    sizes = [30, 40, 50]
    beta_allowed = ["440"]

    for i, size in enumerate(sizes):

        directory_main = f'dati_ising/L{size}'

        for filename_main in os.listdir(directory_main):
            f = os.path.join(directory_main, filename_main)
            if os.path.isfile(f):
                result = re.search('ising_(.*)L' + str(size), f)
                beta = result.group(1)
                if beta in beta_allowed:
                    Thread(target=save_correlation, args=(beta, f, size)).start()


def fit_func(size, a, z):
    return a * size ** z


def get_z():
    size_list = [30, 40, 50]
    beta = "0.4400"
    tau_list = []
    for size in size_list:
        df = pd.read_csv(f"correlation_data/beta={beta}_L={size}.csv")
        corr_list = df["correlation"]
        tau = corr_list.sum()
        tau_list.append(tau)

    return curve_fit(fit_func, size_list, tau_list)

if __name__ == "__main__":
    popt, pcov = get_z()
    size_list = [30, 40, 50]
    beta = "0.4400"
    tau_list = []
    for size in size_list:
        df = pd.read_csv(f"correlation_data/beta={beta}_L={size}.csv")
        corr_list = df["correlation"]
        tau = corr_list.sum()
        tau_list.append(tau)
    print(popt[1], np.sqrt(pcov[1][1]))
    size_array = np.linspace(29, 51, 1000)
    plt.errorbar(size_list, tau_list, fmt='.')
    plt.plot(size_array, fit_func(size_array, *popt))
    plt.title(r"Andamento di $\tau$ a $\beta=0.44$ per diversi $L$")
    plt.xlabel(r"$L$")
    plt.ylabel(r"$\tau$")
    plt.savefig("/Users/marcoparrinello/Desktop/ising_images/tau_L.png")
    # https://arxiv.org/pdf/1404.0209.pdf

