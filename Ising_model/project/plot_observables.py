import numpy as np
import os
import re
import pandas as pd
import matplotlib.pyplot as plt
from bootstrap import variance_array


def chi_function(dataframe, size_lattice):
    volume = size_lattice ** 2

    evm2 = dataframe["m2"].mean()
    evm = dataframe["m"].mean()

    chi = volume * (evm2 - evm ** 2)
    return chi


def m_function(dataframe):
    evm = dataframe["m"].mean()

    return evm


def capacity_function(dataframe, size_lattice):
    volume = size_lattice ** 2

    e = dataframe["e"].mean()
    e2 = dataframe["e2"].mean()

    capacity = volume * (e2 - e ** 2)

    return capacity


def binder(dataframe):
    evm4 = dataframe["m4"].mean()
    evm2 = dataframe["m2"].mean()

    return 1 - evm4 / (3 * evm2 ** 2)


if __name__ == "__main__":

    sizes = [20, 30, 40, 50]

    output = {}

    for i, size in enumerate(sizes):

        directory = f'dati_ising/L{size}'

        betas = []
        # chis = []
        # Cs = []
        binders = []
        dfs = []

        for filename in os.listdir(directory):
            f = os.path.join(directory, filename)
            if os.path.isfile(f):
                result = re.search('ising_(.*)L' + str(size), f)
                data = np.loadtxt(f, delimiter="\t")
                df = pd.DataFrame(data)
                df.columns = ["index", "m", "m2", "m4", "e", "e2"]
                df["m"].apply(np.abs)
                df["index"] = df["index"].astype(int)

                dfs.append(df)
                betas.append(float("0." + result.group(1)))
                # chis.append(chi_function(df, size))
                # Cs.append(capacity_function(df, size))
                binders.append(m_function(df))

        data = {
            "beta": betas,
            # "suscettività magnetica": chis,
            # "capacità termica": Cs,
            # "magnetizzazione": Ms,
            "binder": binders,
            "datas": dfs
        }

        df = pd.DataFrame(data)

        output[f"L{size}"] = df

    # ax = output["L20"].plot.scatter("beta", "magnetizzazione", marker="d", c="r")
    # output["L30"].plot.scatter("beta", "magnetizzazione", marker="v", c="g", ax=ax)
    # output["L50"].plot.scatter("beta", "magnetizzazione", marker="^", c="y", ax=ax)
    # output["L40"].plot.scatter("beta", "magnetizzazione", marker="s", c="b", ax=ax)

    markers = ["d", "v", "^", "s"]
    colours = ["r", "g", "y", "b"]
    for i, size in enumerate(sizes):
        print(size)
        x = output[f"L{size}"]["beta"]
        y = output[f"L{size}"]["binder"]
        dy = variance_array(input_datas=output[f"L{size}"]["datas"], function=binder)
        data = {
            "beta": x,
            "binder": y,
            "delta_binder": dy
        }

        df = pd.DataFrame(data)
        df.to_csv(f"binder_data_L{size}.csv", index=False)

        plt.errorbar(x, y, dy, fmt=markers[i], c=colours[i], markersize=0, capsize=4)

    plt.title("Cumulante di binder per diversi volumi")
    plt.xlabel(r"$\beta$")
    plt.ylabel(r"$B$")
    plt.legend(["L20", "L30", "L40", "L50"])
    plt.savefig("/Users/marcoparrinello/Desktop/ising_images/binder.png")
