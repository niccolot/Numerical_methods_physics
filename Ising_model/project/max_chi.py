import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np

sizes = ["L20", "L30", "L40", "L50"]
widths = [0.008, 0.0080, 0.0070, 0.0060]
dfs = []
max_values = []
points = []

beta_pc = []


def max_location(dataframe):
    max_chi_f = dataframe["suscettivita"].max()
    max_value_f = dataframe.iloc[dataframe["suscettivita"].idxmax()]["beta"]
    return max_value_f, max_chi_f


def fit_function(x, chi_max, b, beta_pc):
    return chi_max + b * (x - beta_pc) ** 2

# for size in sizes:
#     df = pd.read_csv(f"suscettivita_data_{size}.csv")
#
#     df = df[df.beta.between(0.40, 0.44)]
#
#     x = df["beta"]
#     y = df["suscettivita"]
#     dy = df["delta_suscettivita"]
#
#     plt.errorbar(x, y, dy, fmt=".", capsize=4)
#
#     plt.show()

for size, width in zip(sizes, widths):
    df = pd.read_csv(f"suscettivita_data_{size}.csv")
    max_value, max_chi = max_location(df)
    point = df[np.abs(df["beta"] - max_value) < width]
    x_fit = point["beta"]
    y_fit = point["suscettivita"]
    dy_fit = point["delta_suscettivita"]

    plt.errorbar(x_fit, y_fit, dy_fit, fmt="o")
    popt, pcov = curve_fit(f=fit_function,
                           xdata=x_fit,
                           ydata=y_fit,
                           sigma=dy_fit,
                           p0=[max_chi, 5 * 10 ** 4, max_value])
    x_test = np.linspace(min(x_fit), max(x_fit), 1000)
    y_test = fit_function(x_test, *popt)
    plt.plot(x_test, y_test)
    plt.show()

    point = df[np.abs(df["beta"] - popt[-1]) < width]
    x_fit = point["beta"]
    y_fit = point["suscettivita"]
    dy_fit = point["delta_suscettivita"]

    plt.errorbar(x_fit, y_fit, dy_fit, fmt="o")
    popt, pcov = curve_fit(f=fit_function,
                           xdata=x_fit,
                           ydata=y_fit,
                           sigma=dy_fit,
                           p0=[max_chi, 5 * 10 ** 4, max_value])
    x_test = np.linspace(min(x_fit), max(x_fit), 1000)
    y_test = fit_function(x_test, *popt)
    plt.plot(x_test, y_test)
    plt.show()
    new_elem = {
        "beta_pc": popt[-1],
        "delta_beta_pc": np.sqrt(pcov[-1][-1]),
        "chi_max": popt[0],
        "delta_chi_max": np.sqrt(pcov[0][0]),
        "size": int(size[1:])
    }
    beta_pc.append(new_elem)


df = pd.DataFrame(beta_pc)
print(df)
df.to_csv("beta_pc_vs_chi_max.csv", index=False)
