import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

df = pd.read_csv("chi_at_betac.csv")

sizes = [20, 30, 40, 50]


def fit_function(L, a, alpha):
    return a * (L ** alpha)

output = curve_fit(fit_function, sizes, df["chi"], p0=[0.12, 1.75], sigma=df["delta_chi"])
popt = output[0]
pcov = output[1]

print(popt[1], np.sqrt(pcov[1][1]))

chi_square = np.sum((df["chi"] - fit_function(sizes, *popt)) ** 2 / df["delta_chi"] ** 2)
print(chi_square / 2)

plt.errorbar(sizes, df["chi"], df["delta_chi"], fmt='o', capsize=4)
size_list = np.linspace(19, 51, 1000)
plt.plot(size_list, fit_function(size_list, *popt))
plt.title(r"$\chi$ nel punto critico al variare di $L$")
plt.xlabel(r"L")
plt.ylabel(r"$\chi_{max}$")
# plt.savefig("/Users/marcoparrinello/Desktop/ising_images/chibetac_vs_L.png")

# si ottiene 1.743 pm 0.005 e chisq = 0.06
# plt.show()
