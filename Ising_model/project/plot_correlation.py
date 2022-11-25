import pandas as pd
from matplotlib import pyplot as plt

beta_list = ["0.4100", "0.4200", "0.4400"]
size = "50"

for i, beta in enumerate(beta_list):
    df = pd.read_csv(f"correlation_data/beta={beta}_L={size}.csv")
    x = df["k"]
    y = df["correlation"]
    plt.plot(x, y)
    plt.legend([r"$\beta$=" + beta for beta in beta_list])
    plt.title(fr"Critical slowing down per $L={size}$ a diversi $\beta$")
    plt.xlabel(r"k")
    plt.ylabel(r"C(k)")
# plt.savefig(f"/Users/marcoparrinello/Desktop/ising_images/correlation_L{size}.png")
plt.show()