import numpy as np
import os
import re
import pandas as pd
import matplotlib.pyplot as plt

sizes = ["L20", "L30", "L40", "L50"]

fig, axs = plt.subplots(len(sizes), 1, tight_layout=True)

for i, size in enumerate(sizes):

    directory = 'dati_ising/' + size
    betas = []
    for filename in os.listdir(directory):
        f = os.path.join(directory, filename)
        # checking if it is a file
        if os.path.isfile(f):
            result = re.search('ising_(.*)' + size, f)
            betas.append(float("0." + result.group(1)))
            if result.group(1) == "440":
                x = np.loadtxt(f, delimiter="\t")
                df = pd.DataFrame(x)
                df.columns = ["index", "m", "m2", "m4", "e", "e2"]
                df["index"] = df["index"].astype(int)
    axs[i].plot(df.index[:2000], df.m[:2000])
    axs[i].set_ylabel("L=" + size[1:])
    axs[i].set_ylim([0,1])

axs[3].set_xlabel("Passi della catena di Markov")
axs[0].title.set_text("Magnetizzazione per $\\beta=0.44$")
plt.savefig("/Users/marcoparrinello/Desktop/ising_images/magnetization_thermalization.png")
