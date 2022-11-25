import os
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from bootstrap import variance_array


def binder(dataframe):
    evm4 = dataframe["m4"].mean()
    evm2 = dataframe["m2"].mean()

    return 1 - evm4 / (3 * evm2 ** 2)


sizes = ["L20", "L30", "L40", "L50"]

markers = ["d", "v", "^", "s"]
colours = ["r", "g", "y", "b"]


for i, size in enumerate(sizes):

    directory = 'dati_ising/' + size

    betas = []
    binders = []
    delta_binders = []

    for filename in os.listdir(directory):
        f = os.path.join(directory, filename)
        if os.path.isfile(f):
            result = re.search('ising_(.*)' + size, f)
            data = np.loadtxt(f, delimiter="\t")
            df = pd.DataFrame(data)
            df.columns = ["index", "m", "m2", "m4", "e", "e2"]
            betas.append(float("0." + result.group(1)))
            binders.append(binder(df))
            delta_binders.append(variance_array(input_datas=[df], function=binder)[0])

    plt.errorbar(betas, binders, delta_binders, fmt=markers[i], c=colours[i], markersize=0, capsize=4)
plt.show()
