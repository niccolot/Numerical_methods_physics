import numpy as np
import os
import re
import pandas as pd
import matplotlib.pyplot as plt


def binder(mag_var):
    mag_var = np.array(mag_var)
    mag4 = np.mean(mag_var ** 4)
    mag2 = np.mean(mag_var ** 2)
    return 1 - mag4 / (3 * mag2 ** 2)


def variance_func(input_data, function, **kwargs):
    len_blocks = [10, 100, 500, 1000, 5000]
    variances_different_len = []
    for len_block in len_blocks:
        number_extraction = len(input_data) // len_block
        values_new_samples = []
        for j in range(100):
            indexes = np.random.choice(range(0, len(input_data) - len_block + 1),
                                       number_extraction,
                                       replace=True)
            new_indexes = []
            for index in indexes:
                elem = [item for item in range(index, index + len_block)]
                new_indexes += elem

            new_sample = input_data.iloc[new_indexes].reset_index(drop=True)
            values_new_samples.append(function(new_sample, **kwargs))

        values_new_samples = np.array(values_new_samples)
        _variance = values_new_samples.std()
        variances_different_len.append(_variance)
    return max(variances_different_len)


def variance_array(input_datas, function, **kwargs_):
    output_ = []
    print(len(input_datas))
    i = 0
    for input_data in input_datas:
        print(i + 1)
        i = i + 1
        variance = variance_func(input_data, function, **kwargs_)
        output_.append(variance)

    return output_

#
# if __name__ == "__main__":
#     sizes = ["L40"]
#     output = {}
#     for i, size in enumerate(sizes):
#
#         directory = 'dati_ising/' + size
#
#         betas = []
#         chis = []
#         Cs = []
#         dfs = []
#         Ms = []
# #
#         for filename in os.listdir(directory):
#             f = os.path.join(directory, filename)
#             if os.path.isfile(f):
#                 result = re.search('ising_(.*)' + size, f)
#                 data = np.loadtxt(f, delimiter="\t")
#                 df = pd.DataFrame(data)
#                 df.columns = ["index", "m", "m2", "m4", "e", "e2"]
#                 df["m"].apply(np.abs)
#                 df["index"] = df["index"].astype(int)
#                 dfs.append(df)
#                 betas.append(float("0." + result.group(1)))
#
#                 print(betas[-1])
#
#                 chis.append(chi_function(df, size))
#                 Cs.append(capacity_function(df, size))
#                 Ms.append(m_function(df))
#
#             data = {
#                 "beta": betas,
#                 "suscettività magnetica": chis,
#                 "capacità termica": Cs,
#                 "magnetizzazione": Ms,
#                 "datas": dfs
#             }
#
#             df = pd.DataFrame(data)
#
#             output[size] = df
#
#     data = output["L40"]["datas"][0]
#
#     print(data.head())
#     deltas = []
#     pisa = [10, 100, 1000, 2000, 3000, 4000, 6000, 8000, 10000]
#     for len_block in pisa:
#         deltas.append(variance(input_data=data,
#                                function=chi_function,
#                                len_block=len_block,
#                                size_lattice="L40"))
#
#     plt.scatter(pisa, deltas)
#     plt.show()
