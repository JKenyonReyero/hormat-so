import subprocess, os, re, shutil
from datetime import datetime
import numpy as np
import json
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from matplotlib.ticker import MaxNLocator
from scipy.interpolate import interp1d
from matplotlib.ticker import MultipleLocator

def read_int_xsec_n_l(file_path):
    # Currently unused but saving for later if needed
    data = np.load(file_path)

    xsec_data = {
        tuple(k.split("_")[0:3]): data[k]
        for k in data.files
    }
    return xsec_data

cwd = os.getcwd()  # directory of wherever the user ran the command

data_dirpath = "/poster_diff_xsec_20260317_144627/"
pot_type = ["LELO", "LENLO", "NL"]
all_data = dict()

for i in pot_type:
    data = np.loadtxt(cwd+data_dirpath+f"{i}/nodes_4/l_0/21.s2n_15,829")
    all_data[f"{i}_x"] = data[:,0]
    all_data[f"{i}_xsec"] = data[:,1]

plt.rcParams.update({
    "font.size": 30,
    "axes.labelsize": 38,
    "xtick.labelsize": 28,
    "ytick.labelsize": 28,
    "legend.fontsize": 20,
})

fig, ax= plt.subplots(1, 1, figsize=(15, 6))
plt.subplots_adjust(
    left=0.15,
    right=0.971,
    bottom=0.179,
    top=0.975,
    hspace=0.037,
    wspace=0.037
)
ax.plot(all_data["LELO_x"], all_data["LELO_xsec"], c="r", linestyle="dotted", label="LELO", linewidth=3)
ax.plot(all_data["LENLO_x"], all_data["LENLO_xsec"], c="r", linestyle="dashed", label="LENLO", linewidth=3)
ax.plot(all_data["NL_x"], all_data["NL_xsec"], c="r", linestyle="solid", label="NL", linewidth=3)
ax.set_yscale('log')
ax.set_xlim(xmin=0.0, xmax=180)

ax.xaxis.set_major_locator(MultipleLocator(30))

ax.set_ylabel("Cross section\n[mb/Sr]")
ax.set_xlabel("Scattering angle [deg]")
ax.legend(title=r"$E_p=$30MeV")
ax.grid(alpha=0.25)
plt.savefig(f"diff_xsec_Li30.svg", pad_inches=0.02, transparent=True)
plt.show()