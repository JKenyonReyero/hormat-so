import subprocess, os, re, shutil
from datetime import datetime
import numpy as np
import json
import matplotlib.pyplot as plt

def read_int_xsec_n_l(file_path):
    # Currently unused but saving for later if needed
    data = np.load(file_path)

    xsec_data = {
        tuple(k.split("_")[0:3]): data[k]
        for k in data.files
    }
    return xsec_data

def read_s2n(file_path):
    with open(file_path, "r") as f:
        arr = np.array(json.load(f))
    return arr

#--------------------------------------------------------------
# Set up directories
#--------------------------------------------------------------
script_dir = os.path.dirname(os.path.abspath(__file__))

cwd = os.getcwd()  # directory of wherever the user ran the command

# data_dirpath = str(input("directory path from the current working directory include the / but no .\n"))
# if data_dirpath[-1] != "/":
#     data_dirpath = data_dirpath+"/"
data_dirpath = "/dataTPM_20260205_193642/"  # temp
int_xsec_n_l = read_int_xsec_n_l(cwd+data_dirpath+"integrated_xsecs_n_l.npz")
s2n = read_s2n(cwd+data_dirpath+"s2n.txt")

# print(int_xsec_n_l[("0","0","0")])

# Plot
fig, ax = plt.subplots()
for i in range(0,3):
    ax.plot(s2n, int_xsec_n_l[("0",f"{i}","0")], marker="o", label=f"NL n{i} l0")

    ax.plot(s2n, int_xsec_n_l[("1",f"{i}","0")], marker="o", label=f"LENLO n{i} l0")

    ax.plot(s2n, int_xsec_n_l[("2",f"{i}","0")], marker="o", label=f"LELO n{i} l0")

ax.set_xlabel(r"$S_{2n}$ (MeV)")
ax.set_ylabel(r"Integrated cross section $\sigma_R$ (mb)")
ax.legend()
ax.grid(alpha=0.2)

plt.show()


# Plot
fig, ax = plt.subplots()
for i in range(0,6):
    ax.plot(s2n, int_xsec_n_l[("0",f"0",f"{i}")], marker="o", label=f"NL n0 l{i}")

    ax.plot(s2n, int_xsec_n_l[("1",f"0",f"{i}")], marker="o", label=f"LENLO n0 l{i}")

    ax.plot(s2n, int_xsec_n_l[("2",f"0",f"{i}")], marker="o", label=f"LELO n0 l{i}")

ax.set_xlabel(r"$S_{2n}$ (MeV)")
ax.set_ylabel(r"Integrated cross section $\sigma_R$ (mb)")
ax.legend()
ax.grid(alpha=0.2)

plt.show()