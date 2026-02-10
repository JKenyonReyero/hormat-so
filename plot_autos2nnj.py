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
data_dirpath = "/datatpm_20260210_170917/"  # temp
int_xsec_n_l = read_int_xsec_n_l(cwd+data_dirpath+"integrated_xsecs_n_l.npz")
s2n = read_s2n(cwd+data_dirpath+"s2n.txt")

# print(int_xsec_n_l[("0","0","0")])

# Plot
fig, ax = plt.subplots(figsize=(8,6))
for i in range(0,3):
    ax.plot(s2n, int_xsec_n_l[("0",f"{i}","0")], label=f"LELO n={i} l=0", linestyle=":")

    ax.plot(s2n, int_xsec_n_l[("1",f"{i}","0")], label=f"LENLO n={i} l=0", linestyle="--")

    ax.plot(s2n, int_xsec_n_l[("2",f"{i}","0")], label=f"NL n={i} l=0", linestyle="-")

ax.set_xlabel(r"$S_{2n}$ (MeV)", fontsize=16)
ax.set_ylabel(r"Integrated cross section $\sigma_R$ (mb/sr)", fontsize=16)
ax.legend(fontsize=12)
ax.grid(alpha=0.2)

plt.show()


# Plot
fig, ax = plt.subplots(figsize=(8,6))
colours = ["r", "b", "k", "g", "pink", "y"]
for i in range(0,6):
    ax.plot(s2n, int_xsec_n_l[("0",f"0",f"{i}")], label=f"LELO n=0 l={i}", c=colours[i], linestyle=":")

    ax.plot(s2n, int_xsec_n_l[("1",f"0",f"{i}")], label=f"LENLO n=0 l={i}", c=colours[i], linestyle="--")

    ax.plot(s2n, int_xsec_n_l[("2",f"0",f"{i}")], label=f"NL n=0 l={i}", c=colours[i], linestyle="-")

ax.set_xlabel(r"$S_{2n}$ (MeV)", fontsize=16)
ax.set_ylabel(r"Integrated cross section $\sigma_R$ (mb/sr)", fontsize=16)
ax.legend(fontsize=12)
ax.grid(alpha=0.2)

plt.show()


fig, ax = plt.subplots(figsize=(8,6))
for i in range(0,3):
    # ax.plot(s2n, int_xsec_n_l[("0",f"{i}","0")], label=f"NL n{i} l0")

    ax.plot(s2n, int_xsec_n_l[("1",f"{i}","0")]/int_xsec_n_l[("2",f"{i}","0")], label=f"LENLO/NL n={i} l=0", linestyle="--")

    ax.plot(s2n, int_xsec_n_l[("0",f"{i}","0")]/int_xsec_n_l[("2",f"{i}","0")], label=f"LELO/NL n={i} l=0", linestyle=":")

ax.axhline(y=1,c="k",alpha=0.2,linewidth=2)
ax.set_xlabel(r"$S_{2n}$ (MeV)", fontsize=16)
ax.set_ylabel(r"$\sigma_R$ ratio (dimensionless)", fontsize=16)
ax.legend(fontsize=12)
ax.grid(alpha=0.2)

plt.show()


fig, ax = plt.subplots(figsize=(8,6))
for i in range(0,6):
    # ax.plot(s2n, int_xsec_n_l[("0",f"{i}","0")], label=f"NL n{i} l0")

    ax.plot(s2n, int_xsec_n_l[("1","0",f"{i}")]/int_xsec_n_l[("2","0",f"{i}")], label=f"LENLO/NL n=0 l={i}", linestyle="--")

    ax.plot(s2n, int_xsec_n_l[("0","0",f"{i}")]/int_xsec_n_l[("2","0",f"{i}")], label=f"LELO/NL n=0 l={i}", linestyle=":")

ax.axhline(y=1,c="k",alpha=0.2,linewidth=2)
ax.set_xlabel(r"$S_{2n}$ (MeV)", fontsize=16)
ax.set_ylabel(r"$\sigma_R$ ratio (dimensionless)", fontsize=16)
ax.legend(fontsize=12)
ax.grid(alpha=0.2)

plt.show()

