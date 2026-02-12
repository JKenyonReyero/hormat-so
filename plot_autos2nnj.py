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
data_dirpath = "/dataTPMLi_20260211_154611/"  # temp
int_xsec_n_l = read_int_xsec_n_l(cwd+data_dirpath+"integrated_xsecs_n_l.npz")
s2n = read_s2n(cwd+data_dirpath+"s2n.txt")

data_dirpath = "/dataTPMPang_20260211_154353/"  # temp
int_xsec_n_l2 = read_int_xsec_n_l(cwd+data_dirpath+"integrated_xsecs_n_l.npz")
s2n2 = read_s2n(cwd+data_dirpath+"s2n.txt")

# print(int_xsec_n_l[("0","0","0")])
colours = ["#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999"]
# Plot
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8,6), sharex=True, sharey=True)
for i in range(0,3):
    ax1.plot(s2n, int_xsec_n_l[("0",f"{i}","0")], label=f"LELO n={i} l=0", c=colours[i], linestyle=":")

    ax1.plot(s2n, int_xsec_n_l[("1",f"{i}","0")], label=f"LENLO n={i} l=0", c=colours[i], linestyle="--")

    ax1.plot(s2n, int_xsec_n_l[("2",f"{i}","0")], label=f"NL n={i} l=0", c=colours[i], linestyle="-")

    ax2.plot(s2n2, int_xsec_n_l2[("0",f"{i}","0")], label=f"LELO n={i} l=0", c=colours[i], linestyle=":")

    ax2.plot(s2n2, int_xsec_n_l2[("1",f"{i}","0")], label=f"LENLO n={i} l=0", c=colours[i], linestyle="--")

    ax2.plot(s2n2, int_xsec_n_l2[("2",f"{i}","0")], label=f"NL n={i} l=0", c=colours[i], linestyle="-")

ax2.set_xlabel(r"$S_{2n}$ (MeV)", fontsize=16)
ax1.set_ylabel(r"$\sigma_R$ (mb/sr)", fontsize=16)
ax1.legend(fontsize=12)
ax2.legend(fontsize=12)
ax1.grid(alpha=0.2)
ax2.grid(alpha=0.2)

plt.show()


# Plot
fig, ax = plt.subplots(figsize=(8,6))
for i in range(0,6):
    ax.plot(s2n, int_xsec_n_l[("0",f"0",f"{i}")], label=f"LELO n=0 l={i}", c=colours[i], linestyle=":")

    ax.plot(s2n, int_xsec_n_l[("1",f"0",f"{i}")], label=f"LENLO n=0 l={i}", c=colours[i], linestyle="--")

    ax.plot(s2n, int_xsec_n_l[("2",f"0",f"{i}")], label=f"NL n=0 l={i}", c=colours[i], linestyle="-")

ax.set_xlabel(r"$S_{2n}$ (MeV)", fontsize=16)
ax.set_ylabel(r"$\sigma_R$ (mb/sr)", fontsize=16)
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

fig, ax = plt.subplots(figsize=(8,6))
for i in range(0,6):
    # ax.plot(s2n, int_xsec_n_l[("0",f"{i}","0")], label=f"NL n{i} l0")

    ax.plot(s2n, (int_xsec_n_l[("1","0",f"{i}")]-int_xsec_n_l[("2","0",f"{i}")])/int_xsec_n_l[("2","0",f"{i}")], label=f"LENLO/NL n=0 l={i}", linestyle="--")

    ax.plot(s2n, (int_xsec_n_l[("0","0",f"{i}")]-int_xsec_n_l[("2","0",f"{i}")])/int_xsec_n_l[("2","0",f"{i}")], label=f"LELO/NL n=0 l={i}", linestyle=":")

# ax.axhline(y=1,c="k",alpha=0.2,linewidth=2)
ax.set_xlabel(r"$S_{2n}$ (MeV)", fontsize=16)
ax.set_ylabel(r"$\sigma_R$ ratio (dimensionless)", fontsize=16)
ax.legend(fontsize=12)
ax.grid(alpha=0.2)

plt.show()


from matplotlib.lines import Line2D
plt.rcParams.update({
    "font.size": 30,
    "axes.labelsize": 38,
    "xtick.labelsize": 28,
    "ytick.labelsize": 28,
    "legend.fontsize": 28,
})
fig, ax = plt.subplots(figsize=(12,8))
plt.rcParams.update({
    "font.size": 30,
    "axes.labelsize": 38,
    "xtick.labelsize": 28,
    "ytick.labelsize": 28,
    "legend.fontsize": 28,
})

for i in range(6):
    ax.plot(s2n, int_xsec_n_l[("0","0",f"{i}")],
            color=colours[i], linestyle=":", linewidth=3.5)
    ax.plot(s2n, int_xsec_n_l[("1","0",f"{i}")],
            color=colours[i], linestyle="--", linewidth=3.5)
    ax.plot(s2n, int_xsec_n_l[("2","0",f"{i}")],
            color=colours[i], linestyle="-", linewidth=3.5)

# Legend for l (colors)
color_handles = [
    Line2D([0],[0], color=colours[i], lw=2, label=f"l={i}")
    for i in range(6)
]

# Legend for models (linestyles)
style_handles = [
    Line2D([0],[0], color="black", lw=3.5, linestyle=":", label="LELO"),
    Line2D([0],[0], color="black", lw=3.5, linestyle="--", label="LENLO"),
    Line2D([0],[0], color="black", lw=3.5, linestyle="-", label="NL"),
]

leg1 = ax.legend(handles=color_handles, loc="upper left", frameon=False)
leg2 = ax.legend(handles=style_handles, title="Locality", loc="upper right", frameon=False)

ax.add_artist(leg1)

ax.set_xlabel(r"$S_{2n}$ (MeV)")
ax.set_ylabel(r"$\sigma_R$ (mb/sr)")
ax.grid(alpha=0.2)

plt.show()