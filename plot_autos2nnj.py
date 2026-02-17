import subprocess, os, re, shutil
from datetime import datetime
import numpy as np
import json
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


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

def read_pnl(file_path):
    with open(file_path, "r") as f:
        arr = json.load(f)
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
pnl = read_pnl(cwd+data_dirpath+"pnl.txt")

data_dirpath = "/dataTPMPang_20260211_154353/"  # temp
int_xsec_n_l2 = read_int_xsec_n_l(cwd+data_dirpath+"integrated_xsecs_n_l.npz")
s2n2 = read_s2n(cwd+data_dirpath+"s2n.txt")
pnl2 = read_pnl(cwd+data_dirpath+"pnl.txt")


# print(int_xsec_n_l[("0","0","0")])
colours = ["#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999"]


# Legend for models (linestyles)
style_handles = [
    Line2D([0],[0], color="black", lw=3.5, linestyle=":", label="LELO"),
    Line2D([0],[0], color="black", lw=3.5, linestyle="--", label="LENLO"),
    Line2D([0],[0], color="black", lw=3.5, linestyle="-", label="NL"),
]
color_handles_n = [
    Line2D([0],[0], color=colours[i], lw=2, label=f"l={i}")
    for i in range(6)
]

x_offset = 0.02 * (max(s2n) - min(s2n))
x_offset2 = 0.02 * (max(s2n2) - min(s2n2))


# Plot
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8,6), sharex=True, sharey=True)
for i in range(0,3):
    ax1.plot(s2n, int_xsec_n_l[("0",f"{i}","0")], c=colours[i], linestyle=":")

    ax1.plot(s2n, int_xsec_n_l[("1",f"{i}","0")], c=colours[i], linestyle="--")

    NL_line, = ax1.plot(s2n, int_xsec_n_l[("2",f"{i}","0")], c=colours[i], linestyle="-")

    # --- Direct label using NL curve ---
    y_end = NL_line.get_ydata()[-1]

    ax1.text(
        s2n[-1] + x_offset,
        y_end,
        rf"$n={i}$",
        color=colours[i],
        fontsize=28,
        va="center",
    )

    ax2.plot(s2n2, int_xsec_n_l2[("0",f"{i}","0")], c=colours[i], linestyle=":")

    ax2.plot(s2n2, int_xsec_n_l2[("1",f"{i}","0")], c=colours[i], linestyle="--")

    NL_line, = ax2.plot(s2n2, int_xsec_n_l2[("2",f"{i}","0")], c=colours[i], linestyle="-")

    # --- Direct label using NL curve ---
    y_end = NL_line.get_ydata()[-1]

    ax2.text(
        s2n2[-1] + x_offset2,
        y_end,
        rf"$n={i}$",
        color=colours[i],
        fontsize=28,
        va="center",
    )

    

leg = ax1.legend(handles=style_handles, title="Locality", loc="upper right", frameon=False)
ax1.add_artist(leg)

ax2.set_xlabel(r"$S_{2n}$ (MeV)", fontsize=16)
ax1.set_ylabel(r"$\sigma_R$ (mb/sr)", fontsize=16)
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


plt.rcParams.update({
    "font.size": 30,
    "axes.labelsize": 38,
    "xtick.labelsize": 28,
    "ytick.labelsize": 28,
    "legend.fontsize": 18,
})
fig, ax = plt.subplots(figsize=(12,8))

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

leg1 = ax.legend(handles=color_handles, loc="upper left", frameon=False)
leg2 = ax.legend(handles=style_handles, title="Locality", loc="upper right", frameon=False)

ax.add_artist(leg1)

ax.set_xlabel(r"$S_{2n}$ (MeV)")
ax.set_ylabel(r"$\sigma_R$ (mb/sr)")
ax.grid(alpha=0.2)

plt.show()



x_offset = 0.02 * (max(s2n) - min(s2n))

k = 3/(93)  # m_triton / (mtriton+90Zr)
try:
    triton_E = k*(float(pnl["elab"])+s2n-8.48)  # k(E_p + S2n - triton binding energy)
except:
    triton_E = k*(28.53+s2n-8.48)  # Assume eE_p is 28.53MeV

print(triton_E)
fig, (ax1, ax2) = plt.subplots(
    1, 2,
    figsize=(8, 12),
    sharex=True,
    sharey=True
)

for i in range(6):

    # LELO
    line1, = ax1.plot(s2n, int_xsec_n_l[("0","0",f"{i}")], color=colours[i], linestyle=":", linewidth=3)

    # LENLO
    line2, = ax1.plot(s2n, int_xsec_n_l[("1","0",f"{i}")], color=colours[i], linestyle="--", linewidth=3)

    # NL (reference → slightly thicker)
    line3, = ax1.plot(s2n, int_xsec_n_l[("2","0",f"{i}")], color=colours[i], linestyle="-", linewidth=3.5)

    # # --- Direct label using NL curve ---
    # y_end = line3.get_ydata()[-1]
    # vertical_offset = 0.05 * 2

    # ax1.text(
    #     s2n[-1] + x_offset,
    #     y_end+ (i + 2.5) * vertical_offset,
    #     rf"$\ell={i}$",
    #     color=colours[i],
    #     fontsize=28,
    #     va="center",
    # )

        # LELO
    line1, = ax2.plot(s2n2, int_xsec_n_l2[("0","0",f"{i}")], color=colours[i], linestyle=":", linewidth=3)

    # LENLO
    line2, = ax2.plot(s2n2, int_xsec_n_l2[("1","0",f"{i}")], color=colours[i], linestyle="--", linewidth=3)

    # NL (reference → slightly thicker)
    line3, = ax2.plot(s2n2, int_xsec_n_l2[("2","0",f"{i}")], color=colours[i], linestyle="-", linewidth=3.5)

    # # --- Direct label using NL curve ---
    # y_end = line3.get_ydata()[-1]

    # ax2.text(
    #     s2n2[-1] + x_offset2,
    #     y_end+ (i - 2.5) * vertical_offset,
    #     rf"$\ell={i}$",
    #     color=colours[i],
    #     fontsize=28,
    #     va="center",
    # )

# Legend on each axis
ax1.legend(
    handles=style_handles,
    title="Locality",
    loc="upper right",
    frameon=False
)
# Legend on each axis
ax2.legend(
    handles=color_handles,
    loc="upper right",
    frameon=False
)

ax1.text(
    0.5, 0.9,                      # x=middle, y=slightly above top
    "Li",
    transform=ax1.transAxes,        # <- key line
    ha="center",
    va="bottom",
    fontsize=20,
)

ax2.text(
    0.5, 0.9,                      # x=middle, y=slightly above top
    "Pang",
    transform=ax2.transAxes,        # <- key line
    ha="center",
    va="bottom",
    fontsize=20,
)

# Labels only on outer axes (cleaner look)
ax2.set_xlabel(r"$S_{2n}$ (MeV)")
ax1.set_ylabel(r"$\sigma_R$ (mb/sr)")
ax2.set_ylabel(r"$\sigma_R$ (mb/sr)")

ax1.grid(alpha=0.2)
ax2.grid(alpha=0.2)

plt.subplots_adjust(
    left=0.125,
    right=0.95,
    bottom=0.125,
    top=0.9,
    hspace=0.0,
)

import numpy as np

# Add a twin axis on top
ax_top = ax1.twiny()

# Get main axis limits
x_main_min, x_main_max = ax1.get_xlim()

# Original triton_E range
triton_min, triton_max = min(triton_E), max(triton_E)

# Linear map triton_E to span the full main x-axis
triton_mapped = (triton_E - triton_min) / (triton_max - triton_min) * (x_main_max - x_main_min) + x_main_min

# Stretch twin axis
ax_top.set_xlim(ax1.get_xlim())

# Create ticks every 0.25 in triton_E space
tick_values = np.arange(triton_min, triton_max + 0.25, 0.1)

# Map tick_values to the main axis
tick_positions = (tick_values - triton_min) / (triton_max - triton_min) * (x_main_max - x_main_min) + x_main_min

ax_top.set_xticks(tick_positions)
ax_top.set_xticklabels([f"{t:.2f}" for t in tick_values])

# Example vertical line (scaled)
ax_top.axvline(
    x=x_main_min + (7.8 - triton_min) / (triton_max - triton_min) * (x_main_max - x_main_min),
    color='r',
    linestyle='--'
)
plt.show()





fig, (ax1, ax2) = plt.subplots(
    1, 2,
    figsize=(8, 12),
    sharex=True,
    sharey=True
)

for i in range(6):

    # LELO
    line1, = ax1.plot(s2n, int_xsec_n_l[("0","1",f"{i}")], color=colours[i], linestyle=":", linewidth=3)

    # LENLO
    line2, = ax1.plot(s2n, int_xsec_n_l[("1","1",f"{i}")], color=colours[i], linestyle="--", linewidth=3)

    # NL (reference → slightly thicker)
    line3, = ax1.plot(s2n, int_xsec_n_l[("2","1",f"{i}")], color=colours[i], linestyle="-", linewidth=3.5)

    # # --- Direct label using NL curve ---
    # y_end = line3.get_ydata()[-1]
    # vertical_offset = 0.05 * 2

    # ax1.text(
    #     s2n[-1] + x_offset,
    #     y_end+ (i + 2.5) * vertical_offset,
    #     rf"$\ell={i}$",
    #     color=colours[i],
    #     fontsize=28,
    #     va="center",
    # )

        # LELO
    line1, = ax2.plot(s2n2, int_xsec_n_l2[("0","1",f"{i}")], color=colours[i], linestyle=":", linewidth=3)

    # LENLO
    line2, = ax2.plot(s2n2, int_xsec_n_l2[("1","1",f"{i}")], color=colours[i], linestyle="--", linewidth=3)

    # NL (reference → slightly thicker)
    line3, = ax2.plot(s2n2, int_xsec_n_l2[("2","1",f"{i}")], color=colours[i], linestyle="-", linewidth=3.5)

    # # --- Direct label using NL curve ---
    # y_end = line3.get_ydata()[-1]

    # ax2.text(
    #     s2n2[-1] + x_offset2,
    #     y_end+ (i - 2.5) * vertical_offset,
    #     rf"$\ell={i}$",
    #     color=colours[i],
    #     fontsize=28,
    #     va="center",
    # )

# Legend on each axis
ax1.legend(
    handles=style_handles,
    title="Locality",
    loc="upper right",
    frameon=False
)
# Legend on each axis
ax2.legend(
    handles=color_handles,
    loc="upper right",
    frameon=False
)

ax1.text(
    0.5, 0.9,                      # x=middle, y=slightly above top
    "Li",
    transform=ax1.transAxes,        # <- key line
    ha="center",
    va="bottom",
    fontsize=20,
)

ax2.text(
    0.5, 0.9,                      # x=middle, y=slightly above top
    "Pang",
    transform=ax2.transAxes,        # <- key line
    ha="center",
    va="bottom",
    fontsize=20,
)

# Labels only on outer axes (cleaner look)
ax2.set_xlabel(r"$S_{2n}$ (MeV)")
ax1.set_ylabel(r"$\sigma_R$ (mb/sr)")
ax2.set_ylabel(r"$\sigma_R$ (mb/sr)")

ax1.grid(alpha=0.2)
ax2.grid(alpha=0.2)

plt.subplots_adjust(
    left=0.125,
    right=0.95,
    bottom=0.125,
    top=0.9,
    hspace=0.0,
)

import numpy as np

# Add a twin axis on top
ax_top = ax1.twiny()

# Get main axis limits
x_main_min, x_main_max = ax1.get_xlim()

# Original triton_E range
triton_min, triton_max = min(triton_E), max(triton_E)

# Linear map triton_E to span the full main x-axis
triton_mapped = (triton_E - triton_min) / (triton_max - triton_min) * (x_main_max - x_main_min) + x_main_min

# Stretch twin axis
ax_top.set_xlim(ax1.get_xlim())

# Create ticks every 0.25 in triton_E space
tick_values = np.arange(triton_min, triton_max + 0.25, 0.1)

# Map tick_values to the main axis
tick_positions = (tick_values - triton_min) / (triton_max - triton_min) * (x_main_max - x_main_min) + x_main_min

ax_top.set_xticks(tick_positions)
ax_top.set_xticklabels([f"{t:.2f}" for t in tick_values])

# Example vertical line (scaled)
ax_top.axvline(
    x=x_main_min + (7.8 - triton_min) / (triton_max - triton_min) * (x_main_max - x_main_min),
    color='r',
    linestyle='--'
)
plt.show()













