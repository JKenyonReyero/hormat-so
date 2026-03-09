import subprocess, os, re, shutil
from datetime import datetime
import numpy as np
import json
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from matplotlib.ticker import MaxNLocator


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

def forward(x):
    # Map S2n → Triton energy (descending)
    return triton_max - (
        (x - s2n_min) / (s2n_max - s2n_min)
    ) * (triton_max - triton_min)

def inverse(E):
    # Map Triton energy → S2n
    return s2n_min + (
        (triton_max - E) / (triton_max - triton_min)
    ) * (s2n_max - s2n_min)

#--------------------------------------------------------------
# Set up directories
#--------------------------------------------------------------
script_dir = os.path.dirname(os.path.abspath(__file__))

cwd = os.getcwd()  # directory of wherever the user ran the command

# data_dirpath = str(input("directory path from the current working directory include the / but no .\n"))
# if data_dirpath[-1] != "/":
#     data_dirpath = data_dirpath+"/"
temp=int(input("[0] for 28.53 MeV, [1] for 40 MeV, [2] poster data\n"))
if temp == 0:
    data_dirpath = "/dataTPMLi_28,53MeVold_20260225_120138/"  # temp
    data_dirpath = "/dataTPMLi28MeV33s2n_20260302_135816/"  # temp
    data_dirpath = "/trynewh28MeVTPMLi_20260302_164127/"  # temp
    data_dirpath = "/testnospnohor28E40fm_20260303_185318/"  # temp Li pot rmax=40fm
    # data_dirpath = "/testnospnohor28EP40fm_20260303_185715/"  # temp Pang pot rmax=40fm

    data_dirpath = "/testnospnohor28E_20260304_111456/"  # temp Li pot
    data_dirpath = "/testnospnohor28EP_20260304_111845/"  # temp Pang pot


    # data_dirpath = "/testnospnnhor28E40fm_20260303_190126/"  # temp  Li pot new hormat (used to compare 40fm) rmax=40fm
    # data_dirpath = "/testnospnnhor28EP40fm_20260303_190731/"  # temp  Pang pot new hormat (used to compare 40fm) rmax=40fm

    data_dirpath = "/testnospnnhor28E_20260304_112216/"  # temp  Li pot new hormat (used to compare 40fm)
    data_dirpath = "/testnospnnhor28EP_20260304_112614/"  # temp  Pang pot new hormat (used to compare 40fm)

    data_dirpath = "/try2,5horad28MeV_20260306_150837/"  # temp 2.5fm ho rad 28.53 MeV old hormat
elif temp ==1:
    data_dirpath = "/dataTPMLi40MeVold_20260225_122037/"  # temp
    data_dirpath = "/try2,5horad40MeV_20260306_152919/"  # temp 2.5 fm ho rad 40 MeV old hormat
elif temp == 2:
    data_dirpath = "/posteroldh2,5horad30MeVLi_20260306_160213/"
int_xsec_n_l = read_int_xsec_n_l(cwd+data_dirpath+"integrated_xsecs_n_l.npz")
s2n = read_s2n(cwd+data_dirpath+"s2n.txt")
exit_ecm = read_s2n(cwd+data_dirpath+"exit_ecm.txt")
pnl = read_pnl(cwd+data_dirpath+"pnl.txt")
if temp == 0:
    data_dirpath = "/dataTPMpang28,53MeVold_20260225_120426/"  # temp
    data_dirpath = "/dataTPMPang28MeV33s2n_20260302_160934/"  # temp
    data_dirpath = "/trynewh28MeVTPMPang_20260302_164449/"  # temp
    data_dirpath = "/dataTPMLi_28,53MeVold_20260225_120138/"  # temp
    # data_dirpath = "/testnospnnhor28E40fm_20260303_190126/"  # temp  Li pot these actually are all 40fm by mistake
    # # data_dirpath = "/testnospnnhor28EP40fm_20260303_190731/"  # temp  Pang pot these actually are all 40fm by mistake
    # data_dirpath = "/test40fm_20260303_181346/"  # temp testing 40fm no spin Li old hormat
    # data_dirpath = "/test40fmP_20260303_181947/"  # temp testing 40fm no spin Pang old hormat
    # data_dirpath = "/test40fmLinhor_20260303_182750/"  # temp testing 40fm no spin Li new hormat
    # data_dirpath = "/test40fmPangnhor_20260303_183214/"  # temp testing 40fm no spin Pang new hormat
    data_dirpath = "/dataTPMLi28MeV33s2n_20260302_135816/"  # temp testing 40fm no spin Pang new hormat

    data_dirpath = "/testnospnnhor28E_20260304_112216/"  # temp  Li pot new hormat (used to compare 40fm)
    # data_dirpath = "/testnospnnhor28EP_20260304_112614/"  # temp  Pang pot new hormat (used to compare 40fm)

    # data_dirpath = "/testnospnohor28E40fm_20260303_185318/"  # temp Li pot rmax=40fm
    # data_dirpath = "/testnospnohor28EP40fm_20260303_185715/"  # temp Pang pot rmax=40fm

    data_dirpath = "/testnospnnhor28E40fm_20260303_190126/"  # temp  Li pot new hormat (used to compare 40fm) rmax=40fm
    data_dirpath = "/testnospnnhor28EP40fm_20260303_190731/"  # temp  Pang pot new hormat (used to compare 40fm) rmax=40fm

    data_dirpath = "/try3,0horad28MeV_20260306_151826/"  # temp 2.5 fm ho rad 28.53 MeV old hormat
elif temp ==1:
    data_dirpath = "/dataTPMPang40MeVold_20260225_122313/"  # temp
    data_dirpath = "/try3,0horad40MeV_20260306_152513/"  # temp 3.0 fm ho rad 40 MeV old hormat
elif temp == 2:
    data_dirpath = "/posteroldh2,5horad40MeVLi_20260306_160524/"

int_xsec_n_l2 = read_int_xsec_n_l(cwd+data_dirpath+"integrated_xsecs_n_l.npz")
s2n2 = read_s2n(cwd+data_dirpath+"s2n.txt")
exit_ecm2 = read_s2n(cwd+data_dirpath+"exit_ecm.txt")

pnl2 = read_pnl(cwd+data_dirpath+"pnl.txt")
data_dirpath = "/testpchan_20260303_150557/"  # Li pot
# data_dirpath = "/pchanpang_20260303_154952/"  # Pang pot
pchan_int_xsec_n_l = read_int_xsec_n_l(cwd+data_dirpath+"integrated_xsecs_n_l.npz")
pchan_s2n = read_s2n(cwd+data_dirpath+"s2n.txt")
pchan_exit_ecm = read_s2n(cwd+data_dirpath+"exit_ecm.txt")
pchan_pnl = read_pnl(cwd+data_dirpath+"pnl.txt")

# print(int_xsec_n_l[("0","0","0")])
colours = ["#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999"]

plt.rcParams.update({
    "font.size": 30,
    "axes.labelsize": 38,
    "xtick.labelsize": 28,
    "ytick.labelsize": 28,
    "legend.fontsize": 18,
})

# Legend for models (linestyles)
style_handles = [
    Line2D([0],[0], color="black", lw=3.5, linestyle=":", label="LELO"),
    Line2D([0],[0], color="black", lw=3.5, linestyle="--", label="LENLO"),
    Line2D([0],[0], color="black", lw=3.5, linestyle="-", label="NL"),
]
color_handles_n = [
    Line2D([0],[0], color=colours[i], lw=2, label=f"n={i}")
    for i in range(6)
]

# Legend for l (colors)
color_handles = [
    Line2D([0],[0], color=colours[i], lw=2, label=f"l={i}")
    for i in range(6)
]

x_offset = 0.02 * (max(s2n) - min(s2n))
x_offset2 = 0.02 * (max(s2n2) - min(s2n2))



x_offset = 0.02 * (max(s2n) - min(s2n))

k = 90/(93)  # m_triton / (mtriton+90Zr)
try:
    E_cm = k*float(pnl["elab"])
    triton_E = k*(E_cm+s2n-8.48)  # k(E_p + S2n - triton binding energy)
except:
    E_cm = k*28.53
    triton_E = (E_cm-s2n+8.48)  # Assume eE_p is 28.53MeV

print(triton_E)
fig, (ax1, ax2) = plt.subplots(
    1, 2,
    figsize=(10, 10),
    sharex=True,
    sharey=True
)

for i in range(6):

    # LELO
    line1, = ax1.plot(s2n, int_xsec_n_l[("0","0",f"{i}")], color=colours[i], linestyle=":", linewidth=3)

    # LENLO
    line2, = ax1.plot(s2n, int_xsec_n_l[("1","0",f"{i}")], color=colours[i], linestyle="--", linewidth=3)

    # NL (reference → slightly thicker)
    line3, = ax1.plot(s2n, int_xsec_n_l[("2","0",f"{i}")], color=colours[i], linestyle="-", linewidth=3)

    if temp != 2:
        line_pchan, = ax1.plot(pchan_s2n, pchan_int_xsec_n_l[("2","0",f"{i}")], color="k", linestyle="-", linewidth=3)

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
    line3, = ax2.plot(s2n2, int_xsec_n_l2[("2","0",f"{i}")], color=colours[i], linestyle="-", linewidth=3)

    if temp != 2:
        line_pchan, = ax2.plot(pchan_s2n, pchan_int_xsec_n_l[("2","0",f"{i}")], color="k", linestyle="-", linewidth=3)


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
    0.15, 0.9,                      # x=left, y=slightly above top
    "Li",
    transform=ax1.transAxes,        # <- key line
    ha="center",
    va="bottom",
    fontsize=20,
)

ax2.text(
    0.15, 0.9,                      # x=left, y=slightly above top
    "Pang",
    transform=ax2.transAxes,        # <- key line
    ha="center",
    va="bottom",
    fontsize=20,
)

# Labels only on outer axes (cleaner look)
ax2.set_xlabel(r"$S_{2n}$ (MeV)")
ax1.set_ylabel(r"$\sigma_R$ (mb)")
ax1.set_xlabel(r"$S_{2n}$ (MeV)")

ax1.grid(alpha=0.2)
ax2.grid(alpha=0.2)

plt.subplots_adjust(
    left=0.151,
    right=0.974,
    bottom=0.125,
    top=0.884,
    hspace=0.037,
    wspace=0.037
)

# --- Get main axis limits ---
s2n_min = np.min(s2n)
s2n_max = np.max(s2n)
# --- Triton energy limits ---
triton_min = np.min(exit_ecm)
triton_max = np.max(exit_ecm)

for ax in [ax1, ax2]:
    # --- Create secondary axis on top ---
    ax_top = ax.secondary_xaxis("top", functions=(forward, inverse))
    ax_top.set_xlabel("$E_{CoM}$ (MeV)")

    # --- Use MaxNLocator for clean ticks ---
    ax_top.xaxis.set_major_locator(
        MaxNLocator(nbins=4, integer=True)
    )

    # Optional: make ticks point outward for clarity
    ax_top.tick_params(direction="out")

    # Example vertical line (scaled)
    ax.axvline(
        x=inverse(7.8),
        color='r',
        linestyle='--'
    )
plt.show()





fig, (ax1, ax2) = plt.subplots(
    1, 2,
    figsize=(8, 10),
    sharex=True,
    sharey=True
)

for i in range(4):

    # LELO
    line1, = ax1.plot(s2n, int_xsec_n_l[("0",f"{i}",f"0")], color=colours[i], linestyle=":", linewidth=3)

    # LENLO
    line2, = ax1.plot(s2n, int_xsec_n_l[("1",f"{i}",f"0")], color=colours[i], linestyle="--", linewidth=3)

    # NL (reference → slightly thicker)
    line3, = ax1.plot(s2n, int_xsec_n_l[("2",f"{i}",f"0")], color=colours[i], linestyle="-", linewidth=3)

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
    line1, = ax2.plot(s2n2, int_xsec_n_l2[("0",f"{i}",f"0")], color=colours[i], linestyle=":", linewidth=3)

    # LENLO
    line2, = ax2.plot(s2n2, int_xsec_n_l2[("1",f"{i}",f"0")], color=colours[i], linestyle="--", linewidth=3)

    # NL (reference → slightly thicker)
    line3, = ax2.plot(s2n2, int_xsec_n_l2[("2",f"{i}",f"0")], color=colours[i], linestyle="-", linewidth=3)

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
    loc="upper right",
    frameon=False
)
# Legend on each axis
ax2.legend(
    handles=color_handles_n,
    loc="upper right",
    frameon=False
)

ax1.text(
    0.15, 0.9,                      # x=middle, y=slightly above top
    "Li",
    transform=ax1.transAxes,        # <- key line
    ha="center",
    va="bottom",
    fontsize=20,
)

ax2.text(
    0.15, 0.9,                      # x=middle, y=slightly above top
    "Pang",
    transform=ax2.transAxes,        # <- key line
    ha="center",
    va="bottom",
    fontsize=20,
)

# Labels only on outer axes (cleaner look)
ax2.set_xlabel(r"$S_{2n}$ (MeV)")
ax1.set_ylabel(r"$\sigma_R$ (mb)")
ax1.set_xlabel(r"$S_{2n}$ (MeV)")

ax1.grid(alpha=0.2)
ax2.grid(alpha=0.2)

plt.subplots_adjust(
    left=0.151,
    right=0.974,
    bottom=0.125,
    top=0.884,
    hspace=0.037,
    wspace=0.037
)

# --- Get main axis limits ---
s2n_min = np.min(s2n)
s2n_max = np.max(s2n)
# --- Triton energy limits ---
triton_min = np.min(exit_ecm)
triton_max = np.max(exit_ecm)

for ax in [ax1, ax2]:
    # --- Create secondary axis on top ---
    ax_top = ax.secondary_xaxis("top", functions=(forward, inverse))
    ax_top.set_xlabel("$E_{CoM}$ (MeV)")

    # --- Use MaxNLocator for clean ticks ---
    ax_top.xaxis.set_major_locator(
        MaxNLocator(nbins=4, integer=True)
    )

    # Optional: make ticks point outward for clarity
    ax_top.tick_params(direction="out")

    # Example vertical line (scaled)
    ax.axvline(
        x=inverse(7.8),
        color='r',
        linestyle='--'
    )
plt.show()

# Plot
# fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8,6), sharex=True, sharey=True)
# for i in range(0,3):
#     ax1.plot(s2n, int_xsec_n_l[("0",f"{i}","0")], c=colours[i], linestyle=":")

#     ax1.plot(s2n, int_xsec_n_l[("1",f"{i}","0")], c=colours[i], linestyle="--")

#     NL_line, = ax1.plot(s2n, int_xsec_n_l[("2",f"{i}","0")], c=colours[i], linestyle="-")

#     # --- Direct label using NL curve ---
#     y_end = NL_line.get_ydata()[-1]

#     ax1.text(
#         s2n[-1] + x_offset,
#         y_end,
#         rf"$n={i}$",
#         color=colours[i],
#         fontsize=28,
#         va="center",
#     )

#     ax2.plot(s2n2, int_xsec_n_l2[("0",f"{i}","0")], c=colours[i], linestyle=":")

#     ax2.plot(s2n2, int_xsec_n_l2[("1",f"{i}","0")], c=colours[i], linestyle="--")

#     NL_line, = ax2.plot(s2n2, int_xsec_n_l2[("2",f"{i}","0")], c=colours[i], linestyle="-")

#     # --- Direct label using NL curve ---
#     y_end = NL_line.get_ydata()[-1]

#     ax2.text(
#         s2n2[-1] + x_offset2,
#         y_end,
#         rf"$n={i}$",
#         color=colours[i],
#         fontsize=28,
#         va="center",
#     )

    

# leg = ax1.legend(handles=style_handles, title="Locality", loc="upper right", frameon=False)
# ax1.add_artist(leg)

# ax2.set_xlabel(r"$S_{2n}$ (MeV)", fontsize=16)
# ax1.set_ylabel(r"$\sigma_R$ (mb)", fontsize=16)
# ax1.grid(alpha=0.2)
# ax2.grid(alpha=0.2)

# plt.show()


# Plot
# fig, ax = plt.subplots(figsize=(8,6))
# for i in range(0,5):
#     ax.plot(s2n, int_xsec_n_l[("0",f"0",f"{i}")], label=f"LELO n=0 l={i}", c=colours[i], linestyle=":")

#     ax.plot(s2n, int_xsec_n_l[("1",f"0",f"{i}")], label=f"LENLO n=0 l={i}", c=colours[i], linestyle="--")

#     ax.plot(s2n, int_xsec_n_l[("2",f"0",f"{i}")], label=f"NL n=0 l={i}", c=colours[i], linestyle="-")

# ax.set_xlabel(r"$S_{2n}$ (MeV)", fontsize=16)
# ax.set_ylabel(r"$\sigma_R$ (mb)", fontsize=16)
# ax.legend(fontsize=12)
# ax.grid(alpha=0.2)

# plt.show()


fig, (ax1, ax2) = plt.subplots(
    1, 2,
    figsize=(8, 10),
    sharex=True,
    sharey=True
)

for i in range(0,3):
    # ax.plot(s2n, int_xsec_n_l[("0",f"{i}","0")], label=f"NL n{i} l0")

    ax1.plot(s2n, int_xsec_n_l[("1",f"{i}","0")]/int_xsec_n_l[("2",f"{i}","0")], label=f"LENLO/NL n={i} l=0", linestyle="--", linewidth=3)

    ax1.plot(s2n, int_xsec_n_l[("0",f"{i}","0")]/int_xsec_n_l[("2",f"{i}","0")], label=f"LELO/NL n={i} l=0", linestyle=":", linewidth=3)

    #pchan
    ax1.plot(s2n, pchan_int_xsec_n_l[("2",f"{i}","0")]/int_xsec_n_l[("2",f"{i}","0")], label=f"pChan/NL", linestyle="-", linewidth=3)

    ax2.plot(s2n2, int_xsec_n_l2[("1",f"{i}","0")]/int_xsec_n_l2[("2",f"{i}","0")], label=f"LENLO/NL n={i} l=0", linestyle="--", linewidth=3)

    ax2.plot(s2n2, int_xsec_n_l2[("0",f"{i}","0")]/int_xsec_n_l2[("2",f"{i}","0")], label=f"LELO/NL n={i} l=0", linestyle=":", linewidth=3)

    #pchan
    ax2.plot(s2n2, pchan_int_xsec_n_l[("2",f"{i}","0")]/int_xsec_n_l2[("2",f"{i}","0")], label=f"pChan/NL", linestyle="-", linewidth=3)


# Legend on each axis
ax1.legend(
    handles=style_handles,
    loc="upper right",
    frameon=False
)
# Legend on each axis
# ax2.legend(
#     handles=color_handles_n,
#     loc="upper right",
#     frameon=False
# )

ax1.text(
    0.15, 0.9,                      # x=middle, y=slightly above top
    "Li",
    transform=ax1.transAxes,        # <- key line
    ha="center",
    va="bottom",
    fontsize=20,
)

ax2.text(
    0.15, 0.9,                      # x=middle, y=slightly above top
    "Pang",
    transform=ax2.transAxes,        # <- key line
    ha="center",
    va="bottom",
    fontsize=20,
)

# Labels only on outer axes (cleaner look)
ax2.set_xlabel(r"$S_{2n}$ (MeV)")
ax1.set_ylabel(r"$\sigma_R$ ratio (dimensionless)")
ax1.set_xlabel(r"$S_{2n}$ (MeV)")

ax1.grid(alpha=0.2)
ax2.grid(alpha=0.2)

ax1.axhline(y=1.0, c="k", alpha=0.25)
ax2.axhline(y=1.0, c="k", alpha=0.25)

plt.subplots_adjust(
    left=0.198,
    right=0.974,
    bottom=0.125,
    top=0.884,
    hspace=0.037,
    wspace=0.037
)

# --- Get main axis limits ---
s2n_min = np.min(s2n)
s2n_max = np.max(s2n)
# --- Triton energy limits ---
triton_min = np.min(exit_ecm)
triton_max = np.max(exit_ecm)

for ax in [ax1, ax2]:
    # --- Create secondary axis on top ---
    ax_top = ax.secondary_xaxis("top", functions=(forward, inverse))
    ax_top.set_xlabel("$E_{CoM}$ (MeV)")

    # --- Use MaxNLocator for clean ticks ---
    ax_top.xaxis.set_major_locator(
        MaxNLocator(nbins=4, integer=True)
    )

    # Optional: make ticks point outward for clarity
    ax_top.tick_params(direction="out")

    # Example vertical line (scaled)
    ax.axvline(
        x=inverse(7.8),
        color='r',
        linestyle='--'
    )
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


leg1 = ax.legend(handles=color_handles, loc="upper left", frameon=False)
leg2 = ax.legend(handles=style_handles, loc="upper right", frameon=False)

ax.add_artist(leg1)

ax.set_xlabel(r"$S_{2n}$ (MeV)")
ax.set_ylabel(r"$\sigma_R$ (mb)")
ax.grid(alpha=0.2)

plt.show()

















