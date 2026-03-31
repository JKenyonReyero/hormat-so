# Imports
import matplotlib
from matplotlib import pyplot as plt
import numpy as np
from scipy.integrate import simpson, trapezoid
from scipy.interpolate import interp1d
from matplotlib.ticker import MultipleLocator
import os, json

def read_int_xsec(dir, file):
    data = np.loadtxt(dir+file)
    x = data[:,0]
    y = data[:,1]
    return x, y

def calc_spreading_func(E_grid, Ea):
    ep0 = 19.4   # Epsilon 0 for gamma in MeV
    ep1 = 1.40   # Epsilon 1 for gamma in MeV
    E0 = 18.4    # E0 for gamma in MeV
    E1 = 1.60    # E1 for gamma in MeV

    gamma = np.zeros(len(Ea))
    f = np.zeros((len(Ea),len(E_grid)))

    #! spread
    for i in range(0,len(Ea)):
        # gamma func
        first = (ep0*(Ea[i]**2))/(Ea[i]**2 + E0**2)
        second = (ep1*(Ea[i]**2))/(Ea[i]**2 + E1**2)
        gamma[i] = ((first + second)/2)+0.075

        # spread
        for j in range(1,len(E_grid)):
            # 75 kev added to account for experimental error
            bottom = (Ea[i]-E_grid[j])**2 + ((gamma[i])**2)/4
            f[i][j] = (2*np.pi)**-1 * (gamma[i])/bottom
    return f

def spread(xsec, Ea, f):
    #! spread
    for i in range(0,len(Ea)):
        # integrate
        area = simpson(f[i][:], dx=dx)
        # normalise
        n0 = xsec[i]/area
        f[i] = f[i]*n0

        # area_norm = simpson(f[i][:], dx=dx)
        # print(f"A: {area} xsec: {(xsec[i])} n0: {n0} A_norm: {area_norm} Gamma: {gamma[i]}")

    total = np.sum(f, axis=0)
    return f, total

def read_spread(dir, E_grid):
    E_ex, int_xsecs_NL = read_int_xsec(dir+"/NL/", "int_xsecs")
    E_ex, int_xsecs_LENLO = read_int_xsec(dir+"/LENLO/", "int_xsecs")
    E_ex, int_xsecs_LELO = read_int_xsec(dir+"/LELO/", "int_xsecs")

    base_f = calc_spreading_func(E_grid, E_ex)

    f_NL, total_NL = spread(int_xsecs_NL, E_ex, base_f)
    f_LENLO, total_LENLO = spread(int_xsecs_LENLO, E_ex, base_f)
    f_LELO, total_LELO = spread(int_xsecs_LELO, E_ex, base_f)
    return E_ex, (f_NL, total_NL), (f_LENLO, total_LENLO), (f_LELO, total_LELO)

def read_s2n(file_path):
    with open(file_path, "r") as f:
        arr = np.array(json.load(f))
    return arr

def spin_dist_fig(E_grid, total_1=tuple, total_2=tuple, savesuffixE=str):
    fig, (ax1,ax2) = plt.subplots(
        2, 1,
        figsize=(8, 10),
        sharex=True
    )
    plt.subplots_adjust(
        left=0.105,
        right=0.988,
        bottom=0.09,
        top=0.976,
        hspace=0.037,
        wspace=0.037
    )

    # li 30
    ax1.plot(E_grid, total_1[0], c="k", linestyle="-", label="NL", linewidth=3)
    ax1.plot(E_grid, total_1[1], c="k", linestyle="--", label="LENLO", linewidth=3)
    ax1.plot(E_grid, total_1[2], c="k", linestyle=":", label="LELO", linewidth=3)
    # pang 30
    ax1.plot(E_grid, total_2[0], c="b", linestyle="-", label="NL", linewidth=3)
    ax1.plot(E_grid, total_2[1], c="b", linestyle="--", label="LENLO", linewidth=3)
    ax1.plot(E_grid, total_2[2], c="b", linestyle=":", label="LELO", linewidth=3)

    # li 30
    ax2.plot(E_grid, total_1[1]/total_1[0], c="k", linestyle="--", label="LENLO", linewidth=3)
    ax2.plot(E_grid, total_1[2]/total_1[0], c="k", linestyle=":", label="LELO", linewidth=3)
    # pang 30
    ax2.plot(E_grid, total_2[1]/total_2[0], c="b", linestyle="--", label="LENLO", linewidth=3)
    ax2.plot(E_grid, total_2[2]/total_2[0], c="b", linestyle=":", label="LELO", linewidth=3)

    ax2.axhline(y=1.0, color="k", alpha=0.25)

    # Axis labels
    ax1.set_ylabel(r"$\sigma$ [mb]")
    ax2.set_ylabel(r"$\sigma_{LE}/\sigma_{NL}$")
    ax2.set_xlabel(r"$E_{ex}$ (MeV)")

    # if li_or_pang =="Li":
    #     ax1.set_xlim((0,max(E_ex)))
    #     ax1.set_ylim((0,8))
    #     ax2.yaxis.set_major_locator(MultipleLocator(0.01))
    # elif li_or_pang == "Pang":
    #     ax1.set_xlim((0,max(E_ex)))
    #     ax1.set_ylim((0,25))
    #     ax2.yaxis.set_major_locator(MultipleLocator(0.05))

    plt.tight_layout()
    ax1.legend()
    ax2.legend()
    ax1.grid(alpha=0.2)
    ax2.grid(alpha=0.2)
    plt.savefig(f"spin_dist_diff_lipang{savesuffixE}MeV.svg", pad_inches=0.02, transparent=True)
    plt.show()

def spin_dist_fig_big(E_grid, s2n_30, s2n_40, exit_ecm_30, exit_ecm_40, total_1=tuple, total_2=tuple, total_3=tuple, total_4=tuple, savesuffixE=str):
    # fig, axes = plt.subplots(
    #     2, 2,
    #     figsize=(16, 10),
    #     sharex=True
    # )
    fig = plt.figure(figsize=(16,10))

    # Top row (share y with each other)
    ax1 = fig.add_subplot(2,2,1)
    ax2 = fig.add_subplot(2,2,2, sharex=ax1, sharey=ax1)

    # Bottom row (share y with each other, but not with top row)
    ax3 = fig.add_subplot(2,2,3, sharex=ax1)
    ax4 = fig.add_subplot(2,2,4, sharex=ax1, sharey=ax3)

    # Create interpolation functions
    # 'fill_value="extrapolate"' allows zooming outside original x range
    exit_ecm_of_s2n = interp1d(s2n_30, exit_ecm_30, kind='linear', fill_value="extrapolate")
    s2n_of_exit_ecm = interp1d(exit_ecm_30, s2n_30, kind='linear', fill_value="extrapolate")

    # Create a secondary x-axis on top
    def forward(s2nval):
        return exit_ecm_of_s2n(s2nval)

    def inverse(ecmval):
        return s2n_of_exit_ecm(ecmval)

    # --- Create secondary axis on top ---
    ax_top1 = ax1.secondary_xaxis("top", functions=(forward, inverse))
    ax_top1.set_xlabel("$E_t$ (MeV)", labelpad=12)

    # Create interpolation functions
    # 'fill_value="extrapolate"' allows zooming outside original x range
    exit_ecm_of_s2n40 = interp1d(s2n_40, exit_ecm_40, kind='linear', fill_value="extrapolate")
    s2n_of_exit_ecm40 = interp1d(exit_ecm_40, s2n_40, kind='linear', fill_value="extrapolate")
    # Create a secondary x-axis on top
    def forward2(s2nval):
        return exit_ecm_of_s2n40(s2nval)

    def inverse2(ecmval):
        return s2n_of_exit_ecm40(ecmval)

    # --- Create secondary axis on top ---
    ax_top2 = ax2.secondary_xaxis("top", functions=(forward2, inverse2))
    ax_top2.set_xlabel("$E_t$ (MeV)", labelpad=12)

    ax1.tick_params(labelbottom=False)
    ax2.tick_params(labelleft=False, labelbottom=False)
    ax4.tick_params(labelleft=False)   
    plt.subplots_adjust(
        left=0.1,
        right=0.988,
        bottom=0.11,
        top=0.886,
        hspace=0.0,  # 0.037
        wspace=0.0  # 0.037
    )
    # li 30
    ax1.plot(E_grid, total_1[0], c="k", linestyle="-", label="Li NL", linewidth=3)
    ax1.plot(E_grid, total_1[1], c="k", linestyle="--", label="Li LENLO", linewidth=3)
    ax1.plot(E_grid, total_1[2], c="k", linestyle=":", label="Li LELO", linewidth=3)
    # pang 30
    ax1.plot(E_grid, total_2[0], c="b", linestyle="-", label="Pang NL", linewidth=3)
    ax1.plot(E_grid, total_2[1], c="b", linestyle="--", label="Pang LENLO", linewidth=3)
    ax1.plot(E_grid, total_2[2], c="b", linestyle=":", label="Pang LELO", linewidth=3)

    # li 40
    ax2.plot(E_grid, total_3[0], c="r", linestyle="-", label="Li NL", linewidth=3)
    ax2.plot(E_grid, total_3[1], c="r", linestyle="--", label="Li LENLO", linewidth=3)
    ax2.plot(E_grid, total_3[2], c="r", linestyle=":", label="Li LELO", linewidth=3)
    # pang 40
    ax2.plot(E_grid, total_4[0], c="g", linestyle="-", label="Pang NL", linewidth=3)
    ax2.plot(E_grid, total_4[1], c="g", linestyle="--", label="Pang LENLO", linewidth=3)
    ax2.plot(E_grid, total_4[2], c="g", linestyle=":", label="Pang LELO", linewidth=3)

    # li 30
    ax3.plot(E_grid, total_1[1]/total_1[0], c="k", linestyle="--", label="Li LENLO/NL", linewidth=3)
    ax3.plot(E_grid, total_1[2]/total_1[0], c="k", linestyle=":", label="Li LELO/NL", linewidth=3)
    # pang 30
    ax3.plot(E_grid, total_2[1]/total_2[0], c="b", linestyle="--", label="Pang LENLO/NL", linewidth=3)
    ax3.plot(E_grid, total_2[2]/total_2[0], c="b", linestyle=":", label="Pang LELO/NL", linewidth=3)

    # li 40
    ax4.plot(E_grid, total_3[1]/total_3[0], c="r", linestyle="--", label="Li LENLO/NL", linewidth=3)
    ax4.plot(E_grid, total_3[2]/total_3[0], c="r", linestyle=":", label="Li LELO/NL", linewidth=3)
    # pang 40
    ax4.plot(E_grid, total_4[1]/total_4[0], c="g", linestyle="--", label="Pang LENLO/NL", linewidth=3)
    ax4.plot(E_grid, total_4[2]/total_4[0], c="g", linestyle=":", label="Pang LELO/NL", linewidth=3)

    ax3.axhline(y=1.0, color="k", alpha=0.25)
    ax4.axhline(y=1.0, color="k", alpha=0.25)

    ax1.set_ylim(ymin=-0.1, ymax=30)
    ax1.set_xlim(xmin=-0.1, xmax=24)

    # Axis labels
    ax1.set_ylabel(r"$\sigma$ [mb]")
    ax3.set_ylabel(r"$\sigma_{LE}/\sigma_{NL}$")
    ax3.set_xlabel(r"$E_{ex}$ (MeV)")
    ax4.set_xlabel(r"$E_{ex}$ (MeV)")

    # if li_or_pang =="Li":
    #     ax1.set_xlim((0,max(E_ex)))
    #     ax1.set_ylim((0,8))
    #     ax2.yaxis.set_major_locator(MultipleLocator(0.01))
    # elif li_or_pang == "Pang":
    #     ax1.set_xlim((0,max(E_ex)))
    #     ax1.set_ylim((0,25))
    #     ax2.yaxis.set_major_locator(MultipleLocator(0.05))

    ax3.yaxis.set_major_locator(MultipleLocator(0.05))


    # plt.tight_layout()
    ax1.legend(title=r"$E_p = $30MeV")
    ax2.legend(title=r"$E_p = $40MeV")
    ax3.legend(title=r"$E_p = $30MeV", loc="upper left", bbox_to_anchor=(0.6, 0.925))
    ax4.legend(title=r"$E_p = $40MeV", loc="upper left", bbox_to_anchor=(0.6, 0.925))
    ax1.grid(alpha=0.2)
    ax2.grid(alpha=0.2)
    ax3.grid(alpha=0.2)
    ax4.grid(alpha=0.2)
    plt.savefig(f"spin_dist_diff_lipang{savesuffixE}MeV.svg", pad_inches=0.02, transparent=True)
    plt.show()

# variables
jpi = ["0+","2+","4+","6+","3-","0+","5-","1-","3-","1-","3-","2+","4+","6+","8+","5-","1-","7-","3-","5-","3-","3-","5-","5-","0+","0+","2+","2+","2+","4+","1-","4+","0+","2+","3-","2+","4+","5-","2+","7-","2+"]

cwd = os.getcwd()  # directory of wherever the user ran the command
# li_or_pang = int(input("[1] Li or [2] Pang\n"))
# which_E = int(input("[1] 30MeV [2] 40MeV\n"))
# if li_or_pang == 1 and which_E == 1:
#     data_dirpath = "/poster_spin_distLi30MeV_20260310_125333/"
#     li_or_pang = "Li"
# elif li_or_pang == 2 and which_E == 1:
#     data_dirpath = "/poster_spin_distPang30MeV_20260310_153537/"
#     li_or_pang = "Pang"
# elif li_or_pang == 1 and which_E == 2:
#     data_dirpath = "/poster_spin_distLi40MeV_20260311_160735/"
#     li_or_pang = "Li"
# elif li_or_pang == 2 and which_E == 2:
#     data_dirpath = "/poster_spin_distPang40MeV_20260311_160924/"
#     li_or_pang = "Pang"
# else:
#     print("----Wrong----")
# E_ex, int_xsecs_NL = read_int_xsec(cwd+data_dirpath+"/NL/", "int_xsecs")
# E_ex, int_xsecs_LENLO = read_int_xsec(cwd+data_dirpath+"/LENLO/", "int_xsecs")
# E_ex, int_xsecs_LELO = read_int_xsec(cwd+data_dirpath+"/LELO/", "int_xsecs")

dx = 0.001
E_min = -40
E_max = 40
E_grid = np.arange(E_min, E_max+dx, dx)
data_dirpath_li30 = "/poster_spin_distLi30MeV_20260317_121643/"
data_dirpath_pang30 = "/poster_spin_distPang30MeV_20260317_121914/"
data_dirpath_li40 = "/poster_spin_distLi40MeV_20260317_122009/"
data_dirpath_pang40 = "/poster_spin_distPang40MeV_20260317_122054/"

data_dirpath_li30 = "/poster_spin_distLi30MeV_20260317_121643/"
data_dirpath_pang30 = "/poster_spin_distPang30MeV_20260317_121914/"
data_dirpath_li40 = "/nlcorli30_20260320_161228/"
data_dirpath_pang40 = "/nlcorpang30_20260320_161815/"

data_dirpath_li30 = "/poster_spin_distLi40MeV_20260317_122009/"
data_dirpath_pang30 = "/poster_spin_distPang40MeV_20260317_122054/"
data_dirpath_li40 = "/nlcorli40_20260320_162833/"
data_dirpath_pang40 = "/nlcorpang40_20260320_162903/"

E_ex_li30, (f_NL_li30, total_NL_li30), (f_LENLO_li30, total_LENLO_li30), (f_LELO_li30, total_LELO_li30) = read_spread(cwd+data_dirpath_li30, E_grid)

s2n_30 = read_s2n(cwd+data_dirpath_li30+"s2n.txt")
E_t_30 = read_s2n(cwd+data_dirpath_li30+"exit_ecm.txt")
E_ex_pang30, (f_NL_pang30, total_NL_pang30), (f_LENLO_pang30, total_LENLO_pang30), (f_LELO_pang30, total_LELO_pang30) = read_spread(cwd+data_dirpath_pang30, E_grid)

E_ex_li40, (f_NL_li40, total_NL_li40), (f_LENLO_li40, total_LENLO_li40), (f_LELO_li40, total_LELO_li40) = read_spread(cwd+data_dirpath_li40, E_grid)

s2n_40 = read_s2n(cwd+data_dirpath_li40+"s2n.txt")
E_t_40 = read_s2n(cwd+data_dirpath_li40+"exit_ecm.txt")
E_ex_pang40, (f_NL_pang40, total_NL_pang40), (f_LENLO_pang40, total_LENLO_pang40), (f_LELO_pang40, total_LELO_pang40) = read_spread(cwd+data_dirpath_pang40, E_grid)

plt.rcParams["font.size"] = 30
plt.rcParams["axes.labelsize"] = 38
plt.rcParams["axes.titlesize"] = 22
plt.rcParams["xtick.labelsize"] = 28
plt.rcParams["ytick.labelsize"] = 28
plt.rcParams["legend.fontsize"] = 16
plt.rcParams["legend.title_fontsize"] = 18


# fig, (ax1, ax2, ax3) = plt.subplots(
#     3, 1,
#     figsize=(10,8),
#     sharex=True,
#     sharey=True
# )

# for i in range(0,len(E_ex)):
#     ax1.plot(E_grid, f_NL[i], label=jpi[i])
#     ax2.plot(E_grid, f_LENLO[i], label=jpi[i])
#     ax3.plot(E_grid, f_LELO[i], label=jpi[i])

# ax1.plot(E_grid, total_NL, c="k")
# ax2.plot(E_grid, total_LENLO, c="k")
# ax3.plot(E_grid, total_LELO, c="k")

# ax1.set_title("NL")
# ax2.set_title("LENLO")
# ax3.set_title("LELO")

# # Axis labels
# ax1.set_ylabel(r"$\sigma$ [mb]")
# ax1.set_xlabel("Excitation Energy (MeV)")
# ax2.set_xlabel("Excitation Energy (MeV)")
# ax3.set_xlabel("Excitation Energy (MeV)")

# # Optional grids
# ax1.grid(True, alpha=0.25)
# ax2.grid(True, alpha=0.25)
# ax3.grid(True, alpha=0.25)

# # Tight layout
# plt.tight_layout()
# plt.legend()

# plt.show()


# 30 MeV
spin_dist_fig(E_grid, total_1=(total_NL_li30, total_LENLO_li30, total_LELO_li30),
              total_2=(total_NL_pang30, total_LENLO_pang30, total_LELO_pang30), savesuffixE="30")
# 40MeV
spin_dist_fig(E_grid, total_1=(total_NL_li40, total_LENLO_li40, total_LELO_li40),
              total_2=(total_NL_pang40, total_LENLO_pang40, total_LELO_pang40), savesuffixE="30")
# 30 and 40
spin_dist_fig_big(E_grid, s2n_30, s2n_40, E_t_30, E_t_40, total_1=(total_NL_li30, total_LENLO_li30, total_LELO_li30),
              total_2=(total_NL_pang30, total_LENLO_pang30, total_LELO_pang30), 
              total_3=(total_NL_li40, total_LENLO_li40, total_LELO_li40),
              total_4=(total_NL_pang40, total_LENLO_pang40, total_LELO_pang40),
              savesuffixE="30and40")

