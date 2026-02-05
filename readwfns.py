import numpy as np
import re
import matplotlib.pyplot as plt
import itertools

def read_wavefunctions(filename):
    data = {}  # dictionary: L -> array([[r, Re, Im], ...])
    current_L = None
    block = []

    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            # Detect new block
            if line.startswith("# L=") or line.startswith('# l='):
                if current_L is not None and block:
                    data[current_L] = np.array(block, dtype=float)
                    block = []
                current_L = int(re.findall(r"\d+", line)[0])

            # Parse data lines
            else:
                # Convert Fortran 'D' exponents to Python 'E'
                line = line.replace("D", "E")
                try:
                    vals = [float(x) for x in line.split()]
                    if len(vals) == 3:
                        block.append(vals)
                except ValueError:
                    continue

        # Save last block
        if current_L is not None and block:
            data[current_L] = np.array(block, dtype=float)

    return data


# Example usage:
# wavefunctions = read_wavefunctions("readperey_20250916_154400/pChanle/wfns_proton.readperey")
wavefunctions = read_wavefunctions("fort.211")
# wavefunctions2 = read_wavefunctions("readperey_20250916_154400/pChannl/wfns_proton.readperey")
# wavefunctions2 = read_wavefunctions("readpereyn_20250916_172841/pChanle/wfns_proton.readpereyn")
# wavefunctions = read_wavefunctions("readperen2_20250916_173251/pChanle/wfns_proton.readperen2")
# wavefunctions2 = read_wavefunctions("waves_in.testulocn")
wavefunctions2 = read_wavefunctions("wfnsgreg.16")


print(wavefunctions.keys())
print(wavefunctions2.keys())

# Access block for L=0:
print("L=0 data shape:", wavefunctions[0].shape)
print("First rows of L=0:\n", wavefunctions[0][:5])


import matplotlib.pyplot as plt

# Choose L to plot
L_plot = 0

# Extract data
data = wavefunctions[L_plot]
r = data[:, 0]        # radius
Re = data[:, 1]       # real part
Im = data[:, 2]       # imaginary part

# Create the plot
plt.figure(figsize=(8,5))
plt.plot(r, Re, label='Re', color='blue')
plt.plot(r, Im, label='Im', color='red')
plt.xlabel('Radius r')
plt.ylabel('Wavefunction')
plt.title(f'Wavefunction for L={L_plot}')
plt.legend()
plt.grid(True)
plt.show()


import matplotlib.pyplot as plt

# Define which L values to plot (every 3rd from 0 to 18)
Ls_to_plot = list(range(0, 10, 3))

plt.figure(figsize=(10,6))

for L in Ls_to_plot:
    if L in wavefunctions:  # make sure this L exists in your data
        data = wavefunctions[L]
        r = data[:, 0]        # radius
        Re = data[:, 1]       # real part
        plt.plot(r, Re, label=f'L={L}')

plt.xlabel('Radius r')
plt.ylabel('Real part of wavefunction')
plt.title('Real component of wavefunctions for every 3rd angular momentum')
plt.legend()
plt.grid(True)
plt.show()

import matplotlib.pyplot as plt
import numpy as np
import itertools

# Define which L values to plot
Ls_to_plot = list(range(0, 1, 3))

# Create figure with 3 subplots
fig, axes = plt.subplots(3, 1, figsize=(10, 15), sharex=True)

# Get the default color cycle
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
color_iter = itertools.cycle(colors)

for L in Ls_to_plot:
    if L in wavefunctions:
        c = next(color_iter)
        data = wavefunctions[L]
        data2 = wavefunctions2[L]
        r = data[:, 0]
        Re = data[:, 1]
        Im = data[:, 2]
        density = Re**2 + Im**2  # |psi|^2

        r2 = data2[:, 0]
        Re2 = data2[:, 1]
        Im2 = data2[:, 2]
        density2 = Re2**2 + Im2**2  # |psi|^2

        # Normalize the probability density
        norm = np.trapezoid(density, r)
        density_normalized = density / norm

        # Normalize the probability density
        norm2 = np.trapezoid(density2, r)
        density_normalized2 = density2 / norm2

        # Plot real, imaginary, and normalized |psi|^2
        axes[0].plot(r, Re, label=f'L={L} LE NLO', c=c)
        axes[1].plot(r, Im, label=f'L={L} LE NLO', c=c)
        axes[2].plot(r, density_normalized, label=f'L={L} LE NLO', c=c)

        # Plot real, imaginary, and normalized |psi|^2
        axes[0].plot(r, Re2, label=f'L={L} NL', linestyle="dotted", c=c)
        axes[1].plot(r, Im2, label=f'L={L} NL', linestyle="dotted", c=c)
        axes[2].plot(r, density_normalized2, label=f'L={L} NL', linestyle="dotted", c=c)

# Set labels, titles, and legends
axes[0].set_ylabel('Re(ψ)')
axes[0].set_title('Real part of wavefunction')
axes[0].legend()
axes[0].grid(True)

axes[1].set_ylabel('Im(ψ)')
axes[1].set_title('Imaginary part of wavefunction')
axes[1].legend()
axes[1].grid(True)

axes[2].set_xlabel('Radius r')
axes[2].set_ylabel('|ψ|² (normalized)')
axes[2].set_title('Normalized probability density |ψ|²')
axes[2].legend()
axes[2].grid(True)

plt.tight_layout()
plt.show()

