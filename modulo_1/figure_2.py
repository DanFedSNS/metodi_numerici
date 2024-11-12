import numpy as np
import matplotlib.pyplot as plt

# Array of different L values to consider
L_array = np.linspace(10, 30, 5, dtype=int)
model = "ising2d_hex_cluster"
colors = plt.get_cmap('tab10') 

# Create the layout for the subplots
fig, ax = plt.subplots(2, 2, figsize=(12, 10))

for i, L in enumerate(L_array):
    filepath = f'./analysis_{model}/L{L}.dat'
    data = np.loadtxt(filepath, delimiter=",")

    beta = data[:, 0]  
    specific_heat = data[:, 1]  # Specific heat
    susceptibility = data[:, 2]  # Susceptibility
    magn_abs_avg = data[:, 3]  # Average absolute magnetization
    energy_avg = data[:, 4]  # Average energy

    # Plot for specific heat
    ax[0, 0].plot(beta, specific_heat, color=colors(i), label=f'L={L}', marker='o', linestyle='none')
    ax[0, 0].set_title("Specific Heat")
    ax[0, 0].set_xlabel("Beta")
    ax[0, 0].set_ylabel("Specific Heat")
    ax[0, 0].legend()

    # Plot for susceptibility
    ax[0, 1].plot(beta, susceptibility, color=colors(i), label=f'L={L}', marker='o', linestyle='none')
    ax[0, 1].set_title("Susceptibility")
    ax[0, 1].set_xlabel("Beta")
    ax[0, 1].set_ylabel("Susceptibility")
    ax[0, 1].legend()

    # Plot for average absolute magnetization
    ax[1, 0].plot(beta, magn_abs_avg, color=colors(i), label=f'L={L}', marker='o', linestyle='none')
    ax[1, 0].set_title("Average Absolute Magnetization")
    ax[1, 0].set_xlabel("Beta")
    ax[1, 0].set_ylabel("Average Absolute Magnetization")
    ax[1, 0].legend()

    # Plot for average energy
    ax[1, 1].plot(beta, energy_avg, color=colors(i), label=f'L={L}', marker='o', linestyle='none')
    ax[1, 1].set_title("Average Energy")
    ax[1, 1].set_xlabel("Beta")
    ax[1, 1].set_ylabel("Average Energy")
    ax[1, 1].legend()

# Adjust the layout and save the figure
plt.tight_layout()
plt.savefig('./figure_2.pdf', format='pdf')
