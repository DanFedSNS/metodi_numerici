import numpy as np
import matplotlib.pyplot as plt
import os
import sys

# Parametri
L_array = [120]
beta_couple = [
    [0.43, 0.4357, 0.4429, 0.45],   # Per il modello "ising2d_sq_cluster"
    [0.43, 0.4357, 0.4429, 0.45], # Per il modello "ising2d_tri_cluster"
    [0.43, 0.4357, 0.4429, 0.45]    # Per il modello "ising2d_hex_cluster"
]
modelli = ["ising2d_sq_cluster", "ising2d_tri_cluster", "ising2d_hex_cluster"]

def load_params(filepath):
    params = {}
    with open(filepath, 'r') as f:
        for line in f:
            key, value = line.strip().split('=')
            if key in ['fig_width', 'fig_height', 'dpi']:
                params[key] = float(value)
            else:
                params[key] = float(value) if '.' in value else int(value)
    return params

os.chdir(os.path.dirname(os.path.abspath(__file__)))
params = load_params('params.txt')

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# Creiamo una figura con 4 subplot orizzontali
fig, ax = plt.subplots(1, 4, figsize=(4*params['fig_width'], params['fig_height']))

colori = plt.get_cmap('tab10')

# Definire i limiti dell'asse x
x_limits_magnetization = [-1, 1]

# Loop sui valori di beta
for j in range(4):
    # Loop sui modelli
    for k, modello in enumerate(modelli):
        beta = beta_couple[k][j]
        # Loop sui valori di L
        for i in range(len(L_array)):
            L = L_array[i]
            valori_mag = np.linspace(-1, 1, L**2 + 1)

            filepath = f'./{modello}/L{L}_beta{(beta):.4f}.dat'

            data = np.loadtxt(filepath, skiprows=1, delimiter=",")
            magnetization = data[:, 0]

            atol = 60 / (L**2)
            # Calcolo della frequenza delle occorrenze di magnetizzazione
            magnetization_frequencies = [
                np.count_nonzero(np.isclose(magnetization, val, atol=atol)) for val in valori_mag
            ]
            
            # Plot della distribuzione di magnetizzazione
            label = f'{modello.split("_")[1]} (L={L})'
            ax[j].plot(
                valori_mag, magnetization_frequencies,
                color=colori(k), label=label,
                marker='none', linewidth=params['line_width_axes']
            )

    for spine in ax[j].spines.values():
        spine.set_linewidth(params['line_width_axes'])

    ax[j].set_xlabel('Magnetization per site', fontsize=params['font_size_axis'], labelpad=params['label_pad'])
    ax[j].set_ylabel('Occurrences', fontsize=params['font_size_axis'], labelpad=params['label_pad'])
    ax[j].set_xlim(x_limits_magnetization)

    handles, labels = ax[j].get_legend_handles_labels()
    ax[j].legend(handles, labels, loc='best', fontsize=params['font_size_legend'])
    ax[j].grid(True, which='minor', linestyle=':', linewidth=params['line_width_grid_minor'])
    ax[j].grid(True, which='major', linestyle='--', linewidth=params['line_width_grid_major'])

    ax[j].tick_params(axis='x', labelsize=params['font_size_ticks'], 
               width=params['line_width_axes'], direction='in')
    ax[j].tick_params(axis='y', labelsize=params['font_size_ticks'], 
               width=params['line_width_axes'], direction='in')
    ax[j].set_xticks(ticks=[-1, -0.5, 0, 0.5, 1])
    ax[j].margins(x=0.00, y=0.00)

# Migliora la spaziatura tra i subplot
fig.tight_layout(pad=params['pad'])
plt.savefig('./figure_1.pdf', format='pdf', dpi=int(params['dpi']))
plt.close(fig)
