import numpy as np
import matplotlib.pyplot as plt
import os
import sys

# Array of different L values to consider
L_array = np.linspace(70, 150, 5, dtype=int)
models = ["ising2d_tri_cluster", "ising2d_sq_cluster", "ising2d_hex_cluster"]

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

#plt.rc('text', usetex=True)
plt.rc('font', family='serif')

colors = plt.get_cmap('tab10') 
# Dizionario con i valori di beta_c per i tre modelli
beta_c_values_anal = {
    "ising2d_hex_cluster": np.log(2 + np.sqrt(3)) / 2,
    "ising2d_tri_cluster": np.log(np.sqrt(3)) / 2,
    "ising2d_sq_cluster": np.log(1 + np.sqrt(2)) / 2
}
beta_c_values = {
    "ising2d_hex_cluster": 0.65847,
    "ising2d_tri_cluster": 0.27468,
    "ising2d_sq_cluster": 0.44069
}

# Creiamo una figura con 3 subplot verticali
fig, axes = plt.subplots(3, 1, figsize=(params['fig_width'], 2.25*params['fig_height']))

for ax, model in zip(axes, models):
    for i, L in enumerate(L_array):
        filepath = f'./data/analysis_{model}/L{L}.dat'
        data = np.genfromtxt(filepath, delimiter=" ", dtype=float, filling_values=np.nan)


        beta_c = beta_c_values[model]
        beta = data[:, 0]
        beta = (beta - beta_c) * L
        susceptibility = data[:, 7] / L ** 1.75  # Susceptibility

        ax.plot(beta[::2], susceptibility[::2], color=colors(i), label=f'L={L}', marker='s', linestyle='none',
                markerfacecolor='white', markeredgewidth=params['line_width_axes'], zorder=2)

    # Impostazioni estetiche per ogni subplot
    for spine in ax.spines.values():
        spine.set_linewidth(params['line_width_axes'])
    ax.tick_params(axis='x', labelsize=params['font_size_ticks'],
                   width=params['line_width_axes'], direction='in')
    ax.tick_params(axis='y', labelsize=params['font_size_ticks'],
                   width=params['line_width_axes'], direction='in')
    ax.margins(x=0.00, y=0.00)
    ax.grid(True, which='minor', linestyle=':', linewidth=params['line_width_grid_minor'])
    ax.grid(True, which='major', linestyle='--', linewidth=params['line_width_grid_major'])
    ax.set_xlabel("$L^{\\frac{1}{\\nu}}(\\beta - \\beta_c)$", fontsize=params['font_size_axis'], labelpad=params['label_pad'])
    ax.set_ylabel("$\\chi'\\cdot L^{\\frac{\\nu}{\\gamma}}$", fontsize=params['font_size_axis'], labelpad=params['label_pad'])
    ax.legend(loc='best', fontsize=params['font_size_legend'])

plt.tight_layout(pad=params['pad'])
plt.savefig('./figure/figure_7.pdf', format='pdf')
plt.close(fig)
