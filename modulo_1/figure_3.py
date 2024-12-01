import numpy as np
import matplotlib.pyplot as plt
import os
import sys

# Array di valori L da considerare
L_array = np.linspace(10, 200, 20, dtype=int)

# Lista dei modelli
algos = ["ising2d_tri_metro", "ising2d_tri_cluster"]

# Stili per i modelli
styles = ['o', 's']

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

colors = plt.get_cmap('tab10')

fig, ax1 = plt.subplots(figsize=(params['fig_width'], 5 / 6 * params['fig_height']))

# Creare un secondo asse Y
ax2 = ax1.twinx()

for algo_index, algo in enumerate(algos):
    for i, L in enumerate(L_array):
        filepath = f'./data/corr/analysis_{algo}/L{L}.dat'
        data = np.genfromtxt(filepath, delimiter=" ", dtype=float, filling_values=np.nan)
        
        sigma_m = data[6]*1e3
        
        # Plotta i dati sul primo o secondo asse in base al modello
        if algo_index == 0:
            ax1.plot(L, sigma_m, styles[algo_index], color="blue", 
                     label=f'{algo.split("_")[2].capitalize()}' if i == 0 else "",
                     markerfacecolor='white', markeredgewidth=params['line_width_axes'], zorder=2)
        else:
            ax2.plot(L, sigma_m, styles[algo_index], color="red", 
                     label=f'{algo.split("_")[2].capitalize()}' if i == 0 else "",
                     markerfacecolor='white', markeredgewidth=params['line_width_axes'], zorder=2)

# Personalizzazione degli assi
for spine in ax1.spines.values():
    spine.set_linewidth(params['line_width_axes'])
for spine in ax2.spines.values():
    spine.set_linewidth(params['line_width_axes'])

# Assegnare colori agli assi
ax1.spines['left'].set_color("blue") 
ax2.spines['right'].set_color("red")  

# Tick e griglia
ax1.tick_params(axis='x', labelsize=params['font_size_ticks'], 
                width=params['line_width_axes'], direction='in')
ax1.tick_params(axis='y', labelsize=params['font_size_ticks'], 
                width=params['line_width_axes'], direction='in', colors="blue")
ax2.tick_params(axis='y', labelsize=params['font_size_ticks'], 
                width=params['line_width_axes'], direction='in', colors="red")

ax1.grid(True, which='minor', linestyle=':', linewidth=params['line_width_grid_minor'])
ax1.grid(True, which='major', linestyle='--', linewidth=params['line_width_grid_major'])

# Etichette e legenda
ax1.set_xlabel("L", fontsize=params['font_size_axis'], labelpad=params['label_pad'])
ax1.set_ylabel("$\\sigma_{\\langle m \\rangle} \\cdot 10^3$", fontsize=params['font_size_axis'], 
               labelpad=params['label_pad'], color="blue")
ax2.set_ylabel("$\\sigma_{\\langle m \\rangle} \\cdot 10^3$", fontsize=params['font_size_axis'], 
               labelpad=params['label_pad'], color="red")

ax1.legend(loc='upper left', fontsize=params['font_size_legend'])
ax2.legend(loc='upper right', fontsize=params['font_size_legend'])

ax1.set_xticks(ticks=[10, 100, 200])

plt.tight_layout(pad=params['pad'])
plt.savefig('./figure/figure_3.pdf', format='pdf')
plt.close(fig)
