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

fig, ax = plt.subplots(1, 1, figsize=(params['fig_width'], 5/6*params['fig_height']))

for algo_index, algo in enumerate(algos):
    for i, L in enumerate(L_array):
        filepath = f'./data/corr/analysis_{algo}/L{L}.dat'
        data = np.genfromtxt(filepath, delimiter=" ", dtype=float, filling_values=np.nan)
         
        sigma_m = data[6] 

        # Plotta la suscettività per ogni modello al variare di L
        ax.plot(L, sigma_m, styles[algo_index], color=colors(algo_index), label=f'{algo.split("_")[2].capitalize()}' if i == 0 else "", markerfacecolor='white', markeredgewidth = params['line_width_axes'], zorder = 2)
    
for spine in ax.spines.values():
        spine.set_linewidth(params['line_width_axes'])

ax.tick_params(axis='x', labelsize=params['font_size_ticks'], 
                width=params['line_width_axes'], direction='in')
ax.tick_params(axis='y', labelsize=params['font_size_ticks'], 
                width=params['line_width_axes'], direction='in')
ax.margins(x=0.00, y=0.00)
ax.grid(True, which='minor', linestyle=':', linewidth=params['line_width_grid_minor'])
ax.grid(True, which='major', linestyle='--', linewidth=params['line_width_grid_major'])

ax.set_xlabel("L", fontsize=params['font_size_axis'], labelpad=params['label_pad'])
ax.set_ylabel("$\\sigma_{\\langle m \\rangle}$", fontsize=params['font_size_axis'], labelpad=params['label_pad'])
ax.legend(loc='best', fontsize=params['font_size_legend'])
ax.set_xticks(ticks=[10, 50, 80, 100])

plt.tight_layout(pad=params['pad'])
plt.savefig('./figure/figure_3.pdf', format='pdf')
plt.close(fig)