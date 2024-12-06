import numpy as np
import matplotlib.pyplot as plt
import os
import glob

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

# Carichiamo i parametri
os.chdir(os.path.dirname(os.path.abspath(__file__)))
params = load_params('params.txt')

# Configurazione generale di Matplotlib
plt.rc('font', family='serif')
colors = plt.get_cmap('tab10')

# Cartella dei dati
data_folder = './analysis/fig5/'
file_pattern = os.path.join(data_folder, 'wf_Nt*_simbeta*_g*.dat')
files = glob.glob(file_pattern)

# Creiamo una figura
fig, ax = plt.subplots(figsize=(params['fig_width'], 0.75*params['fig_height']))

for idx, file_path in enumerate(sorted(files)):
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('Nt'):
                g_str = line.split(',')[1]
                g = float(g_str.split('=')[1].strip())
                break
    
    # Caricare i dati
    data = np.loadtxt(file_path, skiprows=1)
    positions = data[:, 0]
    occurrences = data[:, 1]
    occurrences /= sum(occurrences)
    
    # Plottare la curva
    ax.plot(positions, occurrences, label=f"$g = {g:.3f}$", 
            color=colors(idx % 10), 
            markerfacecolor='white', 
            markeredgewidth=params['line_width_axes'], 
            zorder=2)


for spine in ax.spines.values():
    spine.set_linewidth(params['line_width_axes'])

ax.tick_params(axis='x', labelsize=params['font_size_ticks'], 
               width=params['line_width_axes'], direction='in')
ax.tick_params(axis='y', labelsize=params['font_size_ticks'], 
               width=params['line_width_axes'], direction='in')
ax.margins(x=0.00, y=0.00)
ax.grid(True, which='minor', linestyle=':', linewidth=params['line_width_grid_minor'])
ax.grid(True, which='major', linestyle='--', linewidth=params['line_width_grid_major'])

# Etichette degli assi
ax.set_xlabel("$x$", fontsize=params['font_size_axis'], labelpad=params['label_pad'])
ax.set_ylabel("$\\psi^2$", fontsize=params['font_size_axis'], labelpad=params['label_pad'])

# Aggiungere la legenda
ax.legend(fontsize=params['font_size_legend'])

# Salvataggio della figura
plt.tight_layout(pad=params['pad'])
plt.savefig('./figure/figure_5.pdf', format='pdf')
plt.close(fig)
