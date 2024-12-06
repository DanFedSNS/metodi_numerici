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

# Creiamo una figura
fig, ax = plt.subplots(figsize=(params['fig_width'], params['fig_height']))

# Percorso ai file
data_dir = './analysis/fig1/'
filepaths = sorted(glob.glob(os.path.join(data_dir, '*')))  # Legge tutti i file nella directory

# Inizializziamo i vettori per i dati
betas = []
curve1 = []
curve2 = []
curve3 = []
errors1 = []
errors2 = []
errors3 = []

# Lettura e salvataggio dei dati
for i, filepath in enumerate(filepaths):
    with open(filepath, 'r') as file:
        # Lettura iniziale del valore di beta
        for line in file:
            if line.startswith('Nt'):
                simbeta_str = line.split(',')[1]  # Secondo campo della riga
                beta = float(simbeta_str.split('=')[1].strip())
                betas.append(beta)
                break
        
        # Legge tutte le righe successive come dati numerici
        data = np.loadtxt(file)

    # Salviamo i dati in vettori
    curve1.append(data[0, 0])
    errors1.append(data[0, 1])
    curve2.append(data[1, 0])
    errors2.append(data[1, 1])
    curve3.append(data[3, 0])
    errors3.append(data[3, 1])

# Convertiamo i vettori in array numpy per il plotting
betas = np.array(betas)
curve1 = np.array(curve1)
errors1 = np.array(errors1)
curve2 = np.array(curve2)
errors2 = np.array(errors2)
curve3 = np.array(curve3)
errors3 = np.array(errors3)

beta_dense = np.linspace(betas.min(), betas.max(), 1000)
y_anal = 0.5*(1 + np.e**(-beta_dense))/(1-np.e**(-beta_dense))

# Eseguiamo il plot
# ax.errorbar(betas, curve1, yerr=errors1, fmt='o', color=colors(0),
#             label="Curve 1", markerfacecolor='white', markeredgewidth=params['line_width_axes'], zorder=2)
ax.errorbar(betas, curve2, yerr=errors2, fmt='o', color='black', markerfacecolor='white', markeredgewidth=params['line_width_axes'], zorder=2)
# ax.errorbar(betas, curve3, yerr=errors3, fmt='o', color=colors(2),
#             label="Curve 3", markerfacecolor='white', markeredgewidth=params['line_width_axes'], zorder=2)
ax.plot(beta_dense, y_anal, color='black', markerfacecolor='white', markeredgewidth=params['line_width_axes'], zorder=2, alpha=0.5)

# Personalizzazione dei subplot
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
ax.set_xlabel("$\\beta$", fontsize=params['font_size_axis'], labelpad=params['label_pad'])
ax.set_ylabel("$\\leftangle H \\rightangle$", fontsize=params['font_size_axis'], labelpad=params['label_pad'])
ax.set_xscale('log')
ax.set_yscale('log')

# Salvataggio della figura
plt.tight_layout(pad=params['pad'])
plt.savefig('./figure/figure_1.pdf', format='pdf')
plt.close(fig)
