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

def sol_anal(x):
    return 0.5*(1 + np.e**(-x))/(1-np.e**(-x))

beta_dense = np.linspace(betas.min(), betas.max(), 1000)
y_anal = sol_anal(beta_dense)

# Creiamo una figura
fig, ax = plt.subplots(2, 1, figsize=(params['fig_width'], 1*params['fig_height']), gridspec_kw={'height_ratios': [3, 1]}, sharex=True)

ax[0].plot(beta_dense, y_anal, color='black', markerfacecolor='white', markeredgewidth=params['line_width_axes'], zorder=2, alpha=0.5)
ax[0].plot(betas, curve2, 's', color='black', markerfacecolor='white', markeredgewidth=params['line_width_axes'], zorder=2)

res = (curve2 - sol_anal(betas))/errors2
ax[1].plot(betas, res, 's', color='black', markerfacecolor='white')


ax[1].axhline(0, color='black', linewidth=0.8, linestyle='--')
ax[1].set_ylabel('Residui', fontsize=params['font_size_axis'], labelpad=params['label_pad'])
ax[1].grid(True)

for ax_ in ax:
    ax_.grid(True)
    ax_.tick_params(axis='x', labelsize=params.get('font_size_ticks', 10), direction='in')
    ax_.tick_params(axis='y', labelsize=params.get('font_size_ticks', 10), direction='in')
    for spine in ax_.spines.values():
        spine.set_linewidth(params.get('line_width_axes', 1))
    ax_.margins(x=0.00, y=0.00)
    ax_.grid(True, which='minor', linestyle=':', linewidth=params['line_width_grid_minor'])
    ax_.grid(True, which='major', linestyle='--', linewidth=params['line_width_grid_major'])

# Etichette degli assi
ax[0].set_xlabel("$\\beta$", fontsize=params['font_size_axis'], labelpad=params['label_pad'])
ax[0].set_ylabel("$\\leftangle H \\rightangle$", fontsize=params['font_size_axis'], labelpad=params['label_pad'])
ax[0].set_xscale('log')
ax[0].set_yscale('log')
ax[1].set_xscale('log')

# Salvataggio della figura
plt.tight_layout(pad=params['pad'])
plt.savefig('./figure/figure_1.pdf', format='pdf')
plt.close(fig)
