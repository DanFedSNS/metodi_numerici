import numpy as np
import matplotlib.pyplot as plt
import os
import sys

# Parametri
L_array = [20]
beta_unificato = [0.21, 0.36, 0.56, 0.79]  # Unico vettore di beta per tutti i modelli
tau = int(5e4)
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
fig, ax = plt.subplots(1, 4, figsize=(2*params['fig_width']+0.25, 0.5*params['fig_height']))

colori = plt.get_cmap('tab10')

# Definire i limiti dell'asse x
x_limits_magnetization = [0, 1]

# Numero di intervalli per il binning
n_intervalli = 401
bin_edges = np.linspace(-1, 1, n_intervalli + 1)

# Loop sui valori di beta
for j, beta in enumerate(beta_unificato):
    # Loop sui modelli
    for k, modello in enumerate(modelli):
        # Loop sui valori di L
        for i in range(len(L_array)):
            L = L_array[i]

            # Costruzione del percorso del file
            filepath = f'./data/{modello}/L{L}_beta{(beta):.5f}.dat'

            # Caricamento dei dati
            data = np.loadtxt(filepath, skiprows=1, delimiter=",")
            magnetization = data[:, 0]
            magnetization = magnetization[tau:]
            
            # Binning classico: calcolo della frequenza per ogni bin
            magnetization_frequencies, _ = np.histogram(magnetization, bins=bin_edges)
            magnetization_frequencies = magnetization_frequencies / len(magnetization)
            # Calcolo del centro di ogni bin per il plotting
            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

            # Plot della distribuzione di magnetizzazione
            label = f'{modello.split("_")[1]}'
            ax[j].plot(
                bin_centers, magnetization_frequencies,
                color=colori(k), label=label,
                marker='none', linewidth=params['line_width_axes']
            )

    # Personalizzazione del grafico
    for spine in ax[j].spines.values():
        spine.set_linewidth(params['line_width_axes'])

    ax[j].set_title('$\\beta = ' + str(beta) + '$', fontsize=params['font_size_axis'])
    ax[j].set_xlabel('m', fontsize=params['font_size_axis'], labelpad=params['label_pad'])
    ax[j].set_ylabel('P(m)', fontsize=params['font_size_axis'], labelpad=params['label_pad'])
    ax[j].set_xlim(x_limits_magnetization)

    handles, labels = ax[j].get_legend_handles_labels()
    ax[j].legend(handles, labels, loc='best', fontsize=params['font_size_legend'])
    ax[j].grid(True, which='minor', linestyle=':', linewidth=params['line_width_grid_minor'])
    ax[j].grid(True, which='major', linestyle='--', linewidth=params['line_width_grid_major'])

    ax[j].tick_params(axis='x', labelsize=params['font_size_ticks'], 
               width=params['line_width_axes'], direction='in')
    ax[j].tick_params(axis='y', labelsize=params['font_size_ticks'], 
               width=params['line_width_axes'], direction='in')
    ax[j].set_xticks(ticks=[0, 0.5, 1])
    ax[j].set_yticks(ticks=[0, 0.05, 0.1])
    ax[j].set_ylim([0, 0.1])
    ax[j].margins(x=0.00, y=0.00)

# Migliora la spaziatura tra i subplot
fig.tight_layout(pad=params['pad'])
plt.savefig('./figure/figure_1.pdf', format='pdf', dpi=int(params['dpi']))
plt.close(fig)
