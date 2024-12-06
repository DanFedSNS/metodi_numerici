import numpy as np
import matplotlib.pyplot as plt
import os

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

def extract_data(filepath):
    with open(filepath, 'r') as file:
        header = file.readline().strip()
        beta = float(header.split('=')[1]) 
        
        data = np.genfromtxt(file, max_rows=3, dtype=float)
    return beta, data

# Carica i parametri globali.
params = load_params('params.txt')

# Configurazione del plot.
plt.rc('font', family='serif')
colors = plt.get_cmap('tab10')

# Directory contenente i file.
data_dir = './fig12/'  # Cambia con il percorso corretto.

# Lista dei file nella directory.
files = [os.path.join(data_dir, f) for f in os.listdir(data_dir) if f.endswith('.txt')]

# Creazione del plot.
fig, ax = plt.subplots(figsize=(params['fig_width'], params['fig_height']))

for filepath in files:
    # Estrazione dei dati.
    beta, data = extract_data(filepath)
    
    # Indici riga e valori corrispondenti.
    x_values = [1, 2, 3]
    y_values = data[:, 0]  # Supponendo che i valori siano nella prima colonna.
    y_errors = data[:, 1]  # Supponendo che gli errori siano nella seconda colonna.
    
    # Aggiunge i dati al grafico.
    ax.errorbar(x_values, y_values, yerr=y_errors, label=f"$\\beta={beta:.2f}$",
                marker='o', linestyle='-', capsize=3)

# Configurazione degli assi.
ax.set_xlabel("Indice riga", fontsize=params['font_size_axis'], labelpad=params['label_pad'])
ax.set_ylabel("Valore", fontsize=params['font_size_axis'], labelpad=params['label_pad'])
ax.legend(loc='best', fontsize=params['font_size_legend'])

# Miglioramenti estetici.
for spine in ax.spines.values():
    spine.set_linewidth(params['line_width_axes'])
ax.tick_params(axis='x', labelsize=params['font_size_ticks'], 
               width=params['line_width_axes'], direction='in')
ax.tick_params(axis='y', labelsize=params['font_size_ticks'], 
               width=params['line_width_axes'], direction='in')
ax.margins(x=0.00, y=0.00)
ax.grid(True, which='minor', linestyle=':', linewidth=params['line_width_grid_minor'])
ax.grid(True, which='major', linestyle='--', linewidth=params['line_width_grid_major'])

plt.tight_layout(pad=params['pad'])

# Salva la figura.
output_path = './figure/figure_1.pdf'
plt.savefig(output_path, format='pdf')
plt.close(fig)

print(f"Figura salvata in {output_path}")
