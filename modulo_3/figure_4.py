import numpy as np
import matplotlib.pyplot as plt
import os
import glob

def load_params(filepath):
    parameters = {}
    with open(filepath, 'r') as file:
        for line in file:
            key, value = line.strip().split('=')
            if key in ['fig_width', 'fig_height', 'dpi']:
                parameters[key] = float(value)
            else:
                parameters[key] = float(value) if '.' in value else int(value)
    return parameters

# Carichiamo i parametri
os.chdir(os.path.dirname(os.path.abspath(__file__)))
plot_params = load_params('params.txt')

# Configurazione generale di Matplotlib
plt.rc('font', family='serif')
color_palette = plt.get_cmap('tab10')

# Creiamo una figura
fig, ax = plt.subplots(figsize=(plot_params['fig_width'], plot_params['fig_height']))

# Percorso ai file
data_directory = './analysis/fig456/'
data_filepaths = sorted(glob.glob(os.path.join(data_directory, '*')))  # Legge tutti i file nella directory

# Inizializziamo i vettori per i dati
beta_values = []
g_values = []
mean_x = []
mean_x_squared = []
mean_H = []
error_x = []
error_x_squared = []
error_H = []

# Lettura e salvataggio dei dati
for index, data_filepath in enumerate(data_filepaths):
    with open(data_filepath, 'r') as data_file:
        # Lettura iniziale del valore di beta
        for line in data_file:
            if line.startswith('Nt'):
                beta_line = line.split(',')[1]  # Secondo campo della riga
                beta_value = float(beta_line.split('=')[1].strip())
                beta_values.append(beta_value)
                g_line = line.split(',')[2]  # Terzo campo della riga
                g_value = float(g_line.split('=')[1].strip())
                g_values.append(g_value)
                break
        
        # Legge tutte le righe successive come dati numerici
        numerical_data = np.loadtxt(data_file)

    # Salviamo i dati in vettori
    mean_x.append(numerical_data[0, 0])
    error_x.append(numerical_data[0, 1])
    mean_x_squared.append(numerical_data[1, 0])
    error_x_squared.append(numerical_data[1, 1])
    mean_H.append(numerical_data[70, 0])
    error_H.append(numerical_data[70, 1])

# Convertiamo i vettori in array numpy per il plotting
beta_values = np.array(beta_values)
g_values = np.array(g_values)
mean_x = np.array(mean_x)
error_x = np.array(error_x)
mean_x_squared = np.array(mean_x_squared)
error_x_squared = np.array(error_x_squared)
entropy_H = np.array(mean_H)
error_H = np.array(error_H)

# Eseguiamo il plot
ax.errorbar(g_values, mean_x, yerr=error_x, fmt='o', color=color_palette(0),
            label="Mean X", markerfacecolor='white', markeredgewidth=plot_params['line_width_axes'], zorder=2)
ax.errorbar(g_values, mean_x_squared, yerr=error_x_squared, fmt='o', color=color_palette(1),
            label="Mean X^2", markerfacecolor='white', markeredgewidth=plot_params['line_width_axes'], zorder=2)
ax.errorbar(g_values, mean_H, yerr=error_H, fmt='o', color=color_palette(2),
            label="Entropy H", markerfacecolor='white', markeredgewidth=plot_params['line_width_axes'], zorder=2)

# Personalizzazione dei subplot
for spine in ax.spines.values():
    spine.set_linewidth(plot_params['line_width_axes'])

ax.tick_params(axis='x', labelsize=plot_params['font_size_ticks'], 
               width=plot_params['line_width_axes'], direction='in')
ax.tick_params(axis='y', labelsize=plot_params['font_size_ticks'], 
               width=plot_params['line_width_axes'], direction='in')
ax.margins(x=0.00, y=0.00)
ax.grid(True, which='minor', linestyle=':', linewidth=plot_params['line_width_grid_minor'])
ax.grid(True, which='major', linestyle='--', linewidth=plot_params['line_width_grid_major'])

# Etichette degli assi
ax.set_xlabel("$g$", fontsize=plot_params['font_size_axis'], labelpad=plot_params['label_pad'])
ax.set_ylabel("$C$", fontsize=plot_params['font_size_axis'], labelpad=plot_params['label_pad'])
ax.set_xscale('log')
ax.legend(loc='best', fontsize=plot_params['font_size_legend'])

# Salvataggio della figura
plt.tight_layout(pad=plot_params['pad'])
plt.savefig('./figure/figure_4.pdf', format='pdf')
plt.close(fig)