import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from scipy.optimize import curve_fit
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

os.chdir(os.path.dirname(os.path.abspath(__file__)))
params = load_params('params.txt')

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

colori = plt.get_cmap('tab10')

data_directory = './analysis/fig9/'
data_filepaths = sorted(glob.glob(os.path.join(data_directory, '*')))  # Legge tutti i file nella directory

fig2, ax2 = plt.subplots(figsize=(params['fig_width'], 0.75*params['fig_height']))

for index, data_filepath in enumerate(data_filepaths):
    with open(data_filepath, 'r') as data_file:
        for line in data_file:
            if line.startswith('Nt'):
                Nt_line = line.split(',')[0] 
                Nt_value = float(Nt_line.split('=')[1].strip())
                simbeta_str = line.split(',')[1]  # Secondo campo della riga
                beta = float(simbeta_str.split('=')[1].strip())
                break
        
        numerical_data = np.loadtxt(data_file)
    
    corr_x_2 = numerical_data[0,0]
    
    ax2.plot(beta, corr_x_2,  's', color='blue', linestyle='none', markerfacecolor='white', markeredgewidth = params['line_width_axes'], zorder = 2)

beta_dense = np.linspace(1, 6, 200)
tau = 1/5
#anal = 1/2 * np.exp(-tau) + np.exp(-beta_dense)*(np.cosh(tau))
anal = 1/2 * (np.exp(tau - beta_dense) + np.exp(-tau)) / (1 - np.exp(-beta_dense))
val_esatto = 1/2 * np.exp(-tau)

ax2.plot(beta_dense, anal, color='blue', marker='none', linewidth=params['line_width_axes'], alpha=0.5, zorder = 0)
ax2.plot(beta_dense, [val_esatto for _ in beta_dense], color='red', marker='none', linestyle = "--", linewidth=params['line_width_axes'], alpha=0.5, zorder = 0)
for spine in ax2.spines.values():
    spine.set_linewidth(params['line_width_axes'])

ax2.tick_params(axis='x', labelsize=params['font_size_ticks'], 
                width=params['line_width_axes'], direction='in')
ax2.tick_params(axis='y', color="black", labelsize=params['font_size_ticks'], 
                width=params['line_width_axes'], direction='in')

ax2.margins(x=0.00, y=0.00)
ax2.grid(True, which='minor', linestyle=':', linewidth=params['line_width_grid_minor'])
ax2.grid(True, which='major', linestyle='--', linewidth=params['line_width_grid_major'])

ax2.set_xlabel("$\\beta$", fontsize=params['font_size_axis'], labelpad=params['label_pad'])
ax2.set_ylabel("$\\langle x(\\tau_0)x(0)\\rangle$", fontsize=params['font_size_axis'], labelpad=params['label_pad'], color = "black")
#ax2.set_xscale('log')
ax2.set_ylim([0.4, 1.1])
# Salva la figura
plt.tight_layout(pad=params['pad'])
plt.savefig('./figure/figure_9.pdf', format='pdf')
plt.close(fig2)