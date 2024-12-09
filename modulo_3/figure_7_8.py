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

fig, ax = plt.subplots(figsize=(params['fig_width'], params['fig_height']))

def exp_decay(x, A, tau):
    return A * np.exp(-x / tau)

def power_law(L, alpha, c):
    return c * L**alpha 

data_directory = './analysis/fig78/'
data_filepaths = sorted(glob.glob(os.path.join(data_directory, '*')))  # Legge tutti i file nella directory

g_values = []
corr_x = []

tau_values = []  
Nt_values = [] 
color_index = 0 

num = 6
x_lim = [0, 2000]
y_lim = 2e-1

for index, data_filepath in enumerate(data_filepaths):
    with open(data_filepath, 'r') as data_file:
        for line in data_file:
            if line.startswith('Nt'):
                Nt_line = line.split(',')[0] 
                Nt_value = float(Nt_line.split('=')[1].strip())
                Nt_values.append(Nt_value)
                break
        
        numerical_data = np.loadtxt(data_file)
    
    autocorr_m = numerical_data
    x_indices = np.arange(len(autocorr_m))

    mask_fit = (autocorr_m >= y_lim) & (x_indices <= x_lim[1]) & (x_indices >= x_lim[0])

    x_fit = x_indices[mask_fit]
    y_fit = autocorr_m[mask_fit]

    label = f'Nt={Nt_value}'

    popt, _ = curve_fit(exp_decay, x_fit, y_fit, p0=[1, 10])
    A_fit, tau_fit = popt

    tau_values.append(tau_fit)
    x_dense = np.linspace(min(x_fit), max(x_fit), 1000)

    x_plot = x_fit
    y_plot = y_fit

    # if L in L_graph:
    ax.plot(x_plot[::num], y_plot[::num], label=label, color=colori(color_index), marker='o', linestyle='none', markerfacecolor='white', markeredgewidth = params['line_width_axes'], zorder = 2)
    # ax.set_xlim(x_plot.min(), x_plot.max())
    # ax.set_ylim(y_plot.min(), y_plot.max()) 
    # ax.plot(x_dense, exp_decay(x_dense, *popt), color=colori(color_index), label=None, marker='none', linewidth=params['line_width_axes'], alpha=0.5, zorder = 0)
    color_index += 1
        

for spine in ax.spines.values():
    spine.set_linewidth(params['line_width_axes'])
ax.tick_params(axis='x', labelsize=params['font_size_ticks'], 
                width=params['line_width_axes'], direction='in')
ax.tick_params(axis='y', labelsize=params['font_size_ticks'], 
                width=params['line_width_axes'], direction='in')
ax.margins(x=0.00, y=0.00)
ax.grid(True, which='minor', linestyle=':', linewidth=params['line_width_grid_minor'])
ax.grid(True, which='major', linestyle='--', linewidth=params['line_width_grid_major'])
ax.set_ylabel("Autocorr", fontsize=params['font_size_axis'], labelpad=params['label_pad'])
#ax.legend(loc="upper right", fontsize=params['font_size_legend'])
ax.set_xlabel("Numero di aggiornamenti", fontsize=params['font_size_axis'], labelpad=params['label_pad'])
ax.set_yscale("log")

plt.tight_layout(pad=params['pad'])
#plt.savefig('./figure/figure_7.pdf', format='pdf')
plt.close(fig)

styles = ['o', 's']

fig2, ax2 = plt.subplots(figsize=(params['fig_width'], 0.75*params['fig_height']))

tau_values =  np.array(tau_values)*10
# Plot dei tempi caratteristici al variare di L
ax2.plot(Nt_values, tau_values, 's', color='blue', linestyle='none', markerfacecolor='white', markeredgewidth = params['line_width_axes'], zorder = 2)

popt_, _ = curve_fit(power_law, Nt_values, tau_values, p0=[2.1, 1])

Nt_fit = np.linspace(1, 800, 1000)

# ax2.plot(Nt_fit, power_law(Nt_fit, *popt_), color='blue', marker='none', linewidth=params['line_width_axes'], alpha=0.5, zorder = 0)
print(popt_)

for spine in ax2.spines.values():
    spine.set_linewidth(params['line_width_axes'])

ax2.tick_params(axis='x', labelsize=params['font_size_ticks'], 
                width=params['line_width_axes'], direction='in')
ax2.tick_params(axis='y', color="black", labelsize=params['font_size_ticks'], 
                width=params['line_width_axes'], direction='in')

ax2.margins(x=0.00, y=0.00)
ax2.grid(True, which='minor', linestyle=':', linewidth=params['line_width_grid_minor'])
ax2.grid(True, which='major', linestyle='--', linewidth=params['line_width_grid_major'])

#ax2.legend(loc="upper left", fontsize=params['font_size_legend'])
ax2.set_xlabel("$N_T$", fontsize=params['font_size_axis'], labelpad=params['label_pad'])
ax2.set_xlim([9, 600])
ax2.set_ylim([100, 6e5])
ax2.set_ylabel("$\\tau_{corr}$", fontsize=params['font_size_axis'], labelpad=params['label_pad'], color = "black")
#ax2.set_xticks(ticks=[0.43, 0.44, 0.45])
ax2.set_yscale("log")
ax2.set_xscale("log")

# Salva la figura
plt.tight_layout(pad=params['pad'])
plt.savefig('./figure/figure_8.pdf', format='pdf')
plt.close(fig)