import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from scipy.optimize import curve_fit

# Array di valori L da considerare
L_array = np.linspace(10, 200, 20, dtype=int)
L_graph = np.linspace(70, 150, 5, dtype=int)
# Lista dei modelli
algos = ["ising2d_tri_metro", "ising2d_tri_cluster"]

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


num = [300, 2]
y_lim_metro = 1e-1
y_lim_cluster = 5e-3
x_lim_metro = [0, 15e3]
x_lim_cluster = [5, 50]

fig, ax = plt.subplots(2, 1, figsize=(params['fig_width'], 1.5*params['fig_height']))

tau_metro = []
tau_cluster = []

def exp_decay(x, A, tau):
    return A * np.exp(-x / tau)

def power_law(L, alpha, c):
    return c * L**alpha

for algo_index, algo in enumerate(algos):
    tau_values = []  # Lista per memorizzare i tau per ogni L
    color_index = 0 
    for i, L in enumerate(L_array):
        filepath = f'./data/corr/analysis_{algo}/L{L}_autocorr.dat'
        data = np.genfromtxt(filepath, delimiter=",", dtype=float, filling_values=np.nan)

        autocorr_m = data
        x_indices = np.arange(len(autocorr_m))

        y_lim = y_lim_metro if "metro" in algo else y_lim_cluster
        x_lim = x_lim_metro if "metro" in algo else x_lim_cluster
        mask_fit = (autocorr_m >= y_lim) & (x_indices <= x_lim[1]) & (x_indices >= x_lim[0])

        x_fit = x_indices[mask_fit]
        y_fit = autocorr_m[mask_fit]

        label = f'L={L}'

        popt, _ = curve_fit(exp_decay, x_fit, y_fit, p0=[1, 10])
        A_fit, tau_fit = popt

        tau_values.append(tau_fit)
        x_dense = np.linspace(min(x_fit), max(x_fit)+1000, 1000)

        x_plot = x_fit[::num[algo_index]]
        y_plot = y_fit[::num[algo_index]]

        if L in L_graph:
            ax[algo_index].plot(x_plot, y_plot, label=label, color=colori(color_index), marker='s', linestyle='none', markerfacecolor='white', markeredgewidth = params['line_width_axes'], zorder = 2)
            ax[algo_index].set_xlim(x_plot.min(), x_plot.max())
            ax[algo_index].set_ylim(y_plot.min(), y_plot.max()) 
            ax[algo_index].plot(x_dense, exp_decay(x_dense, *popt), color=colori(color_index), label=None, marker='none', linewidth=params['line_width_axes'], alpha=0.5, zorder = 0)
            color_index += 1
            
    if "metro" in algo:
        tau_metro = tau_values
    else:
        tau_cluster = tau_values

for ax_ in ax.flat:
    for spine in ax_.spines.values():
        spine.set_linewidth(params['line_width_axes'])
    ax_.tick_params(axis='x', labelsize=params['font_size_ticks'], 
                    width=params['line_width_axes'], direction='in')
    ax_.tick_params(axis='y', labelsize=params['font_size_ticks'], 
                    width=params['line_width_axes'], direction='in')
    ax_.margins(x=0.00, y=0.00)
    ax_.grid(True, which='minor', linestyle=':', linewidth=params['line_width_grid_minor'])
    ax_.grid(True, which='major', linestyle='--', linewidth=params['line_width_grid_major'])
    ax_.set_ylabel("Autocorr di $m$", fontsize=params['font_size_axis'], labelpad=params['label_pad'])
    ax_.legend(loc="upper right", fontsize=params['font_size_legend'])
    ax_.set_xlabel("Numero di aggiornamenti", fontsize=params['font_size_axis'], labelpad=params['label_pad'])
    ax_.set_yscale("log")


plt.tight_layout(pad=params['pad'])
plt.savefig('./figure/figure_4.pdf', format='pdf')
plt.close(fig)

styles = ['s', 's']

fig2, ax = plt.subplots(figsize=(params['fig_width'], 0.75*params['fig_height']))

# Creare un secondo asse Y
ax2 = ax.twinx()

# Plot dei tempi caratteristici al variare di L
ax.plot(L_array, tau_metro, styles[0], label="Metro", color='blue', linestyle='none', markerfacecolor='white', markeredgewidth = params['line_width_axes'], zorder = 2)
ax2.plot(L_array, tau_cluster, styles[1], label="Cluster", color='red', linestyle='none', markerfacecolor='white', markeredgewidth = params['line_width_axes'], zorder = 2)

# Fit per i dati Metro
popt_metro, _ = curve_fit(power_law, L_array, tau_metro, p0=[2.1, 1])
alpha_metro, c_metro = popt_metro

# Fit per i dati Cluster
popt_cluster, _ = curve_fit(power_law, L_array, tau_cluster, p0=[0.2, 1])
alpha_cluster, c_cluster = popt_cluster

print(f"Fit Metro: alpha = {alpha_metro:.4f}, c = {c_metro:.4f}")
print(f"Fit Cluster: alpha = {alpha_cluster:.4f}, c = {c_cluster:.4f}")

# Plot del fit per Metro
L_fit = np.linspace(min(L_array), max(L_array), 100)

ax.plot(L_fit, power_law(L_fit, *popt_metro), color='blue',
        label=f"Fit Metro: $\\tau = L^{{{alpha_metro:.2f}}}$", marker='none', linewidth=params['line_width_axes'], alpha=0.5, zorder = 0)

ax2.plot(L_fit, power_law(L_fit, *popt_cluster), '--', color='red',
        label=f"Fit Cluster: $\\tau = L^{{{alpha_cluster:.2f}}}$", marker='none', linewidth=params['line_width_axes'], alpha=0.5, zorder = 0)


for spine in ax.spines.values():
    spine.set_linewidth(params['line_width_axes'])
for spine in ax2.spines.values():
    spine.set_linewidth(params['line_width_axes'])

# Assegnare colori agli assi
ax.spines['left'].set_color("blue") 
ax2.spines['right'].set_color("red")  

ax.tick_params(axis='x', labelsize=params['font_size_ticks'], 
                width=params['line_width_axes'], direction='in')
ax.tick_params(axis='y', color="blue", labelsize=params['font_size_ticks'], 
                width=params['line_width_axes'], direction='in')
ax2.tick_params(axis='y', color="red", labelsize=params['font_size_ticks'], 
                width=params['line_width_axes'], direction='in')

ax.margins(x=0.00, y=0.00)
ax.grid(True, which='minor', linestyle=':', linewidth=params['line_width_grid_minor'])
ax.grid(True, which='major', linestyle='--', linewidth=params['line_width_grid_major'])

ax.legend(loc="upper left", fontsize=params['font_size_legend'])
ax2.legend(loc="lower right", fontsize=params['font_size_legend'])
ax.set_xlabel("L", fontsize=params['font_size_axis'], labelpad=params['label_pad'])
#ax.set_ylim([0, 100])
ax.set_ylabel("$\\tau$", fontsize=params['font_size_axis'], labelpad=params['label_pad'], color = "blue")
ax2.set_ylabel("$\\tau$", fontsize=params['font_size_axis'], labelpad=params['label_pad'], color = "red")
#ax[0].set_xticks(ticks=[0.43, 0.44, 0.45])

# Salva la figura
plt.tight_layout(pad=params['pad'])
plt.savefig('./figure/figure_5.pdf', format='pdf')
plt.close(fig)