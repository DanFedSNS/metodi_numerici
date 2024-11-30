import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from scipy.optimize import curve_fit

# Array di valori L da considerare
L_array = np.linspace(10, 200, 20, dtype=int)
L_graph = np.linspace(110, 150, 5, dtype=int)
# Lista dei modelli
algos = ["ising2d_tri_cluster", "ising2d_tri_metro"]

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

# Definizione dei limiti per il FIT
x_lim_metro = (0, 100)  # Limiti per il metodo Metro
x_lim_cluster = (0, 30)  # Limiti per il metodo Cluster

fig, ax = plt.subplots(2, 1, figsize=(params['fig_width'], 1.5*params['fig_height']))


tau_metro = []
tau_cluster = []

# Funzione esponenziale per il fit
def exp_decay(x, A, tau):
    return A * np.exp(-x / tau)

# Funzione di potenza per il fit: tau = L^alpha
def power_law(L, alpha, c):
    return c * L**alpha

for algo_index, algo in enumerate(algos):
    tau_values = []  # Lista per memorizzare i tau per ogni L
    for i, L in enumerate(L_array):
        filepath = f'./data/corr/analysis_{algo}/L{L}_autocorr.dat'
        data = np.genfromtxt(filepath, delimiter=",", dtype=float, filling_values=np.nan)

        print(np.isnan(data).any(), np.isinf(data).any())

        autocorr_m = data
        x_indices = np.arange(len(autocorr_m))

        # Definizione del range per il fit
        if "metro" in algo:
            x_lim = x_lim_metro
        else:
            x_lim = x_lim_cluster

        # Filtra i dati per il fit
        mask_fit = (x_indices >= x_lim[0]) & (x_indices <= x_lim[1])
        x_fit = x_indices[mask_fit]
        y_fit = autocorr_m[mask_fit]

        label = f'L={L}'

        popt, _ = curve_fit(exp_decay, x_fit, y_fit, p0=[1, 10])
        A_fit, tau_fit = popt

        # Aggiunta del fit alla lista dei tau
        tau_values.append(tau_fit)

        if L in L_graph:
            ax[algo_index].plot(x_fit, exp_decay(x_fit, *popt), color=colori(algo_index), label=None, marker='none', linewidth=params['line_width_axes'], alpha=0.5, zorder = 0)
            ax[algo_index].plot(x_indices, autocorr_m, label=label, color=colori(algo_index), marker='o', linestyle='none', markerfacecolor='white', markeredgewidth = params['line_width_axes'], zorder = 2)
            
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
    
ax[0].set_xlabel("Numero di aggiornamenti", fontsize=params['font_size_axis'], labelpad=params['label_pad'])
ax[0].set_xlim(x_lim_metro)
ax[0].set_ylim([0, 1])
ax[0].set_ylabel("Autocorr di $m$", fontsize=params['font_size_axis'], labelpad=params['label_pad'])
#ax[0].set_yscale('log')
ax[0].legend(loc="upper right", fontsize=params['font_size_legend'])
#ax[0].set_xticks(ticks=[0.43, 0.44, 0.45])

ax[1].set_xlabel("Numero di aggiornamenti", fontsize=params['font_size_axis'], labelpad=params['label_pad'])
ax[1].set_xlim(x_lim_cluster)
ax[1].set_ylim([0, 1])
ax[1].set_ylabel("Autocorr di $m$", fontsize=params['font_size_axis'], labelpad=params['label_pad'])
#ax[1].set_yscale('log')
ax[1].legend(loc="upper right", fontsize=params['font_size_legend'])
#ax[1].set_xticks(ticks=[0.43, 0.44, 0.45])

plt.tight_layout(pad=params['pad'])
plt.savefig('./figure/figure_4.pdf', format='pdf')
plt.close(fig)

fig2, ax = plt.subplots(figsize=(params['fig_width'], 0.75*params['fig_height']))

# Plot dei tempi caratteristici al variare di L
ax.plot(L_array, tau_metro, label="Metro", color='blue', marker='o', linestyle='none', markerfacecolor='white', markeredgewidth = params['line_width_axes'], zorder = 2)
ax.plot(L_array, tau_cluster, label="Cluster", color='red', marker='o', linestyle='none', markerfacecolor='white', markeredgewidth = params['line_width_axes'], zorder = 2)

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

ax.plot(L_fit, power_law(L_fit, *popt_cluster), '--', color='red',
        label=f"Fit Cluster: $\\tau = L^{{{alpha_cluster:.2f}}}$", marker='none', linewidth=params['line_width_axes'], alpha=0.5, zorder = 0)


for spine in ax.spines.values():
    spine.set_linewidth(params['line_width_axes'])

ax.tick_params(axis='x', labelsize=params['font_size_ticks'], 
                width=params['line_width_axes'], direction='in')
ax.tick_params(axis='y', labelsize=params['font_size_ticks'], 
                width=params['line_width_axes'], direction='in')
ax.margins(x=0.00, y=0.00)
ax.grid(True, which='minor', linestyle=':', linewidth=params['line_width_grid_minor'])
ax.grid(True, which='major', linestyle='--', linewidth=params['line_width_grid_major'])

ax.legend(loc="upper left", fontsize=params['font_size_legend'])
ax.set_xlabel("L", fontsize=params['font_size_axis'], labelpad=params['label_pad'])
#ax.set_ylim([0, 100])
ax.set_ylabel("Tempo Caratteristico", fontsize=params['font_size_axis'], labelpad=params['label_pad'])
#ax[0].set_xticks(ticks=[0.43, 0.44, 0.45])

# Salva la figura
plt.tight_layout(pad=params['pad'])
plt.savefig('./figure/figure_5.pdf', format='pdf')
plt.close(fig)