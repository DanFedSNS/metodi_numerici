import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os
import sys

# Definizione della parabola per il fit
def parabola(x, A, B, C):
    return A * (x - B)**2 + C

# Definizione della legge di potenza con offset
def power_law(x, A, n, C):
    return A * x**n + C

# Parametri
L_array = np.linspace(70, 120, 6, dtype=int)
L_dense = np.linspace(min(L_array), max(L_array), 200)  # 200 punti per una maggiore densità
modelli = ["ising2d_sq_cluster", "ising2d_tri_cluster", "ising2d_hex_cluster"]

# Creazione dei dizionari per salvare i parametri del fit
fit_params = {model: {'A': [], 'B': [], 'C': []} for model in modelli}

# Definizione dei valori iniziali del parametro p0 per ogni modello
p0_values = {
    "ising2d_sq_cluster": [100, 0.4375, 150],  # Valori per il modello ising2d_sq_cluster
    "ising2d_tri_cluster": [100, 0.4375, 150],  # Valori per il modello ising2d_tri_cluster
    "ising2d_hex_cluster": [100, 0.4375, 150]   # Valori per il modello ising2d_hex_cluster
}

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

colors = plt.get_cmap('tab10') 

# Loop sui modelli e sui diversi valori di L
for model in modelli:
    for i, L in enumerate(L_array):
        filepath = f'./analysis_{model}/L{L}.dat'
        # Caricamento dati
        data = np.loadtxt(filepath, delimiter=",")

        beta = data[:, 0]  # Asse x
        susceptibility = data[:, 2]  # Asse y (suscettività)

        try:
            popt, _ = curve_fit(parabola, beta, susceptibility, p0=p0_values[model])
            A, B, C = popt

            # Salvataggio dei parametri
            fit_params[model]['A'].append(A)
            fit_params[model]['B'].append(B)
            fit_params[model]['C'].append(C)

        except RuntimeError:
            print(f"Fit fallito per {model} con L={L}")
            fit_params[model]['A'].append(1)
            fit_params[model]['B'].append(1)
            fit_params[model]['C'].append(1)

# Aggiungere il fit ai parametri B e C con una legge di potenza + offset
fit_params_power_B = {model: None for model in modelli}
fit_params_power_C = {model: None for model in modelli}

for model in modelli:
    # Fit per il parametro B
    try:
        popt_B, _ = curve_fit(power_law, L_array, fit_params[model]['B'], p0=[1, -1, 0])
        A_B, n_B, C_B = popt_B
        fit_params_power_B[model] = (A_B, n_B, C_B)
        print(f"Fit potenza per il parametro B (modello {model}): A = {A_B}, n = {n_B}, C = {C_B}")
    except RuntimeError:
        print(f"Fit potenza fallito per B nel modello {model}")

    # Fit per il parametro C
    try:
        popt_C, _ = curve_fit(power_law, L_array, fit_params[model]['C'], p0=[1, -1, 0])
        A_C, n_C, C_C = popt_C
        fit_params_power_C[model] = (A_C, n_C, C_C)
        print(f"Fit potenza per il parametro C (modello {model}): A = {A_C}, n = {n_C}, C = {C_C}")
    except RuntimeError:
        print(f"Fit potenza fallito per C nel modello {model}")


# Creazione dei subplot
fig, ax = plt.subplots(1, 2, figsize=(2*params['fig_width'], params['fig_height']))

for i, model in enumerate(modelli):
    ax[0].plot(L_array, fit_params[model]['B'], marker='o', color=colors(i), label=model, linestyle='none', markerfacecolor='white', markeredgewidth = params['line_width_axes'], zorder = 2)
    fit_B = power_law(L_dense, *fit_params_power_B[model])
    ax[0].plot(L_dense, fit_B, color=colors(i), label=f'Fit {model}', marker='none', linewidth=params['line_width_axes'], alpha=0.5, zorder = 0)

ax[0].set_xlabel('L')
ax[0].set_ylabel('Parametro B')
ax[0].legend(loc="upper right")
#ax[0].set_xticks(ticks=[0.43, 0.44, 0.45])

for i, model in enumerate(modelli):
    ax[1].plot(L_array, fit_params[model]['C'], marker='o', color=colors(i), label=model, linestyle='none', markerfacecolor='white', markeredgewidth = params['line_width_axes'], zorder = 2)
    fit_C = power_law(L_dense, *fit_params_power_C[model])
    ax[1].plot(L_dense, fit_C, color=colors(i), label=f'Fit {model}', marker='none', linewidth=params['line_width_axes'], alpha=0.5, zorder = 0)

ax[1].set_xlabel('L')
ax[1].set_ylabel('Parametro C')
ax[1].legend(loc="upper right")
#ax[1].set_xticks(ticks=[0.43, 0.44, 0.45])

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

# Impostazioni finali e salvataggio della figura
plt.tight_layout(pad=params['pad'])
plt.savefig('./figure_6.pdf', format='pdf')
plt.close(fig)
