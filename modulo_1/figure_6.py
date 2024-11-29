import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os

def parabola(x, A, B, C):
    return A * (x - B)**2 + C

def power_law(x, A, n, C):
    return A * x**n + C

def load_params(filepath):
    params = {}
    with open(filepath, 'r') as f:
        for line in f:
            key, value = line.strip().split('=')
            params[key] = float(value) if '.' in value else int(value)
    return params

def load_data(filepath):
    data = np.genfromtxt(filepath, delimiter=" ", dtype=float, filling_values=np.nan)
    return data[:, 0], data[:, 7]

def perform_fit(model, x, y, p0):
    try:
        popt, _ = curve_fit(parabola, x, y, p0=p0)
        fit_curve = parabola(x, *popt)  # Calcola la curva fittata
        return popt, fit_curve
    except RuntimeError:
        print(f"Fit fallito per {model}")
        return [1, 1, 1], np.full_like(x, np.nan)


os.chdir(os.path.dirname(os.path.abspath(__file__)))
params = load_params('params.txt')

L_array = np.linspace(60, 150, 10, dtype=int)
L_dense = np.linspace(min(L_array), max(L_array), 200)
modelli = ["ising2d_sq_cluster", "ising2d_tri_cluster", "ising2d_hex_cluster"]
p0_values = {model: [500000, 0.272, 100] for model in modelli}

# Inizializzazione dei dizionari per i parametri di fit
fit_params = {model: {'A': [], 'B': [], 'C': []} for model in modelli}
fit_params_power_B = {}
fit_params_power_C = {}

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
colors = plt.get_cmap('tab10')

for model in modelli:
    for L in L_array:
        filepath = f'./data/analysis_{model}/L{L}.dat'
        beta, susceptibility = load_data(filepath)

        (A, B, C), fit_curve = perform_fit(model, beta, susceptibility, p0_values[model])
        fit_params[model]['A'].append(A)
        fit_params[model]['B'].append(B)
        fit_params[model]['C'].append(C)

for model in modelli:
    popt_B, _ = curve_fit(power_law, L_array, fit_params[model]['B'], p0=[1, -1, 0])
    fit_params_power_B[model] = popt_B
    popt_C, _ = curve_fit(power_law, L_array, fit_params[model]['C'], p0=[1, -1, 0])
    fit_params_power_C[model] = popt_C

    print(f"Esponente nu per {model}: {-1/popt_B[1]}")
    print(f"Esponente gamma per {model}: {-popt_C[1]/popt_B[1]}")

# Creazione dei plot
fig, ax = plt.subplots(2, 1, figsize=(params.get('fig_width', 10), 1.5*params.get('fig_height', 5)))

for i, model in enumerate(modelli):
    label = f'{model.split("_")[1]}'
    ax[0].plot(L_array, fit_params[model]['B'], 'o', color=colors(i), label=label, markerfacecolor='white')
    fit_B = power_law(L_dense, *fit_params_power_B.get(model, [0, 0, 0]))
    ax[0].plot(L_dense, fit_B, color=colors(i), linewidth=params.get('line_width_axes', 1), alpha=0.5)

    ax[1].plot(L_array, fit_params[model]['C'], 'o', color=colors(i), label=label, markerfacecolor='white')
    fit_C = power_law(L_dense, *fit_params_power_C.get(model, [0, 0, 0]))
    ax[1].plot(L_dense, fit_C, color=colors(i), linewidth=params.get('line_width_axes', 1), alpha=0.5)

for idx, label in enumerate(['$\\beta_p$', '$\\chi_p$']):
    ax[idx].set_xlabel('L', fontsize=params['font_size_axis'], labelpad=params['label_pad'])
    ax[idx].set_ylabel(label, fontsize=params['font_size_axis'], labelpad=params['label_pad'])
    ax[idx].legend(loc="upper right", fontsize=params['font_size_legend'])
    ax[idx].grid(True)

    for spine in ax[idx].spines.values():
        spine.set_linewidth(params.get('line_width_axes', 1))

    ax[idx].tick_params(axis='x', labelsize=params.get('font_size_ticks', 10), direction='in')
    ax[idx].tick_params(axis='y', labelsize=params.get('font_size_ticks', 10), direction='in')

plt.tight_layout(pad=params.get('pad', 1))
plt.savefig('./figure/figure_6.pdf', format='pdf')
plt.close(fig)
