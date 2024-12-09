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
    return data[:, 0], data[:, 7], data[:, 8]  # Restituisce anche gli errori

def perform_fit(model, x, y, yerr, p0):
    popt, pcov = curve_fit(parabola, x, y, sigma=yerr, p0=p0)
    fit_curve = parabola(x, *popt)
    perr = np.sqrt(np.diag(pcov))
    return popt, perr, fit_curve

os.chdir(os.path.dirname(os.path.abspath(__file__)))
params = load_params('params.txt')

L_array = np.linspace(60, 150, 10, dtype=int)
L_dense = np.linspace(min(L_array), max(L_array), 200)
# modelli = ["ising2d_tri_cluster", "ising2d_sq_cluster", "ising2d_hex_cluster"]
modelli = ["ising2d_tri_cluster", "ising2d_hex_cluster", "ising2d_sq_cluster"]


fit_params = {model: {'A': [], 'B': [], 'C': [], 'A_err': [], 'B_err': [], 'C_err': []} for model in modelli}
fit_params_power_B = {}
fit_params_power_C = {}
fit_errors_power_B = {}
fit_errors_power_C = {}

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
colors = plt.get_cmap('tab10')

for model in modelli:
    for L in L_array:
        p0_values = {
            "ising2d_tri_cluster": [5000000, np.log(np.sqrt(3)) / 2 - 0.2/L, 0.03*L**1.75], 
            "ising2d_sq_cluster": [5000000, np.log(1 + np.sqrt(2)) / 2 - 0.4/L, 0.05*L**1.75],  
            "ising2d_hex_cluster": [5000000, np.log(2 + np.sqrt(3)) / 2 - 0.5/L, 0.08*L**1.75],     
        }
        filepath = f'./data/analysis_{model}/L{L}.dat'
        beta, susceptibility, susc_err = load_data(filepath)

        (A, B, C), (A_err, B_err, C_err), fit_curve = perform_fit(model, beta, susceptibility, susc_err, p0_values[model])
        fit_params[model]['A'].append(A)
        fit_params[model]['B'].append(B)
        fit_params[model]['C'].append(C)
        fit_params[model]['A_err'].append(A_err)
        fit_params[model]['B_err'].append(B_err)
        fit_params[model]['C_err'].append(C_err)

for model in modelli:
    popt_B, pcov_B = curve_fit(power_law, L_array, fit_params[model]['B'], sigma=fit_params[model]['B_err'], p0=[1, -1, 0])
    fit_params_power_B[model] = popt_B
    fit_errors_power_B[model] = np.sqrt(np.diag(pcov_B))

    popt_C, pcov_C = curve_fit(power_law, L_array, fit_params[model]['C'], sigma=fit_params[model]['C_err'], p0=[1, -1, 0])
    fit_params_power_C[model] = popt_C
    fit_errors_power_C[model] = np.sqrt(np.diag(pcov_C))

    nu = -1 / popt_B[1]
    gamma = -popt_C[1] / popt_B[1]
    nu_err = fit_errors_power_B[model][1] / (popt_B[1]**2)
    gamma_err = np.sqrt((fit_errors_power_C[model][1] / popt_B[1])**2 + 
                        (popt_C[1] * fit_errors_power_B[model][1] / (popt_B[1]**2))**2)
    chi_2 = np.sum((fit_params[model]['C'] - power_law(L_array, *popt_C))**2 / np.array(fit_params[model]['C_err'])**2)
    beta_critico = popt_B[2]
    err_beta_critico = fit_errors_power_B[model][2]
    print(f"Temperatura critica per {model}: {beta_critico:.5f} ± {err_beta_critico:.5f}")
    print(f"Esponente nu per {model}: {nu:.3f} ± {nu_err:.3f}")
    print(f"Esponente gamma per {model}: {gamma:.3f} ± {gamma_err:.3f}")
    print(f"Chi quadro normalizzati per {model}: {chi_2/10:.3f}")

fig, axes = plt.subplots(4, 1, figsize=(params.get('fig_width', 10), 2 * params.get('fig_height', 5)),
                         gridspec_kw={'height_ratios': [3, 1, 3, 1]}, sharex=True)

# Fit e curve per B
for i, model in enumerate(modelli):
    color = colors(i)
    label = f'{model.split("_")[1]}'

    # Grafico di B
    fit_B = power_law(L_dense, *fit_params_power_B.get(model, [0, 0, 0]))
    residuals_B = (np.array(fit_params[model]['B']) - power_law(L_array, *fit_params_power_B[model])) / np.array(fit_params[model]['B_err'])

    axes[0].plot(L_array, fit_params[model]['B'], 'o', color=color, label=label, markerfacecolor='white')
    axes[0].plot(L_dense, fit_B, color=color, linewidth=params.get('line_width_axes', 1), alpha=0.5, zorder=1)
    axes[1].plot(L_array, residuals_B, 'o', color=color, markerfacecolor='white')

# Fit e curve per C
for i, model in enumerate(modelli):
    color = colors(i)
    label = f'{model.split("_")[1]}'

    # Grafico di C
    fit_C = power_law(L_dense, *fit_params_power_C.get(model, [0, 0, 0]))
    residuals_C = (np.array(fit_params[model]['C']) - power_law(L_array, *fit_params_power_C[model])) / np.array(fit_params[model]['C_err'])

    axes[2].plot(L_array, fit_params[model]['C'], 'o', color=color, label=label, markerfacecolor='white')
    axes[2].plot(L_dense, fit_C, color=color, linewidth=params.get('line_width_axes', 1), alpha=0.5, zorder=1)
    axes[3].plot(L_array, residuals_C, 'o', color=color, markerfacecolor='white')

# Personalizzazione grafici principali
axes[0].set_ylabel('$\\beta_p$', fontsize=params['font_size_axis'], labelpad=params['label_pad'])
axes[2].set_ylabel('$\\chi_p$', fontsize=params['font_size_axis'], labelpad=params['label_pad'])
axes[0].legend(loc="upper right", fontsize=params['font_size_legend'])
axes[2].legend(loc="upper right", fontsize=params['font_size_legend'])

# Personalizzazione residui
for ax in [axes[1], axes[3]]:
    ax.axhline(0, color='black', linewidth=0.8, linestyle='--')
    ax.set_ylabel('Residui', fontsize=params['font_size_axis'], labelpad=params['label_pad'])
    ax.grid(True)

# Assi condivisi e configurazione finale
axes[3].set_xlabel('L', fontsize=params['font_size_axis'], labelpad=params['label_pad'])

for ax in axes:
    ax.grid(True)
    ax.tick_params(axis='x', labelsize=params.get('font_size_ticks', 10), direction='in')
    ax.tick_params(axis='y', labelsize=params.get('font_size_ticks', 10), direction='in')
    for spine in ax.spines.values():
        spine.set_linewidth(params.get('line_width_axes', 1))

plt.tight_layout(pad=params.get('pad', 1))
plt.savefig('./figure/figure_6.pdf', format='pdf')
plt.close(fig)
