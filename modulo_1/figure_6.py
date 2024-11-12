import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Definizione della parabola per il fit
def parabola(x, A, B, C):
    return A * (x - B)**2 + C

# Definizione della legge di potenza con offset
def power_law(x, A, n, C):
    return A * x**n + C

# Parametri
L_array = np.linspace(10, 30, 5, dtype=int)
colors = plt.get_cmap('tab10')
modelli = ["ising2d_sq_cluster", "ising2d_tri_cluster", "ising2d_hex_cluster"]

# Creazione dei dizionari per salvare i parametri del fit
fit_params = {model: {'A': [], 'B': [], 'C': []} for model in modelli}

# Definizione dei valori iniziali del parametro p0 per ogni modello
p0_values = {
    "ising2d_sq_cluster": [1, 0.44, 1],  # Valori per il modello ising2d_sq_cluster
    "ising2d_tri_cluster": [1, 0.23, 1],  # Valori per il modello ising2d_tri_cluster
    "ising2d_hex_cluster": [1, 0.64, 1]   # Valori per il modello ising2d_hex_cluster
}

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

# Creazione dei subplot
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Plot del parametro B
for i, model in enumerate(modelli):
    ax1.plot(L_array, fit_params[model]['B'], marker='o', color=colors(i), label=model)
ax1.set_xlabel('L')
ax1.set_ylabel('Parametro B')
ax1.legend()
ax1.grid()

# Plot del parametro C
for i, model in enumerate(modelli):
    ax2.plot(L_array, fit_params[model]['C'], marker='o', color=colors(i), label=model)
ax2.set_xlabel('L')
ax2.set_ylabel('Parametro C')
ax2.legend()
ax2.grid()

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

# Impostazioni finali e salvataggio della figura
plt.tight_layout()
plt.savefig('./figure_6.pdf', format='pdf')
