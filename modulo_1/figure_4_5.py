import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Array di valori L da considerare
L_array = np.linspace(10, 30, 5, dtype=int)

# Lista dei modelli
algos = ["ising2d_tri_metro", "ising2d_tri_cluster"]

# Stili per i modelli
styles = ['o', 's']
colori = plt.get_cmap('tab10')

# Definizione dei limiti per il fit
x_lim_metro = (0, 150)  # Limiti per il metodo Metro
x_lim_cluster = (0, 15)  # Limiti per il metodo Cluster

# Creazione della figura con due subplot (uno sopra l'altro)
fig, (ax_metro, ax_cluster) = plt.subplots(2, 1, figsize=(10, 10), sharex=False)

# Indice della riga di beta fissato
index_beta_fixed = 4

# Liste per memorizzare i tempi caratteristici
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
        filepath = f'./analysis_{algo}/L{L}_autocorr.dat'
        data = np.loadtxt(filepath, delimiter=",")

        autocorr_m = data[index_beta_fixed, :]
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

        color = colori(i)
        label = f'{algo.split("_")[2].capitalize()} - L={L}'

        popt, _ = curve_fit(exp_decay, x_fit, y_fit, p0=[1, 10])
        A_fit, tau_fit = popt

        # Aggiunta del fit alla lista dei tau
        tau_values.append(tau_fit)

        # Plot del fit esponenziale
        if "metro" in algo:
            ax_metro.plot(x_fit, exp_decay(x_fit, *popt), '--', c=color)
        else:
            ax_cluster.plot(x_fit, exp_decay(x_fit, *popt), '--', c=color)

        # Plot dei dati originali
        if "metro" in algo:
            ax_metro.plot(x_indices, autocorr_m, styles[algo_index],
                          c=color,
                          label=label,
                          markerfacecolor='none')
        else:
            ax_cluster.plot(x_indices, autocorr_m, styles[algo_index],
                            c=color,
                            label=label,
                            markerfacecolor='none')

    # Aggiungi i tau alla lista corrispondente
    if "metro" in algo:
        tau_metro = tau_values
    else:
        tau_cluster = tau_values

# Impostazioni del subplot "metro"
ax_metro.set_title("Metodo Metro")
ax_metro.set_xlabel("Indice dei dati")
ax_metro.set_xlim(x_lim_metro)
ax_metro.set_ylabel("Autocorrelazione (M)")
ax_metro.legend(title="Modelli e L", loc="upper right", ncol=2)
ax_metro.grid()

# Impostazioni del subplot "cluster"
ax_cluster.set_title("Metodo Cluster")
ax_cluster.set_xlabel("Indice dei dati")
ax_cluster.set_xlim(x_lim_cluster)
ax_cluster.set_ylabel("Autocorrelazione (M)")
ax_cluster.legend(title="Modelli e L", loc="upper right", ncol=2)
ax_cluster.grid()

# Migliora il layout e salva la figura
plt.tight_layout()
plt.savefig('./figure_4.pdf', format='pdf')

# Creazione della nuova figura con i risultati del fit
fig2, ax = plt.subplots(figsize=(10, 6))

# Plot dei tempi caratteristici al variare di L
ax.plot(L_array, tau_metro, 'o-', label="Metro", color='blue')
ax.plot(L_array, tau_cluster, 's-', label="Cluster", color='red')

# Fit per i dati Metro
popt_metro, _ = curve_fit(power_law, L_array, tau_metro, p0=[1, 1])
alpha_metro, c_metro = popt_metro

# Fit per i dati Cluster
popt_cluster, _ = curve_fit(power_law, L_array, tau_cluster, p0=[1, 1])
alpha_cluster, c_cluster = popt_cluster

print(f"Fit Metro: alpha = {alpha_metro:.4f}, c = {c_metro:.4f}")
print(f"Fit Cluster: alpha = {alpha_cluster:.4f}, c = {c_cluster:.4f}")

# Plot del fit per Metro
L_fit = np.linspace(min(L_array), max(L_array), 100)
ax.plot(L_fit, power_law(L_fit, *popt_metro), '--', color='blue',
        label=f"Fit Metro: $\\tau = L^{{{alpha_metro:.2f}}}$")

# Plot del fit per Cluster
ax.plot(L_fit, power_law(L_fit, *popt_cluster), '--', color='red',
        label=f"Fit Cluster: $\\tau = L^{{{alpha_cluster:.2f}}}$")



ax.legend(title="Modelli e Fit", loc="upper left")

# Impostazioni della figura
ax.set_title("Tempo Caratteristico al variare di L")
ax.set_xlabel("L")
ax.set_ylabel("Tempo Caratteristico (τ)")
ax.legend(title="Modelli", loc="upper left")
ax.grid()

# Salva la figura
plt.tight_layout()
plt.savefig('./figure_5.pdf', format='pdf')

print("Fit esponenziali completati considerando i limiti definiti e figure salvate come 'figure_4.pdf' e 'figure_5.pdf'.")
