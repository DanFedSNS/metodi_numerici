import numpy as np
import matplotlib.pyplot as plt

# Array di valori L da considerare
L_array = np.linspace(70, 120, 6, dtype=int)

# Lista dei modelli
algos = ["ising2d_tri_metro", "ising2d_tri_cluster"]

# Stili per i modelli
styles = ['o', 's']
colors = ['blue', 'red']

# Creazione della figura
fig, ax = plt.subplots(1, 1, figsize=(10, 6))

beta_fixed = 0.275

for algo_index, algo in enumerate(algos):
    for i, L in enumerate(L_array):
        filepath = f'./analysis_{algo}/L{L}.dat'

        data = np.loadtxt(filepath, delimiter=",")
        
        beta = data[:, 0]  
        sigma_m = data[:, 6] 

        # Trova l'indice del valore di beta più vicino a beta_fixed
        beta_index = (np.abs(beta - beta_fixed)).argmin()

        # Plotta la suscettività per ogni modello al variare di L
        ax.plot(L, sigma_m[beta_index], styles[algo_index], color=colors[algo_index],
                label=f'{algo.split("_")[2].capitalize()} - L={L}' if i == 0 else "",
                markerfacecolor='none')

# Impostazioni del grafico
ax.set_xlabel("Dimensione del reticolo (L)")
ax.set_ylabel("Errore su M")
ax.legend()
ax.grid()

# Salva e mostra la figura
plt.tight_layout()
plt.savefig('./figure_3.pdf', format='pdf')
