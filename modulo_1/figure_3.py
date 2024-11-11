import numpy as np
import matplotlib.pyplot as plt

# Array di valori L da considerare
L_array = np.linspace(10, 30, 3, dtype=int)

# Lista dei modelli
models = ["ising2d_sq_cluster", "ising2d_tri_cluster", "ising2d_hex_cluster"]

# Colormap per curve di colore differenti
colors = plt.get_cmap('tab10')  # Colormap

# Crea il layout per i subplot (3 righe e 1 colonna)
fig, ax = plt.subplots(3, 1, figsize=(10, 12))

for model_index, model in enumerate(models):
    # Subplot per ogni modello (ogni modello avrà il suo subplot)
    for i, L in enumerate(L_array):
        filepath = f'./analysis_{model}/L{L}.dat'

        data = np.loadtxt(filepath, delimiter=",")
        
        beta = data[:, 0]  
        specific_heat = data[:, 1]  # Calore specifico
        susceptibility = data[:, 2]  # Suscettività
        magn_abs_avg = data[:, 3]  # Magnetizzazione assoluta media
        energy_avg = data[:, 4]  # Energia media

        # Scegliamo di plottare la suscettività per ogni modello
        ax[model_index].plot(beta, susceptibility, color=colors(i), label=f'L={L}', marker='o', linestyle='-')
        ax[model_index].set_title(f"Suscettività - Modello {model.split('_')[2].upper()}")
        ax[model_index].set_xlabel("Beta")
        ax[model_index].set_ylabel("Suscettività")
        ax[model_index].legend()

# Adjust the layout and save the figure
plt.tight_layout()
plt.savefig('./figure_3.pdf', format='pdf')
