import numpy as np
import matplotlib.pyplot as plt

# Parametri
L_array = [30]
beta_couple = [
    [0.40, 0.43, 0.45, 0.48],   # Per il modello "ising2d_sq_cluster"
    [0.25, 0.2625, 0.2812, 0.3], # Per il modello "ising2d_tri_cluster"
    [0.63, 0.65, 0.69, 0.71]    # Per il modello "ising2d_hex_cluster"
]
modelli = ["ising2d_sq_cluster", "ising2d_tri_cluster", "ising2d_hex_cluster"]
colori = plt.get_cmap('tab10')

# Creiamo una figura con 4 subplot orizzontali
fig, ax = plt.subplots(1, 4, figsize=(20, 5))

# Definire i limiti dell'asse x
x_limits_magnetization = [-1, 1]

# Loop sui valori di beta
for j in range(4):
    
    # Loop sui modelli
    for k, modello in enumerate(modelli):
        beta = beta_couple[k][j]
        # Loop sui valori di L
        for i in range(len(L_array)):
            L = L_array[i]
            valori_mag = np.linspace(-1, 1, L**2 + 1)

            # Legge i dati per il modello e il valore di beta scelto
            filepath = f'./{modello}/L{L}_beta{(beta):.4f}.dat'

            data = np.loadtxt(filepath, skiprows=1, delimiter=",")
            magnetization = data[:, 0]

            atol = 1 / (L**2)
            # Calcolo della frequenza delle occorrenze di magnetizzazione
            magnetization_frequencies = [
                np.count_nonzero(np.isclose(magnetization, val, atol=atol)) for val in valori_mag
            ]

            # Plot della distribuzione di magnetizzazione
            label = f'{modello.split("_")[1]} (L={L})'
            ax[j].plot(
                valori_mag, magnetization_frequencies,
                color=colori(k), label=label,
                marker='none', linewidth=0.5
            )

        # Imposta titolo, etichetta e limiti
        ax[j].set_title(f'Beta={beta}')
        ax[j].set_xlabel('Magnetization per site')
        ax[j].set_ylabel('Occurrences')
        ax[j].set_xlim(x_limits_magnetization)

    # Aggiungi la legenda
    ax[j].legend()

# Migliora la spaziatura tra i subplot
plt.tight_layout()

# Salva la figura in un PDF
plt.savefig('./figure_1.pdf', format='pdf')
plt.close(fig)
