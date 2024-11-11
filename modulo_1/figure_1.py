import numpy as np
import matplotlib.pyplot as plt

L_array = [20]
beta_couple = [0.42, 0.46]  # Choose a beta value for plotting
modello = "ising2d_sq_cluster"
# modello = "ising2d_hex_metro"
# modello = "ising2d_tri_metro"

colors = plt.get_cmap('tab10')

fig, ax = plt.subplots(2, 1, figsize=(8, 12))

# Define x-axis limits
x_limits_magnetization = [-1, 1]  # Modify these limits as needed

for j in range(2): 
    for i in range(len(L_array)):    
        beta = beta_couple[j]
        L = L_array[i]  
        valori_mag = np.linspace(-1, 1, L**2+1)  # esempio di valori conosciuti di magnetizzazione

        # Read the data for the chosen beta value
        filepath = f'./{modello}/L{L}_beta{beta:.3f}.dat'
        data = np.loadtxt(filepath, skiprows=1, delimiter=",")  # Load the data file

        magnetization = data[:, 0]  # Magnetization per site

        # Calcola la frequenza delle occorrenze di magnetizzazione ed energia
        magnetization_frequencies = [np.count_nonzero(np.isclose(magnetization, val, atol=0.0024)) for val in valori_mag]

        # Magnetization distribution (normalizzata)
        ax[j].plot(valori_mag, magnetization_frequencies, color=colors(i), label=f'L={L}', marker='none', linestyle='-', linewidth=0.5)
        ax[j].set_title(f'Magnetization Distribution for Beta={beta}')
        ax[j].set_xlabel('Magnetization per site')
        ax[j].set_ylabel('Occurrences')
        ax[j].legend()
        ax[j].set_xlim(x_limits_magnetization)  # Set x-axis limits for magnetization plot

# Save the figure in a PDF
plt.savefig('./figure_1.pdf', format='pdf')
plt.close(fig)


