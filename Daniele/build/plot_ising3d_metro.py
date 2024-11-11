import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

filename = 'ising3d.dat'
#L = np.loadtxt(filename, max_rows=1)
L = 10;
data = np.loadtxt(filename, skiprows=1)

energy = data[:, 0]  # Prima colonna: energia per sito
magnetization = data[:, 1]  # Seconda colonna: magnetizzazione per sito

# Crea un PDF per salvare i plot
with PdfPages('ising3d_histograms.pdf') as pdf:
    plt.figure()
    plt.hist(energy, bins=int(3*(L**3)+1), color='blue', edgecolor='black', alpha=0.7, range=[-3, +3])
    plt.title('Energy per site')
    plt.xlabel('Energy')
    plt.ylabel('Occurencies')
    plt.xlim(min(energy), max(energy))
    pdf.savefig()  # Salva il grafico nel PDF
    plt.close()

    plt.figure()
    plt.hist(magnetization, bins=int(L**3+1), color='red', edgecolor='black', alpha=0.7, range=[-1,1])
    plt.title('Magnetization per site')
    plt.xlabel('Magnetization')
    plt.ylabel('Occurrencies')
    plt.xlim(min(magnetization), max(magnetization))
    pdf.savefig()  # Salva il grafico nel PDF
    plt.close()

print("Istogrammi salvati nel file 'ising3d_histograms.pdf'.")
