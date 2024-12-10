import numpy as np
import matplotlib.pyplot as plt
import os
import glob

def load_params(filepath):
    parameters = {}
    with open(filepath, 'r') as file:
        for line in file:
            key, value = line.strip().split('=')
            if key in ['fig_width', 'fig_height', 'dpi']:
                parameters[key] = float(value)
            else:
                parameters[key] = float(value) if '.' in value else int(value)
    return parameters

def pert(g, n):
    return 3/2 * g * (n**2 + n + 1/2)

# Carichiamo i parametri
os.chdir(os.path.dirname(os.path.abspath(__file__)))
plot_params = load_params('params.txt')

# Configurazione generale di Matplotlib
plt.rc('font', family='serif')
color_palette = plt.get_cmap('tab10')

# Percorso ai file
data_directory = './analysis/fig456/'
data_filepaths = sorted(glob.glob(os.path.join(data_directory, '*')))  # Legge tutti i file nella directory

Nt = 422
times = np.linspace(1, 70, 70)
simbeta = 6.541
tau_zero = 4
num_gaps = 4
num_mc = 100
eta = simbeta / Nt

# Variabili per raccogliere i risultati
g_values = []
gap_samples = {}
gap_std = {}
ix_for_B = tau_zero - 1

# Lettura e calcolo dei dati
for index, data_filepath in enumerate(data_filepaths):
    with open(data_filepath, 'r') as data_file:
        # Lettura del valore di g
        for line in data_file:
            if line.startswith('Nt'):
                g_line = line.split(',')[2] 
                g_value = float(g_line.split('=')[1].strip())
                if g_value not in gap_samples:
                    gap_samples[g_value] = []
                    gap_std[g_value] = []
                g_values.append(g_value)
                break
        
        numerical_data = np.loadtxt(data_file)

    first_column = numerical_data[:, 0]
    second_column = numerical_data[:, 1]

    x2 = first_column[1]
    x2_err = second_column[1]
    x4 = first_column[3]
    x4_err = second_column[3]

    B = np.array([
        [first_column[4 + 9*ix_for_B], 0, first_column[5 + 9*ix_for_B], 0],
        [0, first_column[7 + 9*ix_for_B] - x2**2, 0, first_column[8 + 9*ix_for_B] - x2 * x4],
        [first_column[5 + 9*ix_for_B], 0, first_column[9 + 9*ix_for_B], 0],
        [0, first_column[8 + 9*ix_for_B] - x2 * x4, 0, first_column[11 + 9*ix_for_B] - x4**2]
    ])
    
    B_err = np.array([
        [second_column[4 + 9*ix_for_B], 0, second_column[5 + 9*ix_for_B], 0],
        [0, second_column[7 + 9*ix_for_B] + 2 * x2*x2_err, 0, second_column[8 + 9*ix_for_B] + x2_err*x4 + x4_err*x2],
        [second_column[5 + 9*ix_for_B], 0, second_column[9 + 9*ix_for_B], 0],
        [0, second_column[8 + 9*ix_for_B] + x2_err*x4 + x4_err*x2, 0, second_column[11 + 9*ix_for_B] + 2 * x4 * x4_err]
    ])

    for ix, t in enumerate(times):
        if ix <= ix_for_B:
            continue

        tau = t - tau_zero
        
        A = np.array([
            [first_column[4 + 9*ix], 0, first_column[5 + 9*ix], 0],
            [0, first_column[7 + 9*ix] - x2**2, 0, first_column[8 + 9*ix] - x2 * x4],
            [first_column[5 + 9*ix], 0, first_column[9 + 9*ix], 0],
            [0, first_column[8 + 9*ix] - x2 * x4, 0, first_column[11 + 9*ix] - x4**2]
        ])

        A_err = np.array([
            [second_column[4 + 9*ix], 0, second_column[5 + 9*ix], 0],
            [0, second_column[7 + 9*ix] + 2 * x2*x2_err, 0, second_column[8 + 9*ix] + x2_err*x4 + x4_err*x2],
            [second_column[5 + 9*ix], 0, second_column[9 + 9*ix], 0],
            [0, second_column[8 + 9*ix] + x2_err*x4 + x4_err*x2, 0, second_column[11 + 9*ix] + 2 * x4 * x4_err]
        ])


        eigenvalue_samples = []

        for _ in range(num_mc):
            A_gauss = np.random.randn(*A.shape)
            A_gauss = (A_gauss + A_gauss.T) / 2

            B_gauss = np.random.randn(*B.shape)
            B_gauss = (B_gauss + B_gauss.T) / 2

            A_perturbed = A + A_gauss * A_err
            B_perturbed = B + B_gauss * B_err

            eigenvalues, _ = np.linalg.eig(np.linalg.inv(B_perturbed) @ A_perturbed)

            if any([np.isreal(e) and e > 0 for e in eigenvalues]):
                eigenvalues = -np.log(eigenvalues) / (tau * eta)
                eigenvalue_samples.append(np.sort(eigenvalues))

        eigenvalue_samples = np.array(eigenvalue_samples)
        gap_samples[g_value].append(np.mean(eigenvalue_samples, axis=0))
        gap_std[g_value].append(np.std(eigenvalue_samples, axis=0))

# Calcolo delle medie per ogni valore di g
g_values_unique = np.unique(g_values)

mean_gaps = {}
std_gaps = {}

for g in g_values_unique:
    gap_array = np.array(gap_samples[g]) 
    gap_std_array = np.array(gap_std[g])  
    
    # Calcola la media e deviazione standard per ciascun gap
    #mask = np.isnan(gap_array) == False
    #if g == 0:
    #    print(gap_array[mask])

    mean_gaps[g] = np.nansum(gap_array * 1/(gap_std_array**2), axis = 0) / np.nansum(1/(gap_std_array**2), axis = 0)
    #np.average(gap_array[mask], weights=1 / gap_std_array[mask]**2, axis=0) 
    std_gaps[g] = np.sqrt(6) * np.sqrt(1 / np.nansum(1/(gap_std_array**2), axis=0))

# Creiamo una figura
fig, ax = plt.subplots(figsize=(plot_params['fig_width'], plot_params['fig_height']))

for col in range(num_gaps - 1):  # Loop sui gap
    mean_values = [mean_gaps[g][col] if np.isfinite(mean_gaps[g][col]) else -10 for g in g_values_unique]  # Media per il gap corrente
    std_values = [std_gaps[g][col] if np.isfinite(std_gaps[g][col]) else 0.01 for g in g_values_unique]  # Std per il gap corrente
    
    ax.errorbar(
        g_values_unique, mean_values, yerr=std_values,
        label=f"$\\Delta E_{col + 1}$", color=color_palette(col),
        linestyle='none', marker='.', markerfacecolor='white',
        markeredgewidth=plot_params['line_width_axes'], zorder=2
    )


g_dense = np.logspace(-2, -1, 100)
pert1 = 1 + pert(g_dense, 1) - pert(g_dense, 0)
pert2 = 2 + pert(g_dense, 2) - pert(g_dense, 0)
pert3 = 3 + pert(g_dense, 3) - pert(g_dense, 0)

ax.plot(g_dense, pert1, linestyle = "--")
ax.plot(g_dense, pert2, linestyle = "--")
ax.plot(g_dense, pert3, linestyle = "--")

for spine in ax.spines.values():
    spine.set_linewidth(plot_params['line_width_axes'])

ax.tick_params(axis='x', labelsize=plot_params['font_size_ticks'], 
               width=plot_params['line_width_axes'], direction='in')
ax.tick_params(axis='y', labelsize=plot_params['font_size_ticks'], 
               width=plot_params['line_width_axes'], direction='in')
ax.margins(x=0.00, y=0.00)
ax.grid(True, which='minor', linestyle=':', linewidth=plot_params['line_width_grid_minor'])
ax.grid(True, which='major', linestyle='--', linewidth=plot_params['line_width_grid_major'])

ax.set_xlabel("g", fontsize=plot_params['font_size_axis'], labelpad=plot_params['label_pad'])
ax.set_ylabel("$\\Delta E$", fontsize=plot_params['font_size_axis'], labelpad=plot_params['label_pad'])
ax.legend(loc='best', fontsize=plot_params['font_size_legend'])
ax.set_xscale('log')
ax.set_ylim([0.5, 7.5])
plt.tight_layout(pad=plot_params['pad'])
plt.savefig('./figure/figure_6.pdf', format='pdf')
plt.close(fig)
