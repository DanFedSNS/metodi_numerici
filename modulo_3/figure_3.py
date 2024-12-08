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

# Carichiamo i parametri
os.chdir(os.path.dirname(os.path.abspath(__file__)))
plot_params = load_params('params.txt')

# Configurazione generale di Matplotlib
plt.rc('font', family='serif')
color_palette = plt.get_cmap('tab10')

# Percorso ai file
data_directory = './analysis/medie/fig456_bis/'
data_filepaths = sorted(glob.glob(os.path.join(data_directory, '*')))  # Legge tutti i file nella directory

#times = np.linspace(10, 25, 16)
#times = np.concatenate((times, [50, 75, 100, 150, 200, 400]))
times = np.linspace(5, 100, 20)
tau_zero = 4
num_gaps = 4
num_mc = 1000

# Variabili per raccogliere i risultati
g_values = []
gap_samples = {}
gap_std = {}

# Lettura e calcolo dei dati
for index, data_filepath in enumerate(data_filepaths):
    if data_filepath.endswith(".c"):
        continue
    
    with open(data_filepath, 'r') as data_file:
        # Lettura del valore di g
        for line in data_file:
            if line.startswith('Nt'):
                Nt_line = line.split(',')[0] 
                Nt_value = float(Nt_line.split('=')[1].strip())
                simbeta_line = line.split(',')[1]
                simbeta_value = float(simbeta_line.split('=')[1].strip())

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
        [first_column[-7], 0, first_column[-6], 0],
        [0, first_column[-5] - x2**2, 0, first_column[-4] - x2 * x4],
        [first_column[-6], 0, first_column[-3], 0],
        [0, first_column[-4] - x2 * x4, 0, first_column[-2] - x4**2]
    ])

    B_err = np.array([
        [second_column[-7], 0, second_column[-6], 0],
        [0, second_column[-5] + 2 * x2_err, 0, second_column[-4] + (x2_err/x2 + x4_err/x4) * x2*x4],
        [second_column[-6], 0, second_column[-3], 0],
        [0, second_column[-4] + (x2_err/x2 + x4_err/x4) * x2*x4, 0, second_column[-2] + 2 * x4_err]
    ])

    for ix, t in enumerate(times):

        tau = t - tau_zero
        
        A = np.array([
            [first_column[4 + 6*ix], 0, first_column[5 + 6*ix], 0],
            [0, first_column[6 + 6*ix] - x2**2, 0, first_column[7 + 6*ix] - x2 * x4],
            [first_column[5 + 6*ix], 0, first_column[8 + 6*ix], 0],
            [0, first_column[7 + 6*ix] - x2 * x4, 0, first_column[9 + 6*ix] - x4**2]
        ])

        A_err = np.array([
            [second_column[4 + 6*ix], 0, second_column[5 + 6*ix], 0],
            [0, second_column[6 + 6*ix] + 2 * x2_err, 0, second_column[7 + 6*ix] + (x2_err/x2 + x4_err/x4) * x2*x4],
            [second_column[5 + 6*ix], 0, second_column[8 + 6*ix], 0],
            [0, second_column[7 + 6*ix] + (x2_err/x2 + x4_err/x4) * x2*x4, 0, second_column[9 + 6*ix] + 2 * x4_err]
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
            
            eta = simbeta_value / Nt_value
            if all([np.isreal(e) and e > 0 for e in eigenvalues]):
                eigenvalues = -np.log(eigenvalues) / (tau * eta)
                eigenvalue_samples.append(np.sort(eigenvalues))

        eigenvalue_samples = np.array(eigenvalue_samples)
        gap_samples[g_value].append(np.mean(eigenvalue_samples, axis=0))
        gap_std[g_value].append(np.std(eigenvalue_samples, axis=0))

# Calcolo delle medie per ogni valore di g
g_values_unique = np.unique(g_values)

# Creiamo una figura
fig, ax = plt.subplots(figsize=(plot_params['fig_width'], plot_params['fig_height']))

jdx = 1
for col in range(num_gaps-1):
    print(g_values_unique[jdx])
    ax.errorbar(times, [gap[col] for gap in gap_samples[g_values_unique[jdx]]], yerr=[gap_err[col] for gap_err in gap_std[g_values_unique[jdx]]], label = f"$\\lambda_{col}$", color=color_palette(col), linestyle='none', marker='s', markerfacecolor='white', markeredgewidth=plot_params['line_width_axes'], zorder=2)

for spine in ax.spines.values():
    spine.set_linewidth(plot_params['line_width_axes'])

ax.tick_params(axis='x', labelsize=plot_params['font_size_ticks'], 
               width=plot_params['line_width_axes'], direction='in')
ax.tick_params(axis='y', labelsize=plot_params['font_size_ticks'], 
               width=plot_params['line_width_axes'], direction='in')
ax.margins(x=0.00, y=0.00)
ax.grid(True, which='minor', linestyle=':', linewidth=plot_params['line_width_grid_minor'])
ax.grid(True, which='major', linestyle='--', linewidth=plot_params['line_width_grid_major'])

ax.set_xlabel("t", fontsize=plot_params['font_size_axis'], labelpad=plot_params['label_pad'])
ax.set_ylabel("$\\lambda$", fontsize=plot_params['font_size_axis'], labelpad=plot_params['label_pad'])
ax.legend(loc='best', fontsize=plot_params['font_size_legend'])
# ax.set_yscale('log')
# ax.set_ylim([2.025, 2.075])

plt.tight_layout(pad=plot_params['pad'])
plt.savefig('./figure/figure_3.pdf', format='pdf')
plt.close(fig)
