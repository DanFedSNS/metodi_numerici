import numpy as np
import matplotlib.pyplot as plt
import os
import sys

# Array of different L values to consider
L_array = np.linspace(70, 150, 5, dtype=int)
#L_array = np.linspace(90, 100, 2, dtype=int)
model = "ising2d_tri_cluster"

def load_params(filepath):
    params = {}
    with open(filepath, 'r') as f:
        for line in f:
            key, value = line.strip().split('=')
            if key in ['fig_width', 'fig_height', 'dpi']:
                params[key] = float(value)
            else:
                params[key] = float(value) if '.' in value else int(value)
    return params

os.chdir(os.path.dirname(os.path.abspath(__file__)))
params = load_params('params.txt')

#plt.rc('text', usetex=True)
plt.rc('font', family='serif')

colors = plt.get_cmap('tab10') 

# Creiamo una figura con 4 subplot orizzontali
fig, ax = plt.subplots(2, 2, figsize=(2*params['fig_width']+0.25, 1.5*params['fig_height']))

ax = ax.flatten()
for i, L in enumerate(L_array):
    filepath = f'./data/analysis_{model}/L{L}.dat'
    data = np.genfromtxt(filepath, delimiter=" ", dtype=float, filling_values=np.nan)

    beta = data[:, 0]  
    specific_heat = data[:, 3]  # Specific heat
    susceptibility = data[:, 7]  # Susceptibility
    magn_abs_avg = data[:, 5]  # Average absolute magnetization
    energy_avg = data[:, 1]  # Average energy
    num = 2
    ax[0].plot(beta[::num], specific_heat[::num], color=colors(i), label=f'L={L}', marker='s', linestyle='none', markerfacecolor='white', 
        markeredgewidth = params['line_width_axes'], zorder = 2)
    ax[1].plot(beta[::num], susceptibility[::num], color=colors(i), label=f'L={L}', marker='s', linestyle='none', markerfacecolor='white', 
        markeredgewidth = params['line_width_axes'], zorder = 2)
    ax[2].plot(beta[::num], magn_abs_avg[::num], color=colors(i), label=f'L={L}', marker='s', linestyle='none', markerfacecolor='white', 
        markeredgewidth = params['line_width_axes'], zorder = 2)
    ax[3].plot(beta[::num], energy_avg[::num], color=colors(i), label=f'L={L}', marker='s', linestyle='none', markerfacecolor='white', 
        markeredgewidth = params['line_width_axes'], zorder = 2)


for ax_ in ax.flat:
    for spine in ax_.spines.values():
        spine.set_linewidth(params['line_width_axes'])
    ax_.tick_params(axis='x', labelsize=params['font_size_ticks'], 
                    width=params['line_width_axes'], direction='in')
    ax_.tick_params(axis='y', labelsize=params['font_size_ticks'], 
                    width=params['line_width_axes'], direction='in')
    ax_.margins(x=0.00, y=0.00)
    ax_.grid(True, which='minor', linestyle=':', linewidth=params['line_width_grid_minor'])
    ax_.grid(True, which='major', linestyle='--', linewidth=params['line_width_grid_major'])
    #ax_.set_xticks(ticks=[0.27, 0.275, 0.28])


ax[0].set_xlabel("$\\beta $", fontsize=params['font_size_axis'], labelpad=params['label_pad'])
ax[0].set_ylabel("$C$", fontsize=params['font_size_axis'], labelpad=params['label_pad'])
ax[0].legend(loc='best', fontsize=params['font_size_legend'])

ax[1].set_xlabel("$\\beta $", fontsize=params['font_size_axis'], labelpad=params['label_pad'])
ax[1].set_ylabel("$\chi'$", fontsize=params['font_size_axis'], labelpad=params['label_pad'])
ax[1].legend(loc='best', fontsize=params['font_size_legend'])

ax[2].set_xlabel("$\\beta $", fontsize=params['font_size_axis'], labelpad=params['label_pad'])
ax[2].set_ylabel("$\\langle |m|\\rangle$", fontsize=params['font_size_axis'], labelpad=params['label_pad'])
ax[2].legend(loc='best', fontsize=params['font_size_legend'])

ax[3].set_xlabel("$\\beta $", fontsize=params['font_size_axis'], labelpad=params['label_pad'])
ax[3].set_ylabel("$\\langle E \\rangle$", fontsize=params['font_size_axis'], labelpad=params['label_pad'])
ax[3].legend(loc='best', fontsize=params['font_size_legend'])

plt.tight_layout(pad=params['pad'])
plt.savefig('./figure/figure_2.pdf', format='pdf')
plt.close(fig)