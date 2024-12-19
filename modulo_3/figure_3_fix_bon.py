import numpy as np
import matplotlib.pyplot as plt
import os
import glob

ix_skip = 3
times = np.linspace(1,70,70)

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


def computejack(datajack, filepath, numberofbins, binsize, eta, sampleeff, num_vars, skip_lines, num_gaps, dj_eig):
    # Allocate memory for totals, subtotals, and temporary row data
    totals = np.zeros(num_vars)
    subtotals = np.zeros(num_vars)
    data_row = np.zeros(num_vars)
    

    # Read the file
    with open(filepath, 'r') as fp:
        # Skip the header
        header = fp.readline()

        # Skip the first 'skip_lines' lines
        for _ in range(skip_lines):
            if not fp.readline():
                raise ValueError(f"The value of skip_lines = {skip_lines} is greater than the file length.")

        # Calculate the total sums for all variables
        for _ in range(sampleeff):
            line = fp.readline().split()
            if len(line) < num_vars:
                raise ValueError(f"Insufficient data in line: {line}")
            data_row = np.array([float(line[var]) for var in range(num_vars)])
            totals += data_row

        # Reset file pointer and re-read header to restart bin processing
        fp.seek(0)
        header = fp.readline()

        # Skip the first 'skip_lines' lines again
        for _ in range(skip_lines):
            if not fp.readline():
                raise ValueError(f"The value of skip_lines = {skip_lines} is greater than the file length.")

        # Loop over bins to compute jackknife subtotals
        for i in range(numberofbins):
            # Copy totals to subtotals
            subtotals[:] = totals

            # Subtract bin data from subtotals
            for _ in range(binsize):
                line = fp.readline().split()
                if len(line) < num_vars:
                    raise ValueError(f"Insufficient data in line: {line}")
                data_row = np.array([float(line[var]) for var in range(num_vars)])
                subtotals -= data_row

            # Normalize and assign values to `datajack`
            for var in range(num_vars):
                subtotals[var] /= (numberofbins - 1) * binsize
                datajack[i, var] = subtotals[var]


            eigenvalue_samples = []
            for ix, t in enumerate(times):
                if ix <= ix_skip:
                    continue
                tau = t - 4
                
                A = np.array([
                    [subtotals[4 + 9*ix], 0, subtotals[5 + 9*ix], 0],
                    [0, subtotals[7 + 9*ix] - subtotals[1]**2, 0, subtotals[8 + 9*ix] - subtotals[1] * subtotals[3]],
                    [subtotals[5 + 9*ix], 0, subtotals[9 + 9*ix], 0],
                    [0, subtotals[8 + 9*ix] - subtotals[1] * subtotals[3], 0, subtotals[11 + 9*ix] - subtotals[3]**2]
                ])

                B = np.array([
                    [subtotals[4 + 9*70], 0, subtotals[5 + 9*70], 0],
                    [0, subtotals[7 + 9*70] - subtotals[1]**2, 0, subtotals[8 + 9*70] - subtotals[1] * subtotals[3]],
                    [subtotals[5 + 9*70], 0, subtotals[9 + 9*70], 0],
                    [0, subtotals[8 + 9*70] - subtotals[1] * subtotals[3], 0, subtotals[11 + 9*70] - subtotals[3]**2]
                ])

                eigenvalues, _ = np.linalg.eig(np.linalg.inv(B) @ A)
                
                if any([np.isreal(e) and e > 0 for e in eigenvalues]):
                    eigenvalues = -np.log(eigenvalues) / (tau * eta)
                    eigenvalues = np.sort(eigenvalues)
                    eigenvalue_samples.append(eigenvalues)

                for kk in range(num_gaps):
                    dj_eig[t][i, kk] = eigenvalues[kk]   #dj_eig[i, kk] = np.nanmean([row[kk] for row in eigenvalue_samples])
            
def count_lines(filepath):
    with open(filepath, 'r') as file:
        return sum(1 for _ in file)
    
def compute_statistics(datajack, numberofbins, num_vars):
    num_res = num_vars
    ris = np.zeros(num_res)
    err = np.zeros(num_res)

    # Compute average
    for j in range(num_res):
        ris[j] = np.sum(datajack[:, j]) / numberofbins

    # Compute error
    for j in range(num_res):
        err[j] = np.sqrt(((ris[j] - datajack[:, j]) ** 2).sum() * (numberofbins - 1) / numberofbins)

    return ris, err

def process_file_header(filepath, skip_lines):
    sample = count_lines(filepath)

    with open(filepath, 'r') as fp:
        header = fp.readline()  # Read the first line

        # Extract sample and measevery from the header
        import re

        # Update sample to account for measevery and skip lines
        
        sample -= skip_lines

        if sample <= 0:
            raise ValueError("Invalid sample size after processing skip_lines.")

        # Define eta, binsize, numberofbins, and sampleeff
        simbeta = float(re.search(r"simbeta = (\d+\.\d+)", header).group(1))
        Nt = int(re.search(r"Nt = (\d+)", header).group(1))
        eta = simbeta / Nt
        binsize = sample // 100
        numberofbins = sample // binsize
        sampleeff = numberofbins * binsize

        return eta, binsize, numberofbins, sampleeff
    
if __name__ == "__main__":
    num_gaps = 4

    # Example usage
    skip_lines = 10
    filepath = './misure/Nt422_simbeta6.541_g0.000.dat'

    eta, binsize, numberofbins, sampleeff = process_file_header(filepath, skip_lines)

    num_vars = 4 + 71 * 9
    datajack = np.zeros((numberofbins, num_vars + num_gaps))
    dj_eig = {tau: np.zeros((numberofbins, num_gaps)) for tau in np.linspace(1,70,70)}

    computejack(datajack, filepath, numberofbins, binsize, eta, sampleeff, num_vars, skip_lines, num_gaps, dj_eig)

    ris = np.zeros((70, num_gaps))
    err = np.zeros((70, num_gaps))
    for i in range(ix_skip+1, 70):
        ris[i, :], err[i, :] = compute_statistics(dj_eig[i], numberofbins, num_gaps)
    
    
    
    fig, ax = plt.subplots(2, 1, figsize=(plot_params['fig_width'], 1*plot_params['fig_height']), gridspec_kw={'height_ratios': [3, 1]}, sharex=True)

    for col in range(num_gaps):
        """valid_indices = [i for i, (gap, gap_err) in enumerate(zip(gap_samples[g_values_unique[jdx]], gap_std[g_values_unique[jdx]]))
            if not (np.isscalar(gap) or np.isnan(gap[col]) or np.isnan(gap_err[col]))]
        times_filtered = [times[i] * eta for i in valid_indices]
        gaps_filtered = [gap_samples[g_values_unique[jdx]][i][col] for i in valid_indices]
        gap_errs_filtered = [gap_std[g_values_unique[jdx]][i][col] for i in valid_indices]"""

        # Grafico con i dati filtrati
        ax[0].errorbar(times[ix_skip+1::3], ris[ix_skip+1::3, col], yerr = err[ix_skip+1::3, col],
                    label=f"$\\Delta E_{col}$", color=color_palette(col),
                    linestyle='none', marker='x', markerfacecolor='white',
                    markeredgewidth=plot_params['line_width_axes'], zorder=2)
        
        residui = (ris[:, col] - (col+1)) / err[:, col]
        ax[1].errorbar(times[ix_skip+1::3], residui[ix_skip+1::3], yerr = err[ix_skip+1::3, col],
                    label=f"$\\Delta E_{col}$", color=color_palette(col),
                    linestyle='none', marker='x', markerfacecolor='white',
                    markeredgewidth=plot_params['line_width_axes'], zorder=2)
        # ax.errorbar(times, [gap[col] for gap in gap_samples[g_values_unique[jdx]]], yerr=[gap_err[col] for gap_err in gap_std[g_values_unique[jdx]]], label = f"$\\lambda_{col}$", color=color_palette(col), linestyle='none', marker='.', markerfacecolor='white', markeredgewidth=plot_params['line_width_axes'], zorder=2)

    for spine in ax[0].spines.values():
        spine.set_linewidth(plot_params['line_width_axes'])

    #ax[0].set_title('g = ' + str(g_values_unique[jdx]))
    ax[0].tick_params(axis='x', labelsize=plot_params['font_size_ticks'], 
                width=plot_params['line_width_axes'], direction='in')
    ax[0].tick_params(axis='y', labelsize=plot_params['font_size_ticks'], 
                width=plot_params['line_width_axes'], direction='in')
    ax[0].margins(x=0.00, y=0.00)
    ax[0].grid(True, which='minor', linestyle=':', linewidth=plot_params['line_width_grid_minor'])
    ax[0].grid(True, which='major', linestyle='--', linewidth=plot_params['line_width_grid_major'])

    ax[0].set_xlabel("$\\tau$", fontsize=plot_params['font_size_axis'], labelpad=plot_params['label_pad'])
    ax[0].set_ylabel("$\\Delta E$", fontsize=plot_params['font_size_axis'], labelpad=plot_params['label_pad'])
    ax[0].legend(loc='best', fontsize=plot_params['font_size_legend'])
    #ax[0].set_yscale('log')
    ax[0].set_ylim([0.5, 4.5])

    plt.tight_layout(pad=plot_params['pad'])
    plt.savefig(f'./analysis/figure_3_fix.pdf', format='pdf')
    plt.close(fig)
