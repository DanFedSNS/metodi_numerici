import numpy as np
import matplotlib.pyplot as plt
import os
import glob

def computejack(datajack, filepath, numberofbins, binsize, eta, sampleeff, num_vars, skip_lines, num_gaps, dj_eig):
    # Allocate memory for totals, subtotals, and temporary row data
    totals = np.zeros(num_vars)
    subtotals = np.zeros(num_vars)
    data_row = np.zeros(num_vars)
    times = np.linspace(1,70,70)
    ix_skip = 3

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
                    eigenvalue_samples.append(np.sort(eigenvalues))

            for kk in range(num_gaps):
                dj_eig[i, kk] = np.nanmean([row[kk] for row in eigenvalue_samples])

            
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
    dj_eig = np.zeros((numberofbins, num_gaps))

    computejack(datajack, filepath, numberofbins, binsize, eta, sampleeff, num_vars, skip_lines, num_gaps, dj_eig)

    ris, err = compute_statistics(dj_eig, numberofbins, num_gaps)

    print("Averages:", ris)
    print("Errors:", err)