import numpy as np
import matplotlib.pyplot as plt
import os

# Parameters
Nt = 200
times = np.linspace(20, 100, 8)
simbeta = 10.0
num_gaps = 5
num_mc = 1000

# Initialize gaps array
gaps = np.zeros((len(times), num_gaps))
gaps_std = np.zeros((len(times), num_gaps))
eta = simbeta / Nt
filename = f'Nt{int(Nt)}_simbeta{simbeta:.1f}.dat'

# Change directory to where files are stored
data_path = r'C:\Users\marco\Documents\metodi_numerici\modulo_3\analysis'
filepath = os.path.join(data_path, filename)
    
if not os.path.exists(filepath):
    raise FileNotFoundError(f"File {filepath} could not be opened")
    
# Read the file data
with open(filepath, 'r') as file:
    header = file.readline()  # Read and ignore header
    data = np.loadtxt(file)  # Read numerical data

# Extract the first column
first_column = data[:, 0]
second_column = data[:, 1]
    
# Extract necessary values
x2 = first_column[1]
x2_err = second_column[1]
x4 = first_column[3]
x4_err = second_column[3]
#x6 = first_column[4]


B = np.array([
        [first_column[-10], 0, first_column[-9], 0, first_column[-8]],
        [0, first_column[-7] - x2**2, 0, first_column[-6] - x2 * x4, 0],
        [first_column[-9], 0, first_column[-5], 0, first_column[-4]],
        [0, first_column[-6] - x2 * x4, 0, first_column[-3] - x4**2, 0],
        [first_column[-8], 0, first_column[-4], 0, first_column[-2]]
])

B_err = np.array([
        [second_column[-10], 0, second_column[-9], 0, second_column[-8]],
        [0, second_column[-7] + 2 * x2_err, 0, second_column[-6] + (x2_err/x2 + x4_err/x4) * x2*x4, 0],
        [second_column[-9], 0, second_column[-5], 0, second_column[-4]],
        [0, second_column[-6] + (x2_err/x2 + x4_err/x4) * x2*x4, 0, second_column[-3] + 2 * x4_err, 0],
        [second_column[-8], 0, second_column[-4], 0, second_column[-2]]
])
"""
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
"""
# Loop over Nt_array
for ix, t in enumerate(times):
    tau = t - 4
    """
    A = np.array([
        [first_column[5], 0, first_column[6]],
        [0, first_column[7] - x2**2, 0],
        [first_column[6], 0, first_column[8]]
    ])
    
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
    """
    A = np.array([
        [first_column[4 + 9*ix], 0, first_column[5 + 9*ix], 0, first_column[6 + 9*ix]],
        [0, first_column[7 + 9*ix] - x2**2, 0, first_column[8 + 9*ix] - x2 * x4, 0],
        [first_column[5 + 9*ix], 0, first_column[9 + 9*ix], 0, first_column[10 + 9*ix]],
        [0, first_column[8 + 9*ix] - x2 * x4, 0, first_column[11 + 9*ix] - x4**2, 0],
        [first_column[6 + 9*ix], 0, first_column[10 + 9*ix], 0, first_column[12 + 9*ix]]
    ])

    A_err = np.array([
        [second_column[4 + 9*ix], 0, second_column[5 + 9*ix], 0, second_column[6 + 9*ix]],
        [0, second_column[7 + 9*ix]  + 2 * x2_err, 0, second_column[8 + 9*ix] + (x2_err/x2 + x4_err/x4) * x2*x4, 0],
        [second_column[5 + 9*ix], 0, second_column[9 + 9*ix], 0, second_column[10 + 9*ix]],
        [0, second_column[8 + 9*ix] + (x2_err/x2 + x4_err/x4), 0, second_column[11 + 9*ix] + 2 * x4_err, 0],
        [second_column[6 + 9*ix], 0, second_column[10 + 9*ix], 0, second_column[12 + 9*ix]]
    ])
    # Loop over Monte Carlo simulations
    eigenvalue_samples = []

    for _ in range(num_mc):
        # Add random perturbation to A and B based on error matrix
        A_gauss = np.random.randn(*A.shape)
        A_gauss = (A_gauss + A_gauss.T) / 2

        B_gauss = np.random.randn(*B.shape)
        B_gauss = (B_gauss + B_gauss.T) / 2

        A_perturbed = A + A_gauss * A_err
        B_perturbed = B + B_gauss * B_err

        # Solve the generalized eigenvalue problem
        eigenvalues, _ = np.linalg.eig(np.linalg.inv(B_perturbed) @ A_perturbed)

        # Process eigenvalues and store them
        if all([np.isreal(e) and e > 0 for e in eigenvalues]):
            eigenvalues = -np.log(eigenvalues) / (tau * eta)
            eigenvalue_samples.append(np.sort(eigenvalues))

    # Compute the mean and std dev of the eigenvalues across the Monte Carlo simulations
    eigenvalue_samples = np.array(eigenvalue_samples)
    gaps[ix, :] = np.mean(eigenvalue_samples, axis=0)
    gaps_std[ix, :] = np.std(eigenvalue_samples, axis=0)
    
# Plot each column of gaps
plt.figure()
for col in range(num_gaps):
    plt.errorbar(times, gaps[:, col], yerr = gaps_std[:, col], capsize = 4)

# Customize plot
plt.xlabel('tau')
plt.ylabel('Gaps')
plt.title('Gaps')
plt.grid(True)
plt.show()
