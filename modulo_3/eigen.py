import numpy as np
import matplotlib.pyplot as plt
import os

# Parameters
Nt_array = np.linspace(80, 130, 11)
simbeta = 10.0
num_gaps = 3

# Initialize gaps array
gaps = np.zeros((len(Nt_array), num_gaps))

# Loop over Nt_array
for ix, Nt in enumerate(Nt_array):
    tau = (Nt / 4 - 1) - 4
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
    
    # Ensure sufficient data length
    if len(first_column) < 13:
        raise ValueError(f"Not enough data in file {filename}")
    
    # Extract necessary values
    x2 = first_column[1]
    x4 = first_column[3]
    x6 = first_column[4]
    
    # Define matrices A and B
    A = np.array([
        [first_column[5], 0, first_column[6]],
        [0, first_column[7] - x2**2, 0],
        [first_column[6], 0, first_column[8]]
    ])
    B = np.array([
        [first_column[9], 0, first_column[10]],
        [0, first_column[11] - x2**2, 0],
        [first_column[10], 0, first_column[12]]
    ])
    
    # Solve the generalized eigenvalue problem
    eigenvalues, _ = np.linalg.eig(np.linalg.inv(B) @ A)
    
    # Process eigenvalues
    eigenvalues = -np.log(eigenvalues) / (tau * eta)
    gaps[ix, :] = np.sort(eigenvalues)

# Print gaps
print(gaps)

# Plot each column of gaps
plt.figure()
for col in range(num_gaps):
    plt.plot(Nt_array, gaps[:, col], '-x')

# Customize plot
plt.xlabel('Nt')
plt.ylabel('Gaps')
plt.title('Gaps')
plt.grid(True)
plt.show()
