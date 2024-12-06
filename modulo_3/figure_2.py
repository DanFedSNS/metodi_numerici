import numpy as np
import matplotlib.pyplot as plt

# File path to the data
Nt = 200
simbeta = .714
g = 0

file_path = f'./analysis/fig2/wf_Nt{Nt}_simbeta{simbeta:.3f}_g{g:.3f}.dat'

# Load the data, skipping the header row
data = np.loadtxt(file_path, skiprows=1)

# Extract positions and occurrences
positions = data[:, 0]  # First column: positions
occurrences = data[:, 1]  # Second column: occurrences

occurrences /= sum(occurrences)

dx = (positions[-1] - positions[0]) / (len(positions) - 1)
mask = occurrences > 0

# Plot the histogram
plt.figure(figsize=(10, 6))
#plt.bar(positions, occurrences, width=dx, align='edge', color='blue', edgecolor='blue')
plt.plot(positions, occurrences)
plt.xlabel('Position')
plt.ylabel('Occurrences')
plt.title('Histogram of Occurrences by Position')
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.show()
