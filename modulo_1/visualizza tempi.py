import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

file_path_list = ["ising2d_tri_cluster.dat", "ising2d_sq_cluster.dat", "ising2d_hex_cluster.dat"]

# Step 1: Parse the file to extract relevant columns
pattern = r"time = ([\d.]+) min, beta = ([\d.]+), L = (\d+),"
beta_c_array = [np.log(np.sqrt(3))/2, np.log(1 + np.sqrt(2)) / 2, np.log(2 + np.sqrt(3))/2]

fig, ax = plt.subplots()
for ix, file_path in enumerate(file_path_list):
    beta_c = beta_c_array[ix]
    # Read the file and extract data
    data = []
    with open(f"./data/time/{file_path}", 'r') as file:
        for line in file:
            match = re.search(pattern, line)
            if match:
                time = float(match.group(1))  # Time in minutes
                beta = float(match.group(2))  # Beta value
                L = int(match.group(3))  # System size L
                data.append((time, beta, L))

    # Step 2: Convert to a DataFrame
    df = pd.DataFrame(data, columns=["time", "beta", "L"])

    L_values = df["L"].unique()

    # Step 4: Plot the data
    avg_time = []
    for L in L_values:
        filtered_data = df[df["L"] == L]
        filtered_data["beta_diff"] = abs(filtered_data["beta"] - beta_c)
        closest_row = filtered_data.loc[filtered_data["beta_diff"].idxmin()]
    
        avg_time.append(closest_row["time"])

    ax.plot(L_values, avg_time, marker='o', linestyle="none", label = f"{file_path}")


plt.ylabel("Time (min)", fontsize=12)
plt.xlabel("L", fontsize=12)
plt.grid(True)
plt.legend()
    
plt.show()
