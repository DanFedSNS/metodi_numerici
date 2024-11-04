import matplotlib.pyplot as plt
import numpy as np
import ast 

modello = "ising2d_metro"
L_array = np.loadtxt("L_array.txt", dtype = int)
beta_array = np.loadtxt("beta_array.txt", dtype = float)

specific_heat_dict = {(L, beta): 0 for L in L_array for beta in beta_array}
susceptibility_dict = {(L, beta): 0 for L in L_array for beta in beta_array}
magn_abs_avg_dict = {(L, beta): 0 for L in L_array for beta in beta_array}
energy_avg_dict = {(L, beta): 0 for L in L_array for beta in beta_array}
binder_cum_dict = {(L, beta): 0 for L in L_array for beta in beta_array}
sigma_energy_bin_dict = {(L, beta): 0 for L in L_array for beta in beta_array}
sigma_magn_bin_dict = {(L, beta): 0 for L in L_array for beta in beta_array}
autocorrelation_energy_dict = {(L, beta): 0 for L in L_array for beta in beta_array}
autocorrelation_magn_dict = {(L, beta): 0 for L in L_array for beta in beta_array}

for L in L_array:
    for beta in beta_array:
        variables = {}
        filepath = f'./analysis_{modello}/L{L}_beta{beta:.2f}.dat'

        with open(filepath, 'r') as file:  
            next(file)  #salta la prima
            for line in file:
                if "=" in line:
                    name, value = line.split('=')
                    variables[name.strip()] = value.strip()

        for key, value in variables.items():
            val = [float(x) for x in value.split(",")]
            if len(val) == 1:   #se è un solo valore non lo salva come array
                locals()[f"{key}_dict"][(L, beta)] = val[0]
            else:
                locals()[f"{key}_dict"][(L, beta)] = val
            
fig, ax = plt.subplots(ncols = 2, nrows = 2)
for L in L_array:
    ax[0,0].set_title("<|m|>")
    ax[0,0].plot(beta_array, [magn_abs_avg_dict[(L, beta)] for beta in beta_array], label = f"L = {L}")

    ax[0,1].set_title("<E>")
    ax[0,1].plot(beta_array, [energy_avg_dict[(L, beta)] for beta in beta_array])

    ax[1,0].set_title("chi")
    ax[1,0].plot(beta_array, [susceptibility_dict[(L, beta)] for beta in beta_array])

    ax[1,1].set_title("c_v")
    ax[1,1].plot(beta_array, [specific_heat_dict[(L, beta)] for beta in beta_array])


fig.legend()
plt.show()


"""
N_mis = m.shape[0]
m_avg = np.average(m)
m_sd = np.std(m)
en_avg = np.average(energy)
en_sd = np.std(energy)

plt.yscale("log")
plt.plot([C(m, m_avg, m_sd, i) for i in range(1000)])
plt.plot([C(energy, en_avg, en_sd, i) for i in range(1000)])
plt.show()


print(f"\nSono state effettuate {N_mis:.1e} misure:\n")
#print(f"m = {m_avg:.3f} +- {m_sd:.2e}\nE = {en_avg:.3f} +- {en_sd:.2e}\n")


fig, ax = plt.subplots(1, 2)
ax[0].hist(m, bins = 101)
ax[1].hist(energy, bins = 101)

plt.show()




plt.plot(m) #time series
plt.show()

config = np.loadtxt('config.dat', skiprows = 1, delimiter = ",", unpack=True)
config = config.T
L = int(np.sqrt(config.shape[1]))

for i in range(config.shape[0]):
    plt.imshow(config[i].reshape(L, L), interpolation='none', cmap = "hot")
    plt.title("m = %.2f, E = %.2f" % (m[i], energy[i]))
    plt.pause(0.0001)
  
plt.show()"""