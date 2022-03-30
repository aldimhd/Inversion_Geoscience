import numpy as np
import matplotlib.pyplot as plt
import random

# ======================================================

g = 6.673e-11 # Gravitional Constant
Rho = 2000
p = 90
l = 20
t = 40
iterasi = 10
V = p * l * t # Volume
Rho_est = random.sample(range(1000, 3000),iterasi) #Mc Sample untuk rho
Z_est = random.sample(range(100, 500),iterasi) #Mc Sample untuk kedalaman
X_est = random.sample(range(100, 750),iterasi) #Mc Sample untuk x
print('Nilai Rho Random',Rho_est)
print('Nilai Z random', Z_est)
print('Nilai x random', X_est)
# print(Rhof)
# print(len(Rhof))

# ======================================================

# Luas Daerah Pengamatan
Weight = 1500
Height = 500
space = 20
xo = np.arange(0, Weight, space)
# print(len(xo))
zo = np.zeros(len(xo))
# print(zo)
y = np.arange(10, Height, space)
x = 750
z = 300

# ======================================================

gobs = (10 ** 5) * g * V * Rho * (z - zo)/(((x - xo) ** 2 + (z - zo) ** 2) ** (3/2)) # Nilai G Observasi

# Pemodelan Inversi

error = np.infty
it = 0
for i in range (len(Rho_est)) :
    it = it + 1 # Banyaknya iterasi
    E = []
    gcal = (10 ** 5) * g * V * Rho_est[i] * (Z_est[i]-zo)/ ((X_est[i]-xo) ** 2 + (Z_est[i]-zo) ** 2) ** (3/2) # Nilai G Kalkulasi
    E.append((gcal-gobs)**2)
    rms = (np.sqrt(np.mean(E))) # Error Estimasi
    if error > rms:
        error = rms
        rho_est = Rho_est[i]
        z_est = Z_est[i]
        x_est = X_est[i]
    gcal = (10 ** 5) * g * V * rho_est * (z_est - zo) / ((x_est - xo) ** 2 + (z_est - zo) ** 2) ** (3 / 2)

    print('-' * 50)
    print('Iterasi Ke -',it)
    print('rho calculasi :',rho_est)
    print('Nilai Z estimasi :',z_est)
    print('Nilai X estimasi :',x_est)
    # print('rms error :',rms)
    print('nilai error :',error)

# ======================================================

    FS = np.zeros(shape=(len(zo), len(xo)))
    FS[40:50, 25:50] = rho_est # letak anomali dari hasil inversi

    # FS2 = np.zeros(shape=(len(zo), len(xo)))
    # FS2[30:45, 25:50] = rho_est  # letak anomali dari hasil inversi

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
    ax1.plot(xo, gobs, "r-", label='G Observasi')
    ax1.plot(xo, gcal, '.', label='G Kalkulasi (Monte Carlo)')
    ax1.legend(loc=1, prop={'size': 8})
    ax1.grid()
    ax1.set_ylabel("Gravity Anomaly (mGal)")
    ax1.set_title("Model Gravity - Anomaly Finite Slab")

    ax2.imshow(FS, aspect='auto', extent=[min(xo), max(xo), max(y), min(y)], cmap="Greys")
    # ax2.imshow(FS2, aspect='auto', extent=[min(xo), max(xo), max(y), min(y)], cmap="Greys")
    ax2.set_xlabel("Distance (m)")
    ax2.set_ylabel("Depth (m)")
    ax2.grid()
    plt.show()
