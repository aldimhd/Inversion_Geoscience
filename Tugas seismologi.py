import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import inv
#Nama : M Aldi
#NIM  : 101117011
#kELAS : GP-A
#Tugas : Pemrogaman Seismologi

data = np.genfromtxt('datagempa.dat')
STA = ['Stasiun1','Stasiun2','Stasiun3','Stasiun4']
xi = data [0:,1]
yi = data [0:,2]
ti = data [0:,3]
to = 0
vp = 2
data1=len(ti)

J = []
K = []
error_data = []
M_data = []
M = [50, 50] #ruang model untuk linearized inversion
error = 0

#Linearized_inversion
for j in range(0, 7): #jumlah iterasi
    error = 0
    for k in range(0, data1):
        J_1 = (np.abs(M[0] - xi[k])) / (vp * ((M[0] - xi[k]) ** 2 + (M[1] - yi[k]) ** 2) ** 0.5)
        J_2 = (np.abs(M[1] - yi[k])) / (vp * ((M[0] - xi[k]) ** 2 + (M[1] - yi[k]) ** 2) ** 0.5)
        J_3 = (to + (((M[0] - xi[k]) ** 2 + (M[1] - yi[k]) ** 2) ** 0.5) / vp)

        J.append([J_1, J_2])
        K.append(np.abs(J_3 - ti[k]))
        error += (J_3 - ti[k]) ** 2

    M = np.transpose(M) + inv(np.transpose(J).dot(J)).dot(np.transpose(J)).dot(K)
    M_data.append(M)
    error_data.append((1/data1) * (error ** 0.5))

result = np.where(error_data == np.min(error_data))
index = result[0][0]
print(M_data)
xo_real = M_data[index][0]
yo_real = M_data[index][1]

print("==" * 50)

print("Nilai koordinat episenter  xo =", xo_real, "dan yo =",yo_real )
# Hasilnya xo = 67.90393529970181 dan yo = 51.918037201807174
print("Nilai error sebesar", np.min(error_data))
# Hasilnya3.555309225172922


xo = np.arange(0,500,50)
yo = np.arange(0,500,50)
#Grid_Search
for j in range(0, len(xo)):
    for h in range (0, len(yo)):
        error = 0
        M = [xo[j],yo[h]] #ruang model untuk grid search
        for k in range(0, data1):
            Jac_1 = (np.abs(M[0] - xi[k])) / (vp * ((M[0] - xi[k]) ** 2 + (M[1] - yi[k]) ** 2) ** 0.5)
            Jac_2 = (np.abs(M[1] - yi[k])) / (vp * ((M[0] - xi[k]) ** 2 + (M[1] - yi[k]) ** 2) ** 0.5)
            Jac_3 = (to + (((M[0] - xi[k]) ** 2 + (M[1] - yi[k]) ** 2) ** 0.5) / vp)

            J.append([Jac_1, Jac_2])
            K.append(np.abs(Jac_3 - ti[k]))
            error += (Jac_3 - ti[k]) ** 2

            M_data.append(M)
            error_data.append((1/data1) * (error ** 0.5))

result = np.where(error_data == np.min(error_data))
index = result[0][0]
xo_real = M_data[index][0]
yo_real = M_data[index][1]

print("==" * 50)
print("Nilai koordinat episenter xo =", xo_real, "dan yo =", yo_real)
# Hasilnya xo = 50 dan yo = 50
print("Nilai error ", np.min(error_data))
# Hasilnya 0.5901750000000003

fig = plt.figure()
warna = fig.patch.set_facecolor("#FFFFFF")
ax1 = fig.add_subplot(1, 1, 1, facecolor='white')
for i in range(0, data1):
    a = np.arange(0, 360, 1)
    b = np.array(a * np.pi / 180)
    r = ((xi[i] - xo_real) ** 2 + (yi[i] - yo_real) ** 2) ** 0.5
    ax1.plot(r * np.cos(b) + xi[i], r * np.sin(b) + yi[i], "-")
    ax1.scatter(xi[i], yi[i], alpha=0.1, cmap='binary')
ax1.scatter(xo_real, yo_real, alpha=0.3, cmap='binary')
ax1.set_xlim(-2.5 * np.max(xi), 2.5 * np.max(xi), 0.1)
ax1.set_ylim(-2.5 * np.max(yi), 2.5 * np.max(yi), 0.1)
ax1.set_title("Coordinate Epicenter")
ax1.set_xlabel("X-Axis")
ax1.set_ylabel("Y-Axis")
ax1.grid()
plt.show()