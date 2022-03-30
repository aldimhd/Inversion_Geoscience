import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import inv
# Modelling Self Potential

k = int(input('Input Momen Dielektrik : ')) #model asli
z = int(input('Input Kedalaman Anomali : ')) #model asli
x0 = int(input('Input Titik Awal Pengukuran : ')) #model asli
i = int(input('Input Interval Pengukuran : '))
x = np.arange(-i, i, 5)
teta = int(input('Input Sudut Polarisasi : '))
t = teta * (np.pi / 180)
r = int(input('Input jari-jari Model Bola : '))
th = np.arange(0, 360) #pemodelan sudut bola


n= len(x)
v=np.zeros((n)) #nilai pengukuran
for i in range(len(x)):
    v[i]+= k * (((((x[i] - x0) * np.cos(t)) + (z * np.sin(t))) / ((((x[i] - x0) ** 2) + (z ** 2)) ** 1.5)))

#model awal
xm = int(input('Input Titik Awal Pengukuran Model Awal : '))
zm = int(input('Input Kedalaman Anomali Model Awal: '))
km = int(input('Input Nilai Momen Dielektrik Model Awal : '))
M  = [xm,zm,km]
it=0
galat= 1 #%
rms= np.infty

while rms>=galat:
# while m0!=xm:
    it = it + 1  # jumlah iterasi
    # Model awal yang diiterasi
    k = km
    z = zm
    x0 = xm
    v1 = np.zeros((n))
    for i in range(len(x)):
        v1[i]+= k*(((((x[i] - x0) * np.cos(t)) + (z * np.sin(t))) / ((((x[i] - x0) ** 2) + (z ** 2)) ** 1.5)))

    J = np.zeros((n,3))
    for i in range(n):
        J[i][0]+= (-3*(x[i]-x0)*k*((x[i]-x0)*np.cos(t)+(z*np.sin(t))))/((((x[i] - x0) ** 2) + (z ** 2)) ** 2.5)
        J[i][1]+= (-3 * z * k * ((x[i] - x0) * np.cos(t) + (z * np.sin(t)))) / ( (((x[i] - x0) ** 2) + (z ** 2)) ** 2.5)
        J[i][2]+= (((x[i]-x0)*np.cos(t)+(z*np.sin(t)))/(((x[i]-x0)**2)+(z**2))**1.5)

    Jt = np.transpose(J)
    a = inv(np.matmul(Jt, J))
    b = (v-v1)
    c = np.matmul(a, Jt)
    m = M + np.matmul(c, b)
    x0 = m[0]
    z  = m[1]
    rms = (np.sqrt(np.mean(b ** 2))) #rms calculasi
    print('iteration=', it)
    print('rms',rms)
    print('x0',m[0])
    print('z' ,m[1])
    print('k' ,m[2])
    xunit = r * np.cos(th) + x0
    yunit = r * np.sin(th) + z

plt.subplot(311)
plt.plot(x, v, '.')
plt.plot(x, v1, 'r-')
plt.xlim(min(x), max(x))
plt.xlabel('Position (m)')
plt.ylabel('Delta v (volt)')
plt.title('Modeling SP - Model Sphere')
plt.subplot(313)
plt.fill(xunit, yunit, 'r')
plt.xlim(min(x), max(x))
plt.ylim(0, z + 20)
plt.gca().invert_yaxis()
plt.xlabel('Position (m)')
plt.ylabel('Depth (m)')
plt.show()