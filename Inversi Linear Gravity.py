import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import inv

x0=[100,300,650,950]
z0=[150,200,100,200]
r=[100,100,100,100]
rho=[2000,9000,2000,3000]
x=np.arange(0,1000,50)

n=len(x)
m=len(rho)


k=6.67*10**-11
vol=[]
for i in range (m):
    vol.append(4/3*np.pi*r[i]**3)


G=np.zeros((n,m))

for i in range (n):
    for j in range (m):
        G[i][j]+= k*vol[j]*z0[j]/(((x[i]-x0[j])**2)+(z0[j]**2))

d=np.matmul(G,rho)*10e5

Gt=np.transpose(G)
a=inv(np.matmul(Gt,G))
b=np.matmul(Gt,d)
m=np.matmul(a,b)

dc=np.matmul(G,m)


p = np.arange(0, 360, 1)  # angle in degree
q = np.array(p * np.pi / 180)  # angle in radian

fig, (ax1, ax2) = plt.subplots(2, 1, sharex = True)
ax1.plot(x, d, "-")
ax1.plot(x,dc,'.')
ax1.grid()
ax1.set_ylabel("Gravity Anomaly (mGal)")
ax1.set_title("Gravity Anomaly of Burried Sphere Model")

ax2.fill(r[0] * np.cos(q) + x0[0], r[0] * np.sin(q) + z0[0], "-")
ax2.fill(r[1] * np.cos(q) + x0[1], r[1] * np.sin(q) + z0[1], "-")
ax2.fill(r[2] * np.cos(q) + x0[2], r[2] * np.sin(q) + z0[2], "-")
ax2.fill(r[3] * np.cos(q) + x0[3], r[3] * np.sin(q) + z0[3], "-")
ax2.axis('equal')
ax2.invert_yaxis()
ax2.grid()
ax2.set_xlabel("Distance (m)")
ax2.set_ylabel("Depth (m)")

plt.show()