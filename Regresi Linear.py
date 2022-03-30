import numpy as np
from numpy.linalg import inv
import csv
import matplotlib.pyplot as plt
import pandas as pd

print('Nama Kelompok :')
print('M.aldi                   (101117011)')
print('Indra Rivaldi Siregar    (101117045)')
print('Akbar Tanjung            (101117055)')
print('Afif Fakhri Sasti        (101116014)')
print(37 * '=')

#Baca_data_1
# dataset = pd.read_csv ('Data_Inversi.csv')
# x = dataset.iloc[:, :-1].values
# y = dataset.iloc[:, 1].values

# Baca_data_2
dataset = np.genfromtxt ('data_t_inv_1.txt')
x = dataset[ :,0]
y = dataset[ :,1]

#Matriks_Kernel
m = 2
G = np.zeros((len(x),m))
for j in range(len(x)):
    G [j][0]= 1
    G [j][1]= x[j]
print('Matriks kernel:')
print(G)
print(30 * '=')

#Model
Gt = np.transpose(G)
Ginv = inv(np.matmul(Gt,G))
c = np.matmul(Ginv,Gt)
m = np.matmul(c,y)
print('Model:')
print(m)
print(30 * '=')

#Regresi liniear
A = np.vstack([x, np.ones(len(x))]).T
b, a = np.linalg.lstsq(A, y, rcond=None)[0]
print('Nilai Gradein:',b)
print('Nilai Intercept',a)

#Sketsa Titik
fig = plt.figure()
plt.scatter(x,y,s=15,label= 'Original Data', color= 'yellow')
plt.plot(x, b*x + a ,label= 'Garis Regresi', color= 'blue')
fig.suptitle('Sketsa Titik')
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.show()

