
#import numpy
import numpy as np
import matplotlib.pyplot as plt


# Number of elements and nodes.
e = int(input("Enter number of Elements: "))
n = int(input("Enter number of nodes: "))
k1 = np.zeros(shape=(n, n), dtype=float)
K = np.matrix(k1)
D = np.matrix([[1, -1], [-1, 1]])
B = np.array([-1, 1])
F = []
L1 = []
A1 = []
E1 = []
S1 = []
S2 = []

# given Data and stiffness matrices and global stiffness matrix
for i in range(e):
    el = i+1
    L = int(input("Enter Length of Element %d in mm : " % el))
    A = int(input("Enter Area of Element %d in mm^2 : " % el))
    E = int(input("Enter Young's Modulus of Element %d in GPa : " % el))
    L1.append(L)
    A1.append(A)
    E1.append(E)
    C = (A*E*10**3)/L
    ke = np.multiply(C, D)
    print("\n k%d =  \n" % el, ke, "\n")
    K[i:i+2, i:i+2] += ke
print("\n k= \n", K)
print("\n")
kr = K[1:, 1:]


# forces acting on nodes and gives nodal displacements
# forces are taken from 2nd nodes because first node is at fixed support of bar
m = 1
for j in range(n-1):
    m += 1
    f = int(input("Enter Force in kN on node %d: " % m))
    f = f*10**3
    F.append(f)
X = np.linalg.solve(kr, F)
X = np.insert(X, 0, 0)
X = np.matrix(X).T
F.insert(0, 0)
for l in range(n):
    o = l+1
    print("\n q%d = " % o, X[l, 0], "mm\n")


# elemental stresses
for k in range(e):
    P = E1[k]*10**3*(1/L1[k])*B
    P = np.matrix(P)
    S = np.matmul(P, X[k:k+2])
    S1.append(S[0, 0])
    print("\n Stress %d = " % (k+1), S[0, 0], "MPa")


# reaction at fixed support
R1 = (np.matmul(K[0], X)-F[0])/10**3
print("\n Reaction at Fixed node(R1)=  ", R1[0, 0], "kN")


# Claculation For Strain
for h in range(e):
    g = h+1
    Strain = S1[h]/E1[h]
    print("\n  Strain at Element %d =" % g, Strain, "\n")
    S2.append(Strain)

 #Plotting Stress-Strain Variation Curve
plt.title("STRESS vs STRAIN Variation Curve ")
plt.xlabel("STRAIN")
plt.ylabel("STRESS")
plt.grid(True)
plt.plot(S2, S1, 'r-', S2, S1, 'bo')
plt.show()
