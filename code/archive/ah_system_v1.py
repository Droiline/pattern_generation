import numpy as np
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import matplotlib
import sys

size = 1000
time_steps = 3000

activ = [[0]*size for _ in range(time_steps)]
inhib = [[0]*size for _ in range(time_steps)]
working_row = [[0,0]]*size

#initial values for activator and inhibitor
innito = 0.01
innita = 2
inniti = 2

#values for diffusion, production and removal rates and source density
#what is source density again?
dt = 0.5
dx = 0.4
diffa = 0.1
diffi = 0.0

stability = (dx*dx) / (2*dt)
print("diffusion rates should be <= ", stability)
print("Diffa: ", diffa)
print("Diffi: ", diffi)

if diffa >= stability or diffi >= stability:
    sys.exit("The system is not stable")

proda = 6
prodi = 4
remove_a = 0.02
remove_i = 0.0075
sdens = 0.015

rho_0 = 0.02
k = 0.0004
mu = 0.05
nu = 0.03
rho = rho_0 + np.random.rand(size)*0.14

#the reaction equations
def simpleMeinhardt_a(oldA, oldI, s, r_a):
    dA = s * ((oldA * oldA) / oldI) - r_a * oldA
    return dA

def simpleMeinhardt_i(oldA, oldI, s, r_i, p_i):
    dI = s * oldA * oldA - r_i * oldI + p_i
    return dI

def meinhardt1_a(oldA, oldI, k, r_a, rho, mu):
    aStar2 = ((oldA*oldA)/(1 + k * oldA*oldA)) + r_a
    dA = (rho * oldI * aStar2) - mu * oldA
    return dA

def meinhardt1_i(oldA, oldI, k, r_a, s, rho, nu):
    aStar2 = ((oldA*oldA)/(1 + k * oldA*oldA)) + r_a
    dI = s - rho * oldI * aStar2 - nu * oldI
    return dI

def meinhardt2_a(oldA, oldI, k, rho_0, rho, mu):
    aStar2 = (oldA*oldA)/(1 + k * oldA*oldA)
    dA = (rho * (aStar2 + rho_0))/oldI - mu * oldA
    return dA

def meinhardt2_i(oldA, oldI, k, r_i, s, rho, nu):
    aStar2 = (oldA*oldA)/(1 + k * oldA*oldA);
    dI = rho * aStar2 - nu * oldI + r_i;
    return dI

#set initial values
activ[0] = [innita] * size
inhib[0] = [inniti] * size
activ[0][int(size/2)] = innita + 1
activ[0][int(size/2)+1] = innita + 1
inhib[0][int(size/2)] = inniti + 1
inhib[0][int(size/2)+1] = inniti + 1

#start iterating from after initial conditions
for t in range(1,time_steps):
    for c in range(size):
        #deal with edge cases using the neutral flow method
        if c == 0:
            old_activ_l = activ[t-1][c]
            old_inhib_l = inhib[t-1][c]
        else:
            old_activ_l = activ[t-1][c-1]
            old_inhib_l = inhib[t-1][c-1]

        if c == size-1:
            old_activ_r = activ[t-1][c]
            old_inhib_r = inhib[t-1][c]
        else:
            old_activ_r = activ[t-1][c+1]
            old_inhib_r = inhib[t-1][c+1]

        #find the change in activator and inhibitor due to diffusion
        diffused_a = diffa * ((old_activ_l + old_activ_r - 2*activ[t-1][c]) / (dx*dx))
        diffused_i = diffi * ((old_inhib_l + old_inhib_r - 2*inhib[t-1][c]) / (dx*dx))
        #find the change in activator and inhibitor due to reaction
#        reacted_a = meinhardt1_a(activ[t-1][c], inhib[t-1][c], k, remove_a, rho[c], mu)
#        reacted_i = meinhardt1_i(activ[t-1][c], inhib[t-1][c], k, remove_a, sdens, rho[c], nu)
        #print("reacted a and i: (", reacted_a, ", ", reacted_i, ")")
#        reacted_a = simpleMeinhardt_a(activ[t-1][c], inhib[t-1][c], sdens, remove_a)
#        reacted_i = simpleMeinhardt_i(activ[t-1][c], inhib[t-1][c], sdens, remove_i, prodi)
        reacted_a = meinhardt2_a(activ[t-1][c], inhib[t-1][c], k, rho_0, rho[c], mu)
        reacted_i = meinhardt2_i(activ[t-1][c], inhib[t-1][c], k, remove_i, sdens, rho[c], nu)
        #find new activator value
        activ[t][c] = activ[t-1][c] + dt * (diffused_a + reacted_a)
        inhib[t][c] = inhib[t-1][c] + dt * (diffused_i + reacted_i)
        #test the diffusion
#        activ[t][c] = activ[t-1][c] + dt * (diffused_a + 0)
#        inhib[t][c] = inhib[t-1][c] + dt * (diffused_i + 0)

#plt.plot(activ[2],inhib[2])
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#ax.set_xlabel("x axis")
#ax.set_ylabel("y axis")
#ax.set_zlabel("z axis")

#ax.plot_wireframe(activ, activ, activ)

fig, axes = plt.subplots(1, 2)
axes[0].set_title('Activator levels')
im1 = axes[0].imshow(activ, interpolation='nearest', cmap='afmhot')
axes[1].set_title('Inhibitor levels')
im2 = axes[1].imshow(inhib, interpolation='nearest', cmap='afmhot')
#cm1 = fig.colorbar(im1)
cm2 = fig.colorbar(im2)
plt.show()


#print("(activator, inhibitor)")
#for t in range(time_steps):
#    print("t: ", t)
#    for c in range(size):
#        print("(", activ[t][c], ", ", inhib[t][c], ")", end = " ")
#    print()
