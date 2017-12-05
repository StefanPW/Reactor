import matplotlib.pyplot as plt
import numpy as np

# calculation parameters
t = 100             # time [s]
nt = -1                #negative time
dt = 0.1           # timestep of calculation

# calculation variables
listN = []
listR = []
listC = []
listTf = []
listTc = []
listTcIN = []
listW = []

# reactor parameters
# initial conditions
PCM = 0
rho1 = PCM * 10E-5
L = 0.001 # generation lifetime
rho = 0
rho0 = 0

ne = 5*10E14#-S*L/rho # initial neutron population

n = ne
# neutronic parameters
hl = [0 ,0.0124 ,0.0305,0.111, 0.301, 1.14, 3.01]
beta =[0, 0.000215, 0.001424, 0.001274, 0.002568, 0.000748, 0.000273]
betaS = sum(beta)
# initial precursor density
# e means initial equilibrium value


C1e = beta[1]*ne/(hl[1]*L)
C2e = beta[2]*ne/(hl[2]*L)
C3e = beta[3]*ne/(hl[3]*L)
C4e = beta[4]*ne/(hl[4]*L)
C5e = beta[5]*ne/(hl[5]*L)
C6e = beta[6]*ne/(hl[6]*L)


# thermal hydraulic parameters

#feedback constants
alphaF = -5*10E-7
alphaC = -5*10E-6

aF = 10E-6#E-11 # energy per neutron [J] 
mF = 85000 #mass of fuel [kg]
CpF = 300 #UO2 heat capacity (@300 K BTW) [J/kgK]
mC = 10000 #mass of coolant in RPV (assumed 10t)
CpC = 5000 #H2O heat capacity (@) [j/kgK]
h = 10E7 #heat conductivity @ fuel coolant interface [W/K]

Wce = 50000 #coolant mass flux [kg/s]
W = Wce

TcINe = 275+273 #[K]
TcIN = TcINe

Tfe = TcIN + aF*ne/(2*W*CpC) + aF*ne/h #fuel temperature [K]
Tce = TcIN + (aF*ne)/(2*W*CpC)#coolant temperature [K]

print(TcIN, Tce, Tfe)

Tf = Tfe
Tc = Tce
print (aF*ne*10E-9)
i = nt/dt
while (i<t/dt):
  
    dn = ((rho-betaS)/L*n + hl[1]*C1 + hl[2]*C2 + hl[3]*C3 \
        + hl[4]*C4 + hl[5]*C5 + hl[6]*C6)*dt
 
    dC1 = (beta[1]*n/L - hl[1] * C1)*dt
    dC2 = (beta[2]*n/L - hl[2] * C2)*dt
    dC3 = (beta[3]*n/L - hl[3] * C3)*dt
    dC4 = (beta[4]*n/L - hl[4] * C4)*dt
    dC5 = (beta[5]*n/L - hl[5] * C5)*dt
    dC6 = (beta[6]*n/L - hl[6] * C6)*dt
   
    dTf = ((aF*n - h*(Tf-Tc))/(mF*CpF))*dt
    
    dTc = ((h*(Tf-Tc)-2*W*CpC*(Tc-TcIN))/(mC*CpC))*dt
   

    n = n + dn
    C1 = C1 + dC1   
    C2 = C2 + dC2 
    C3 = C3 + dC3 
    C4 = C4 + dC4 
    C5 = C5 + dC5 
    C6 = C6 + dC6 
    
    Tf = Tf + dTf
    Tc = Tc + dTc
    
    listN.append(n*10E-5)
    listTf.append(Tf-273)
    listTc.append(Tc-273)
    listTcIN.append(TcIN-273)
    listW.append(W)
    #listC.append(0.005*C)            #values are for plot scaling
    listR.append(5*10E03*rho)

    i = i + 1

ax = np.arange(nt, t, dt)

print ('')
fig, ax1 = plt.subplots()
ax1.plot(ax, listN, 'r', markersize = 1)
ax2 = ax1.twinx()
ax2.plot(ax, listTf, 'g', markersize = 1)
ax2.plot(ax, listTc, 'b')
ax2.plot(ax, listR, 'm')
#ax2.plot(ax, listTcIN, 'b', markersize = 1)
#ax2.plot(ax, listW, 'r', markersize = 1)
#plt.hold(True)
plt.xlabel('time [s]',fontsize = 16)
ax1.set_ylabel('N', fontsize = 16)
ax2.set_ylabel('Temperatures [deg. C]', fontsize = 16)
plt.xlim(xmin=-1)
plt.ylim(ymin=-10, ymax = 1300)
#plt.yscale('log')
plt.show()
