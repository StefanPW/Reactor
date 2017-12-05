import matplotlib.pyplot as plt
import numpy as np
#import math
 
# calculation parameters
t = 100             # time [s]
nt = -1             #negative time
dt = 0.01           # timestep of calculation
 
# calculation variables
listN = []
listR = []
listC = []
listTf = []
listTc = []
listTcIN = []
listW = []
listDT = []
 
# reactor parameters
# initial conditions
PCM = 8
rho1 = PCM * 10E-5
L = 0.001 # generation lifetime
rho0 = -1*10E-5
S = 10E5
if rho0 < 0:
    ne = -S*L/rho0 # initial neutron population 3*10E17
ne = 3*10E17
n = ne
 
# neutronic parameters
hl = [0 ,0.0124 ,0.0305,0.111, 0.301, 1.14, 3.01]
beta =[0, 0.000215, 0.001424, 0.001274, 0.002568, 0.000748, 0.000273]
betaS = sum(beta)
print(betaS)
 
# initial precursor density
# e means initial equilibrium value
 
 
C1e = beta[1]*ne/(hl[1]*L)
C2e = beta[2]*ne/(hl[2]*L)
C3e = beta[3]*ne/(hl[3]*L)
C4e = beta[4]*ne/(hl[4]*L)
C5e = beta[5]*ne/(hl[5]*L)
C6e = beta[6]*ne/(hl[6]*L)
 
C1 = C1e
C2 = C2e
C3 = C3e
C4 = C4e
C5 = C5e
C6 = C6e
 
# thermal hydraulic parameters
 
#feedback constants
alphaF = -1*10E-6
alphaC = -1*10E-6
 
aF = 1*10E-11#E-11 # energy per neutron [J] 
mF = 50000 #mass of fuel [kg]
CpF = 30 #UO2 heat capacity (@300 K BTW) [J/kgK]
mC = 30000 #mass of coolant in RPV (assumed 10t)
CpC = 100 #H2O heat capacity (@) [j/kgK]
h = 1.5*10E5 #3 heat conductivity @ fuel coolant interface [W/K]
 
Wce = 50000 #coolant mass flux [kg/s]
W = Wce
 
TcINe = 273+27#5#275 #[K]
TcIN = TcINe
 
Tfe = TcIN + aF*ne/(2*W*CpC) + aF*ne/h #fuel temperature [K]
Tce = TcIN + (aF*ne)/(2*W*CpC)#coolant temperature [K]

print('')
print('init. feedwatr temp. =  ' + str(round(TcIN,2)) + ' K')
print('init. coolant temp. =   ' + str(round(Tce,2)) + ' K')
print('init. fuel temp. =      ' + str(round(Tfe,2)) + ' K')

Tf = Tfe
Tc = Tce

print('Initial equilibrium power =  ' + str(round(aF*ne*10E-9,2)) + ' GW')
print('...')
i = nt/dt
while (i<t/dt):
    #dusturbance introduction
    #reactitivity
    if i == 0:
         
        W = W #* 0.9
        TcIN = TcIN# - 10

    if i < 0:
         rho = rho0 + alphaF*(Tf-Tfe) + alphaC*(Tc-Tce)

    else:
        rho = rho1 + alphaF*(Tf-Tfe) + alphaC*(Tc-Tce)
     
    dn = ((rho-betaS)/L*n + hl[1]*C1 + hl[2]*C2\
    + hl[3]*C3 + hl[4]*C4 + hl[5]*C5 + hl[6]*C6 + S)*dt
  
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
    if rho == 0:
        TD = 0
    else:
        TD = L/rho + (betaS-rho)/(0.08*rho+rho)
    
    #if i%500 == 0:
    #    print('Td = ' + str(round(TD,0)) + ' s')
   # print (k)
    listN.append(n)
    listDT.append(TD*10)
    listTf.append(Tf)
    listTc.append(Tc)
    listW.append(W)
    #listC.append(0.005*C)            #values are for plot scaling
    listR.append(1*10E4*rho)
 
    i = i + 1

ax = np.arange(nt, t, dt)

print('final coolant temp =  ' + str(round(Tc,2)) + ' K;       deltaTc = ' + str(round(Tc-Tce,2)) + ' K')
print('final fuel temp =     ' + str(round(Tf,2)) + ' K;      deltaTf = ' + str(round(Tf-Tfe)) + ' K')
print('final power output =  ' + str(round(aF*n*10E-9,3)) + ' GW,    which is ' + str(round(n/ne,2)) + " [P1/P0]")
print ('')
if Tf > 1400:
    print('meltdown')
fig, ax1 = plt.subplots()
ax1.plot(ax, listTf, 'g', markersize = 1)
ax1.plot(ax, listTc, 'b')
ax1.plot(ax, listR, 'm')
ax1.plot(ax,listDT,'y')
ax2 = ax1.twinx()
ax2.plot(ax, listN, 'r', markersize = 1)


#ax1.plot(ax, listTcIN, 'b', markersize = 1)
#ax1	.plot(ax, listW, 'r', markersize = 1)
#plt.hold(True)
plt.xlabel('time [s]',fontsize = 16)
ax2.set_ylabel('N', fontsize = 16)
ax1.set_ylabel('Temperatures [deg. C], reactivity [pcm]', fontsize = 16)
plt.xlim(xmin=-1)
ax1.set_ylim(ymin=-10, ymax = 1200)
#ax2.set_yscale('log')
plt.show()
