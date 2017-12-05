import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
from matplotlib import style
from matplotlib.lines import Line2D

class Reactor(object):
    """This class represents a nuclear reactor"""
    
    def __init__(self, PCM = 8, L = 0.001, rho0 = -1*10E-5, S = 10E5, ne = 3*10E17, dt = 0.01, dataSize = 400):
        """Constructor"""
         #      TD, Tf, Tc, dTf, dTc, dC1-6, dn, rho, TcIN, W, C1-6  sie zmienia
        #Lot of additional constants etc.
        # calculation parameters (one chyba bd i tak przesylane)
        #self.t = 100             # time [s]
        self.nt = -1             #negative time
        self.dt = dt
        self.dataSize = dataSize
        # reactor parameters
        # initial conditions
        self.PCM = PCM
        self.rho1 = PCM * 10E-5
        self.L = L # generation lifetime
        self.rho0 = rho0
        self.S = S

        if rho0 < 0: #Ta linijka wlasicwi enie ma znaczneia, jezeli nizej i tak jest na stale przypisywana
            self.n = -S*L/rho0 # initial neutron population 3*10E17
        #self.ne = ne zbedne
        self.n = ne
        
        self.neutron = Neutron()
        self.cooler = Cooler(self.n) 
        
        # initial precursor density
        # e means initial equilibrium value
        self.C = self.n/self.L * (self.neutron.beta[1::]/self.neutron.hl[1::])
        
        # calculation variables
        self.listN = []
        self.listR = []
        self.listC = []
        self.listTf = []
        self.listTc = []
        self.listTcIN = []
        self.listW = []
        self.listDT = []

        self.i = self.nt/self.dt
        #PLOT
        self.fig, self.ax1 = plt.subplots()
        self.ax2 = self.ax1.twinx()

        self.startPoint = self.nt
        self.x = []
        self.lineTf = Line2D(self.x, self.listTf, color = 'g', animated=True, markersize = 1, linewidth = 2)
        self.lineTc = Line2D(self.x, self.listTc, color = 'b', animated=True)
        self.lineR = Line2D(self.x, self.listR, color = 'm',animated=True)
        self.lineDT = Line2D(self.x, self.listDT, color = 'y',animated=True)
        self.lineN = Line2D(self.x, self.listN, color = 'r',animated=True, markersize = 1, linewidth = 2)
        self.ax1.add_line(self.lineTf)
        self.ax1.add_line(self.lineTc)
        self.ax1.add_line(self.lineR)
        self.ax1.add_line(self.lineDT)
        self.ax2.add_line(self.lineN)

        style.use('fivethirtyeight')
        #self.fig = plt.figure()
        #self.ax1 = self.fig.add_subplot(1,1,1)
        #self.ax2 = self.ax1.twinx()
        #plt.xlabel('time [s]',fontsize = 16)
        self.ax1.set_xlabel('time [s]',fontsize = 16)
        self.ax2.set_ylabel('N', fontsize = 16)
        self.ax1.set_ylabel('Temperatures [deg. C], reactivity [pcm]', fontsize = 16)
        self.ax1.set_xlim(xmin=-1)
        self.ax1.set_ylim(ymin=-10, ymax = 1200)
        self.ax2.set_ylim(ymin=0, ymax = 10E18)

    def simulate(self, fps = 60):
        #dusturbance introduction
        #reactitivity     
        rho = self.cooler.calc_rho(self.i, self.rho0, self.rho1) 
        if rho == 0:
            self.TD = 0
        else:
            self.TD = self.L/rho + (self.neutron.betaS-rho)/(1.08*rho) #Tutaj zmieniam to: 0.08*self.rho+self.rho na 1.08 rho xd       

        self.update(rho)
   
        self.appendToLists(rho)
     
        self.i += 1
        if(self.i>self.dataSize): 
            self.popFromLists()
            self.startPoint += self.dt
        #if(self.i%100 == 0):
        #    print(len(self.listN))

    def update(self, rho):
        #Updating values #W sumie te updaty moznaby cale w 1 linijce
        self.n += ((rho-self.neutron.betaS)/self.L*self.n + sum(self.neutron.hl[1::]*self.C) + self.S)*self.dt
        self.C += ((self.n/self.L)*self.neutron.beta[1::] - self.neutron.hl[1::]*self.C)*self.dt
        self.cooler.update(self.n, self.dt)
    
    def appendToLists(self, rho):
        self.listN.append(self.n)
        self.listDT.append(self.TD*10)
        self.listTf.append(self.cooler.Tf)
        self.listTc.append(self.cooler.Tc)
        self.listW.append(self.cooler.W)
        #self.listC.append(0.005*C)            #values are for plot scaling
        self.listR.append(1*10E4*rho) 
    
    def popFromLists(self):
        self.listN.pop(0)
        self.listDT.pop(0)
        self.listTf.pop(0)
        self.listTc.pop(0)
        self.listW.pop(0)
        #self.listC.append(0.005*C)            #values are for plot scaling
        self.listR.pop(0)

    def init(self):
        style.use('fivethirtyeight')
        #self.fig = plt.figure()
        #self.ax1 = self.fig.add_subplot(1,1,1)
        #self.ax2 = self.ax1.twinx()
        #plt.xlabel('time [s]',fontsize = 16)
        self.ax1.set_xlabel('time [s]',fontsize = 16)
        self.ax2.set_ylabel('N', fontsize = 16)
        self.ax1.set_ylabel('Temperatures [deg. C], reactivity [pcm]', fontsize = 16)
        #self.ax1.set_xlim(xmin=-1)
        plt.xlim(xmin = -1)
        self.ax1.set_ylim(ymin=-10, ymax = 1200)
        #Plot options
        #self.ax1.add_line(self.lineTf)
        #self.ax1.add_line(self.lineTc)
        #self.ax1.add_line(self.lineR)
        #self.ax1.add_line(self.lineDT)
        #self.ax2.add_line(self.lineN)
        #print("Rozmiary: ",len(self.listN), len(self.x))
        '''
        self.lineTf, = self.ax1.plot([], [], lw=2)
        self.lineTc, = self.ax1.plot([], [], lw=2)
        self.lineR, = self.ax1.plot([], [], lw=2)
        self.lineDT, = self.ax1.plot([], [], lw=2)
        self.lineN, = self.ax2.plot([], [], lw=2)'''
        #return self.lineN,self.lineR,self.lineTf,self.lineTc,self.lineDT

    def animate(self, i):
        #self.startPoint = self.nt + (i>40)*(i*self.dt)
        for i in range(0,10):
            self.simulate()
            
        #speed = (1 - (np.e**(-(self.dt/1)*(self.i + 99))))*((self.i+99)*self.dt)/2
        #print("Speed: ",speed, (99+self.i)*self.dt, self.i)
        #self.ax1.set_xlim(xmin= -1 + (1)*self.i*self.dt, xmax = 9 + (1)*self.i*self.dt)
        #if(self.i%100==0):print("sf: ",(self.i>10/self.dt)*self.i*self.dt)
        #self.ax1.set_ylim(ymin=-10, ymax = 1200)
        self.x = np.arange(self.startPoint, (self.i-0.5)*self.dt, self.dt)
        self.ax1.set_xlim(xmin = min(self.x), xmax = max(self.x))
        self.ax2.set_ylim(ymin=0, ymax = max(self.listN)*2)
        self.lineTf.set_data(self.x,self.listTf)
        self.lineTc.set_data(self.x,self.listTc)
        self.lineR.set_data(self.x,self.listR)
        self.lineDT.set_data(self.x,self.listDT)
        self.lineN.set_data(self.x,self.listN)
        #self.ax2.plot(self.x, self.listN)
        self.ax1.figure.canvas.draw()
        #if(self.i%10==0): self.fig.canvas.draw()
        #plt.axis(xmin = -1+speed, xmax = 9+speed)
        return self.lineN,self.lineR,self.lineTf,self.lineTc,self.lineDT

    '''def init(self):
        self.lineTf.set_data([],[])
        self.lineTc.set_data([],[])
        self.lineR.set_data([],[])
        self.lineDT.set_data([],[])
        self.lineN.set_data([],[])
        return self.lineTf,self.lineTc,self.lineR, self.lineDT, self.lineN,
        ''self.ax1.plot([], [], 'g', markersize = 1)
        self.ax1.plot([], [], 'b')
        self.ax1.plot([], [], 'm')
        self.ax1.plot([], [],'y')
        self.ax2 = self.ax1.twinx()
        self.ax2.plot([], [], 'r', markersize = 1)
        return self.ax1, self.ax2,''

    def animate(self,i):

        for i in range(0,2):
            self.simulate()
        self.ax1.clear()
        self.ax2.clear()

        x = np.arange(self.nt, (self.i-0.5)*self.dt, self.dt)
        ''self.lineTf.set_data(x,self.listTf)
        self.lineTc.set_data(x,self.listTc)
        self.lineR.set_data(x,self.listR)
        self.lineDT.set_data(x,self.listDT)
        self.lineN.set_data(x,self.listN)
        return self.lineTf,self.lineTc,self.lineR, self.lineDT, self.lineN,
        ''
        #print('final coolant temp =  ' + str(round(Tc,2)) + ' K;       deltaTc = ' + str(round(Tc-Tce,2)) + ' K')
        #print('final fuel temp =     ' + str(round(Tf,2)) + ' K;      deltaTf = ' + str(round(Tf-Tfe)) + ' K')
        #print('final power output =  ' + str(round(aF*n*10E-9,3)) + ' GW,    which is ' + str(round(n/ne,2)) + " [P1/P0]")
        #print ('')
        #if Tf > 1400:
        #    print('meltdown')

        self.ax1.plot(x, self.listTf, 'g', markersize = 1)
        self.ax1.plot(x, self.listTc, 'b')
        self.ax1.plot(x, self.listR, 'm')
        self.ax1.plot(x,self.listDT,'y')
        self.ax2 = self.ax1.twinx()
        self.ax2.plot(x, self.listN, 'r', markersize = 1)
        #ax1.plot(ax, listTcIN, 'b', markersize = 1)
        #ax1	.plot(ax, listW, 'r', markersize = 1)
        #plt.hold(True)
        
        #ax2.set_yscale('log')
        #plt.show() 
        return self.ax1, self.ax2,  #'''

class Neutron(object):
    """This class represents neutron and its properties"""
    #Chyba nie ma sensu dawac ten 0 element, skoro nigdy sie go nie oblicza
    def __init__(self, hl=(0 ,0.0124 ,0.0305,0.111, 0.301, 1.14, 3.01), beta=(0, 0.000215, 0.001424, 0.001274, 0.002568, 0.000748, 0.000273), **kwargs):
        # neutronic parameters
        #self.hl = hl
        self.beta = np.array(beta)
        self.hl = np.array(hl)
        self.betaS = sum(self.beta)

        print("Sum of betas: ", self.betaS)
        return super().__init__(**kwargs)
        
class Cooler(object):
    """This class represents reactor cooler"""

    def __init__(self, ne, alphaF= -1*10E-6, alphaC= -1*10E-6, aF = 1*10E-11, mF = 50000,
    CpF = 30, mC = 30000, CpC = 100, h = 1.5*10E5, Wce = 50000, TcINe = 273+27, **kwargs):
        # thermal hydraulic parameters
 #Tutaj tez mozna te stale przemnozyc zeby nie mnozono ich ciagle
        #feedback constants
        self.alphaF = alphaF #te zmienne maja identyczne wartosci
        self.alphaC = alphaC
 
        self.aF = aF#E-11 # energy per neutron [J] 
        self.mF = mF#mass of fuel [kg]
        self.CpF = CpF #UO2 heat capacity (@300 K BTW) [J/kgK]
        self.mC = mC #mass of coolant in RPV (assumed 10t)
        self.CpC = CpC #H2O heat capacity (@) [j/kgK]
        self.h = h #3 heat conductivity @ fuel coolant interface [W/K]
 
        #self.Wce = Wce #coolant mass flux [kg/s]
        self.W = Wce
 
        #self.TcINe = TcINe#5#275 #[K]
        self.TcIN = TcINe
 
        #Mocna zastapic, to bd szybsze
        self.mF_mul_CpF = self.mF*self.CpF
        self.mC_mul_CpC = self.mC*self.CpC
        self.W_mul_CpC = self.W*self.CpC<<1
        #Tu chyba byl blad CpF zamiast CpC
        self.Tfe = self.TcIN + self.aF*ne/(2*self.W*self.CpF) + self.aF*ne/self.h #fuel temperature [K]
        self.Tce = self.TcIN + (self.aF*ne)/(2*self.W*self.CpC)#coolant temperature [K]

        print('')
        print('init. feedwatr temp. =  ' + str(round(self.TcIN,2)) + ' K')
        print('init. coolant temp. =   ' + str(round(self.Tce,2)) + ' K')
        print('init. fuel temp. =      ' + str(round(self.Tfe,2)) + ' K')

        self.Tf = self.Tfe
        self.Tc = self.Tce
        return super().__init__(**kwargs)

    def calc_rho(self, i, rho0, rho1):
        """Calculates rho"""
        if i == 0:      
            self.W = self.W #* 0.9
            self.TcIN = self.TcIN# - 10
            
        rho = 0
        if i < 0:
            rho = rho0 + self.alphaF*(self.Tf-self.Tfe) + self.alphaC*(self.Tc-self.Tce)
        else:
            rho = rho1 + self.alphaF*(self.Tf-self.Tfe) + self.alphaC*(self.Tc-self.Tce)
        
        return rho

    def update(self, n, dt):
        self.Tf += ((self.aF*n - self.h*(self.Tf-self.Tc))/(self.mF*self.CpF))*dt
        self.Tc += ((self.h*(self.Tf-self.Tc)-2*self.W*self.CpC*(self.Tc-self.TcIN))/(self.mC*self.CpC))*dt


#    zmienne:
#moc reaktora
#liczba neutronów
#temperatura paliwa
#temperatura wody na wejściu
#temperatura wody na wyjściu
#strumień wody (ew moc pompy)
#pozycja prętów sterujących