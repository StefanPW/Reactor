
#Main program
from reactor import *
import tkinter as tk
from tkinter import *
from tkinter.messagebox import showinfo
from PIL import ImageTk, Image
import time

class MainApplication(tk.Tk):
    def __init__(self,parent):
        tk.Tk.__init__(self,parent)
        self.parent = parent
        
        self.initialize()
        
        #self.win = tk.Tk()
        #width = win.winfo_screenwidth()*0.8
        #height = win.winfo_screenheight()*0.8
        #win.geometry(width+"x"+height)
        #win.title("Nuclear Reactor Simulator")

    def initialize(self):
    #    self.grid()

        self.entry = tk.Entry(self)
   #     self.entry.grid(column=1,row=0,sticky='EW')

        self.bar_menu = tk.Menu(self)
        self.menu_about = tk.Menu(self.bar_menu, tearoff=0)
        self.bar_menu.add_cascade(label="About", menu=self.menu_about)
        self.menu_about.add_command(label="About", command=self.callAbout)
        self.config(menu=self.bar_menu)
        self.labelText = tk.StringVar()
        self.label = tk.Label(self, anchor= "w",fg="black",bg="white", textvariable = self.labelText)
 #       self.label.grid(column=3,row=0,columnspan=2)
        
        self.closeButton = tk.Button(self, text=u"Close", command = self.destroy)
  #      self.closeButton.grid(column = 1, row = 3, columnspan = 2)

  #      self.grid_rowconfigure([0,1,2,3],weight = 1)
  #      self.grid_columnconfigure([0,1,2,3],weight = 1)

        self.reactorImage(name = "reactor.gif")
        self.reactorWindow.pack(anchor=NW, side = TOP)
        self.closeButton.pack(anchor=S, fill = BOTH, side = BOTTOM)

    def reactorImage(self, name):

        self.reactorWindow = ReactorWindow(self, name, width = 1500, height = 540, cursor = "plus", relief = RAISED)
       # self.reactorWindow.grid(column = 0, row = 0, rowspan = 3, columnspan = 3)
     	
        self.reactorWindow.configure(height = self.reactorWindow.fileName.height()>>1, width = self.reactorWindow.fileName.width()>>1)
        #self.fileName = tk.PhotoImage(file = name)
        
        #self.B = tk.Button(self, text ="Hello", command = self.helloCallBack)
        #self.C = tk.Canvas(self, bg="blue", height=250, width=300)
        #self.oval = tk.Canvas(self, bg="blue", height=250, width=300)
        #self.oval.create_polygon([10, 50, 45, 80, 99, 120], fill="gray")
        #self.coord = 10, 50, 240, 210
        #self.arc = self.C.create_arc(self.coord, start=0, extent=150, fill="red")

        #self.C.grid(column = 1, row = 1, sticky='N')
        #self.oval.grid(column = 1, row = 2, sticky='N')
        #self.config(menu=self.bar_menu)
        #self.B.grid(column = 3, columnspan=2, row = 2, sticky='N')


    def helloCallBack(self):
        showinfo( "Hello Python", "Hello World")

    def callAbout(self):
        showinfo(title="About", message="My Window")

class ReactorWindow(tk.Canvas):

    def __init__(self, parent, name, cnf = {}, **kw):
        tk.Canvas.__init__(self,parent)
        self.parent = parent

        self.initialize(name)

    def initialize(self, name):

        self.fileName = tk.PhotoImage(file = name)
        self.nowy = self.fileName.subsample(2,2)
        self.background = self.create_image(0,0, image = self.nowy, anchor = NW)

'''fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
ax2 = ax1.twinx()
reactor = Reactor()

def animate(i):
    
    for i in range(0,10):
        reactor.simulate()
    ax1.clear()
    ax2.clear()
    

    ax = np.arange(reactor.nt, (reactor.i-0.5)*reactor.dt, reactor.dt)

    #print('final coolant temp =  ' + str(round(Tc,2)) + ' K;       deltaTc = ' + str(round(Tc-Tce,2)) + ' K')
    #print('final fuel temp =     ' + str(round(Tf,2)) + ' K;      deltaTf = ' + str(round(Tf-Tfe)) + ' K')
    #print('final power output =  ' + str(round(aF*n*10E-9,3)) + ' GW,    which is ' + str(round(n/ne,2)) + " [P1/P0]")
    #print ('')
    #if Tf > 1400:
    #    print('meltdown')
    #print("Liczba ax: ", len(ax), self.nt, self.i*self.dt, self.dt)
        
    #print("Rozmiary tablic: ", len(self.listTf), len(self.listTc), len(self.listR), len(self.listDT), len(self.listN))
    ax1.plot(ax, reactor.listTf, 'g', markersize = 1)
    ax1.plot(ax, reactor.listTc, 'b')
    ax1.plot(ax, reactor.listR, 'm')
    ax1.plot(ax,reactor.listDT,'y')
    
    
    ax2.plot(ax, reactor.listN, 'r', markersize = 1)
    
    #ax1.plot(ax, listTcIN, 'b', markersize = 1)
    #ax1	.plot(ax, listW, 'r', markersize = 1)
    #plt.hold(True)

    #ax2.set_yscale('log')
    #plt.show() '''  

if __name__ == "__main__":
 
    #win = MainApplication(None)
    #win.title("Nuclear Reactor Simulator")
    #width = win.winfo_screenwidth()*0.8
    #height = win.winfo_screenheight()*0.8
    #win.geometry(str(int(width))+"x"+str(int(height)))

    #Teraz to idzie do obiektu graficznego no ale tymczasowo zostawie to tutaj;
    reactor = Reactor()
    
    i = 0
    ani = animation.FuncAnimation(reactor.fig, reactor.animate, interval=int(50), blit = True)#
    #ani=animation.FuncAnimation(reactor.fig, reactor.animate, interval=(1000/60))
    plt.show()
    #while (True):
        #reactor.simulate()
        #ani = animation.FuncAnimation(reactor.fig, reactor.animate, interval=(1000/60))
        #plt.show()
        #i += 1
        #time.sleep(1)

    '''ani = animation.FuncAnimation(fig, animate, interval=40)
    plt.xlabel('time [s]',fontsize = 16)
    ax2.set_ylabel('N', fontsize = 16)
    ax1.set_ylabel('Temperatures [deg. C], reactivity [pcm]', fontsize = 16)
    plt.xlim(xmin=-1)
    ax1.set_ylim(ymin=-10, ymax = 1200)
    plt.show()'''
    #win.mainloop()
    #ani = animation.FuncAnimation(reactor.fig, reactor.animate, interval=(1000/60))
    #plt.show()


#while (0):
    

#s = input("fdf");
    #print("Nuclear Reactor Simulator")

    #myFuelRod = FuelRod(10,0)
    #print(myFuelRod.fuel_count)
    #myReactor = Reactor(myFuelRod)
    #dt = 0.1#100ms simulation step
    #t = 0