#Main program
from reactor import *
import tkinter as tk
from tkinter import *
from tkinter.messagebox import showinfo
from PIL import ImageTk, Image
import time

'''class MainApplication(tk.Tk):
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
        self.background = self.create_image(0,0, image = self.nowy, anchor = NW)'''

if __name__ == "__main__":
 
    #win = MainApplication(None)
    #win.title("Nuclear Reactor Simulator")
    #width = win.winfo_screenwidth()*0.8
    #height = win.winfo_screenheight()*0.8
    #win.geometry(str(int(width))+"x"+str(int(height)))

    #Teraz to idzie do obiektu graficznego no ale tymczasowo zostawie to tutaj;
    inter = 50#Animation will run with inter ms gap
    reactor = Reactor(inter)
    
    i = 0
    ani = animation.FuncAnimation(reactor.fig, reactor.animate, interval=inter, blit = True)#
    #ani=animation.FuncAnimation(reactor.fig, reactor.animate, interval=(1000/60))
    plt.show()
    #while (True):
        #reactor.simulate()
        #ani = animation.FuncAnimation(reactor.fig, reactor.animate, interval=(1000/60))
        #plt.show()
        #i += 1
        #time.sleep(1)

    #win.mainloop()
    #ani = animation.FuncAnimation(reactor.fig, reactor.animate, interval=(1000/60))
    #plt.show()


#while (0):
    

#s = input("fdf");
