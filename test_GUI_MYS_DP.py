#!/usr/bin/python
# -*- coding: utf-8 -*-
from Tkinter import Tk, Text, TOP, BOTH, X, N, LEFT, Canvas, END
from ttk import Frame, Label, Entry, Button
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from subprocess import call
from numpy import loadtxt
import matplotlib.pyplot as plt
import os.path as op
import shlex as sh
class Example(Frame):

    def __init__(self, parent):
        Frame.__init__(self, parent)
        self.outfileNum = 0
        self.parent = parent
        self.initUI()

    def initUI(self):
        self.parent.title("Test the Gauss point DP MYS")
        self.pack(fill=BOTH, expand=True)
        self.fields = \
                'bulk_modulus', \
                'scale_hardening', \
                'max_strain_in', \
                'increment_strain', \
                'Nloop', \
                'initial_confinement', \
                'reference_pressure', \
                'modulus_n', \
                'cohesion', \
                'RMC_shape_k', \
                'dilation_angle_eta', \
                'diletion_scale'

        default_values = \
                        '1E7', \
                        '1E3', \
                        '1E-1', \
                        '1E-3', \
                        '1', \
                        '1E5', \
                        '1E5', \
                        '0.5', \
                        '0.0', \
                        '1.0', \
                        '1.0', \
                        '1.0'
        # ==================
        # Entries for User input:
        self.entries = []
        for idx,field in enumerate(self.fields):
            row = Frame(self)
            row.pack(fill=X)
            labl = Label(row, text=field, width=30)
            labl.pack(side=LEFT, padx=5, pady=5)
            entry = Entry(row)
            entry.insert(END, default_values[idx])
            entry.pack(fill=X, padx=5, expand=True)
            self.entries.append((field,entry))
            # print field

        # ==================
        # Button for calculation
        frameButtonCalc = Frame(self)
        frameButtonCalc.pack(fill=X)
        calcButton = Button(frameButtonCalc, text="calculate",
                        command=self.calculate)
        calcButton.pack(side=LEFT,padx=5, pady=5)
              
        # ==================
        # Raw Frame for plot
        self.canvasFrame = Frame(self)
        self.canvasFrame.pack(fill=BOTH, expand=True)

    def calculate(self):
        # Get the User input into the entries
        argv = []
        for entry in self.entries:
            field_value = float(entry[1].get())
            # print entry[0], '=', field_value
            argv.append(str(field_value))

        # Call the executable to run Gauss point
        print "Start calculation! "
        arg = ' '.join([str(x) for x in argv])
        command = "script -c './multi_ys_DP " + arg + " ' log" 
        call(sh.split(command))

        print "Start Plotting! "
        strain_stress = loadtxt('strain_stress.txt')
        strain = strain_stress[:,0]
        stress = strain_stress[:,1]

        # Refresh for the next calculate (changed the parameter)
        try:
            self.line.set_xdata(strain)
            self.line.set_ydata(stress)
            minY = min(stress)*1.05
            maxY = max(stress)*1.05
            minX = min(strain)*1.05
            maxX = max(strain)*1.05
            self.ax.set_ylim([minY, maxY])
            self.ax.set_xlim([minX, maxX])
            self.canvas.draw()
            self.outfileNum = self.outfileNum + 1
            self.figSave(argv)
            print "Delete the old plot and Re-Draw!"
        # In the first plot:
        except AttributeError:
            print "First time plot!"
            self.fig = Figure()
            self.ax = self.fig.add_subplot(111)
            self.line, = self.ax.plot(strain,stress)
            self.canvas = FigureCanvasTkAgg(self.fig,master=self.canvasFrame)
            self.canvas.draw()
            self.canvas.get_tk_widget().pack(side='top', fill='both', expand=1)
            self.figSave(argv)
            pass

    def figSave(self, values):
        figtitle = ""
        for idx,field in enumerate(self.fields):
            figtitle = figtitle + field +" = " + "{:.2E}".format(float(values[idx])) +'\n'  
        figtitle = figtitle + '\n\n'
        self.ax.set_title(figtitle)
        self.fig.savefig("strain_stress" + str(self.outfileNum) + ".png", bbox_inches='tight')




def main():
    root = Tk()
    root.geometry("600x600+2600+300")
    app = Example(root)
    root.mainloop()  


if __name__ == '__main__':
    main()  