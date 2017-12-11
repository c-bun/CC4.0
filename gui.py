'''
Initial framework adapted from:
https://sukhbinder.wordpress.com/2014/12/25/an-example-of-model-view-controller-design-pattern-with-tkinter-python/
'''

import tkinter as Tk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
import pandas as pd
# TODO importing this will run it. need to put multiprocess call in
# orthogonal_set_finder file??
from run_OSF import *


class Model():

    def __init__(self):
        self.xpoint = 200
        self.ypoint = 200
        self.res = None

    def calculate(self):
        x, y = np.meshgrid(np.linspace(-5, 5, self.xpoint),
                           np.linspace(-5, 5, self.ypoint))
        z = np.cos(x**2 * y**3)
        self.res = {"x": x, "y": y, "z": z}


class View():

    def __init__(self, master):
        self.frame = Tk.Frame(master)
        self.label = Tk.Label(
            master, text="Select a CSV file for processing:")
        self.label.grid(columnspan=2, sticky=Tk.W)

        self.pathBox = Tk.Text(master, height=1, width=20)
        self.pathBox.insert(Tk.END, 'path/to/CSV')
        self.pathBox.grid(columnspan=2, row=1, sticky=Tk.W)

        self.dfBox = Tk.Text(master)
        # self.dfBox.insert(END, str(self.df))
        self.dfBox.grid(columnspan=3, row=2)

        self.load_button = Tk.Button(
            # master, text="Load", command=self.load)
            master, text="Load")
        self.load_button.grid(row=1, column=3)

        self.close_button = Tk.Button(
            master, text="Close", command=master.quit)
        self.close_button.grid(row=3, column=1)
        # self.close_button.pack()

        # self.frame = Tk.Frame(master)
        # self.fig = Figure(figsize=(7.5, 4), dpi=80)
        # self.ax0 = self.fig.add_axes(
        #     (0.05, .05, .90, .90), axisbg=(.75, .75, .75), frameon=False)
        # self.frame.pack(side=Tk.LEFT, fill=Tk.BOTH, expand=1)
        # self.sidepanel = SidePanel(master)
        # self.canvas = FigureCanvasTkAgg(self.fig, master=self.frame)
        # self.canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
        # self.canvas.show()


# class SidePanel():
#
#     def __init__(self, root):
#         self.frame2 = Tk.Frame(root)
#         self.frame2.pack(side=Tk.LEFT, fill=Tk.BOTH, expand=1)
#         self.plotBut = Tk.Button(self.frame2, text="Plot ")
#         self.plotBut.pack(side="top", fill=Tk.BOTH)
#         self.clearButton = Tk.Button(self.frame2, text="Clear")
#         self.clearButton.pack(side="top", fill=Tk.BOTH)


class Controller():

    def __init__(self):
        self.root = Tk.Tk()
        self.model = Model()
        self.view = View(self.root)
        self.view.load_button.bind('<Button>', self.loaddf)
        self.df = None
        # self.view.sidepanel.plotBut.bind("<Button>", self.my_plot)
        # self.view.sidepanel.clearButton.bind("<Button>", self.clear)

    def loaddf(self, event):
        filename = Tk.filedialog.askopenfilename(
            # initialdir='./', title='Select CSV', filetypes=(('CSV files',
            # '*.csv')))
            initialdir='./', title='Select a CSV')
        try:
            df = pd.read_csv(filename, index_col=0)
        except:
            # TODO this should be a dialog box?
            print("Something went wrong with the import of {}."
                  "Please check the file/path.".format(filename))
            raise
        self.view.pathBox.insert(Tk.END, filename)
        self.view.dfBox.insert(Tk.END, str(df))
        self.df = df

    def runOSF(self, event):
        if self.df == None:
            # TODO this should be a dialog.
            print("you must first load a CSV file.")
        else:
            result = run_multiprocess(self.df, (2, 2),
                                      numProcesses=2,
                                      threshold=1,
                                      buffer_length=1000000, fxn=check_RMSs_from_submatrix)
            self.view.dfBox.insert(Tk.END, str(result))

    def run(self):
        self.root.title("CrossCompare")
        self.root.deiconify()
        self.root.mainloop()


if __name__ == '__main__':
    c = Controller()
    c.run()
