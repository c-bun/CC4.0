'''
Initial framework adapted from:
https://sukhbinder.wordpress.com/2014/12/25/an-example-of-model-view-controller-design-pattern-with-tkinter-python/
'''

import tkinter as Tk
from tkinter import messagebox
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
import pandas as pd
from orthogonal_set_finder import *


class Model():

    def __init__(self):
        self.df = None
        self.resultDf = None
        self.dimension = (2, 2)
        self.processes = 2
        self.seq_addition = False

    def run_CC(self):
        self.resultDf = run_multiprocess(
            self.df, self.dimension, numProcesses=self.processes, seq_addition=self.seq_addition)
        print('Done!: {}'.format(self.resultDf[:10]))
        return self.resultDf

    def open_CSV(self, filename):
        try:
            df = pd.read_csv(filename, index_col=0)
        except:
            print("Something went wrong with the import of {}."
                  "Please check the file/path.".format(filename))
            raise
        self.df = df
        return df


class View():

    def __init__(self, master):
        # Path box and load button frame
        loadFrame = Tk.Frame(master)
        loadFrame.grid(row=0, columnspan=4)
        self.label = Tk.Label(
            loadFrame, text="Select a CSV file for processing:")
        self.label.pack(side=Tk.LEFT)

        self.pathBox = Tk.Text(loadFrame, height=1,
                               relief=Tk.RIDGE, borderwidth=2)
        self.pathBox.insert(Tk.END, 'path/to/CSV')
        self.pathBox.pack(side=Tk.LEFT)

        self.load_button = Tk.Button(
            loadFrame, text="Load")
        self.load_button.pack(side=Tk.RIGHT)

        # Area for printing input file contents
        self.dfBox = Tk.Text(master, relief=Tk.RIDGE, borderwidth=2)
        self.dfBox.grid(columnspan=4, row=1)

        # Area for selecting options for the script run
        optionsFrame = Tk.Frame(master)
        optionsFrame.grid(row=2, columnspan=4)

        dimensionLabel = Tk.Label(optionsFrame, text="Dimension:")
        dimensionLabel.pack(side=Tk.LEFT)

        self.dimensionEntry = Tk.Entry(optionsFrame, width=2)
        self.dimensionEntry.insert(Tk.END, '02')
        self.dimensionEntry.pack(side=Tk.LEFT)

        spacer1 = Tk.LabelFrame(optionsFrame, width=10)
        spacer1.pack(side=Tk.LEFT)

        self.processLabel = Tk.Label(optionsFrame, text="Processes:")
        self.processLabel.pack(side=Tk.LEFT)

        self.processEntry = Tk.Entry(optionsFrame, width=2)
        self.processEntry.insert(Tk.END, '02')
        self.processEntry.pack(side=Tk.LEFT)

        spacer2 = Tk.LabelFrame(optionsFrame, width=10)
        spacer2.pack(side=Tk.LEFT)

        self.run_button = Tk.Button(optionsFrame, text="Run")
        self.run_button.pack(side=Tk.LEFT)

        # Area for printing result
        self.resultBox = Tk.Text(master, relief=Tk.RIDGE, borderwidth=2)
        self.resultBox.grid(columnspan=4, row=3)

        # Area for save and close buttons
        self.save_button = Tk.Button(master, text="Save As")
        self.save_button.grid(row=4, column=1)

        self.close_button = Tk.Button(
            master, text="Close", command=master.quit)
        self.close_button.grid(row=4, column=2)


class Controller():

    def __init__(self):
        self.root = Tk.Tk()
        self.model = Model()
        self.view = View(self.root)
        self.view.load_button.bind('<Button>', self.loaddf)
        self.view.run_button.bind('<Button>', self.runOSF)
        self.view.save_button.bind('<Button>', self.saveResult)
        # self.view.sidepanel.plotBut.bind("<Button>", self.my_plot)
        # self.view.sidepanel.clearButton.bind("<Button>", self.clear)

    def loaddf(self, event):
        filename = Tk.filedialog.askopenfilename(
            # initialdir='./', title='Select CSV', filetypes=(('CSV files',
            # '*.csv')))
            initialdir='./', title='Select a CSV')
        try:
            df = self.model.open_CSV(filename)
        except:
            messagebox.showerror(
                "Open file", "Something went wrong with the import of {}.".format(filename))
            raise

        self.view.pathBox.delete("1.0", Tk.END)
        self.view.pathBox.insert(Tk.END, filename)
        self.view.dfBox.delete("1.0", Tk.END)
        self.view.dfBox.insert(Tk.END, str(df))

    def runOSF(self, event):
        if self.model.df is None:
            # TODO this should be a dialog.
            print("you must first load a CSV file.")
        else:
            result = self.model.run_CC()
            self.view.resultBox.delete("1.0", Tk.END)
            self.view.resultBox.insert(Tk.END, str(result[:10]))

    def saveResult(self, event):
        if self.model.resultDf is None:
            messagebox.showwarning(
                "No data to save.", "You must first run an analysis.")
        else:
            savePath = Tk.filedialog.asksaveasfilename(
                initialdir='./', title='Select a save location.')
            try:
                # FIXME resultDf is a list, not a pd df. Need to fix that.
                self.model.resultDf.to_csv(savePath, index=False)
            except:
                messagebox.showerror(
                    "Save file", "Something went wrong when saving {}.".format(savePath))
                raise

    def run(self):
        self.root.title("CrossCompare")
        self.root.deiconify()
        self.root.mainloop()


if __name__ == '__main__':
    c = Controller()
    c.run()
