from tkinter import Tk, Label, Button, W, Text, END
import pandas as pd


class MyFirstGUI:

    def __init__(self, master):
        self.master = master
        master.title("A simple GUI")
        self.df = pd.read_csv('./sample_data/terminal_test.csv', index_col=0)
        self.createWidgets()

    def createWidgets(self):
        self.label = Label(
            self.master, text="Select a CSV file for processing:")
        self.label.grid(columnspan=2, sticky=W)

        self.pathBox = Text(self.master, height=1, width=20)
        self.pathBox.insert(END, 'path/to/CSV')
        self.pathBox.grid(columnspan=2, row=1, sticky=W)

        self.dfBox = Text(self.master)
        #self.dfBox.insert(END, str(self.df.iloc[:6, 1:2]))
        self.dfBox.insert(END, str(self.df))
        self.dfBox.grid(columnspan=3, row=2)

        self.greet_button = Button(
            self.master, text="Load", command=self.greet)
        self.greet_button.grid(row=1, column=3)
        # self.greet_button.pack()

        self.close_button = Button(
            self.master, text="Close", command=self.master.quit)
        self.close_button.grid(row=3, column=1)
        # self.close_button.pack()

    def greet(self):
        print("Greetings!")

root = Tk()
my_gui = MyFirstGUI(root)
root.mainloop()
