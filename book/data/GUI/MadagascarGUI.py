from Tkinter import *
import os

#global SOURCES
#global CWD
#global FILES

CWD = os.getcwd()
SOURCES=["amoco", "bpait", "marmousi", "marmousi2","sigsbee","pluto","other"]
#SELECTFILES = ['No files selected','Please Select']                #Global SELECTFILES
class Madagascar:

    def __init__(self, master):
    
###########################
###     Labels          ###
###########################
        Label(master, text="Sources",fg="blue",font=("Helvetica", 14)).grid(row=1,columnspan=3) 
        Label(master, text="Files",fg="blue",font=("Helvetica", 14)).grid(row=1,column=5,columnspan=5)
        Label(master, text="Header Info",fg="blue",font=("Helvetica", 14)).grid(row=20,column=0,columnspan=3)

        Label(master, text="n1").grid(row=21)
        Label(master, text="n2").grid(row=22)
        Label(master, text="n3").grid(row=23)
        
        Label(master, text="o1").grid(row=21,column=2)
        Label(master, text="o2").grid(row=22,column=2)
        Label(master, text="o3").grid(row=23,column=2)

        Label(master, text="d1").grid(row=21,column=4)
        Label(master, text="d2").grid(row=22,column=4)
        Label(master, text="d3").grid(row=23,column=4)
        
        Label(master, text="Label 1").grid(row=21,column=6)
        Label(master, text="Label 2").grid(row=22,column=6)
        Label(master, text="Label 3").grid(row=23,column=6)

        Label(master, text="Units 1").grid(row=21,column=8)
        Label(master, text="Units 2").grid(row=22,column=8)
        Label(master, text="Units 3").grid(row=23,column=8)

#        displayRow=7
#        for file in SELECTFILES:
#            Label(master,text=file,fg="red",font=("Helvetica",10,"italic")).grid(row=displayRow)
#            displayRow=displayRow+1
        
###########################
###   Header Entry      ###
###########################    
        self.n1 = Entry(master,width=5)
        self.n2 = Entry(master,width=5)
        self.n3 = Entry(master,width=5)

        self.n1.grid(row=21, column=1)
        self.n2.grid(row=22, column=1)
        self.n3.grid(row=23, column=1)

        self.o1 = Entry(master,width=5)
        self.o2 = Entry(master,width=5)
        self.o3 = Entry(master,width=5)

        self.o1.grid(row=21, column=3)
        self.o2.grid(row=22, column=3)
        self.o3.grid(row=23, column=3)
        
        self.d1 = Entry(master,width=5)
        self.d2 = Entry(master,width=5)
        self.d3 = Entry(master,width=5)

        self.d1.grid(row=21, column=5)
        self.d2.grid(row=22, column=5)
        self.d3.grid(row=23, column=5)
    
        self.label1 = Entry(master,width=9)
        self.label2 = Entry(master,width=9)
        self.label3 = Entry(master,width=9)
        
        self.label1.grid(row=21, column=7)
        self.label2.grid(row=22, column=7)
        self.label3.grid(row=23, column=7)
        
        self.delta1 = Entry(master,width=5)
        self.delta2 = Entry(master,width=5)
        self.delta3 = Entry(master,width=5)
        
        self.delta1.grid(row=21, column=9)
        self.delta2.grid(row=22, column=9)
        self.delta3.grid(row=23, column=9)

###########################
###     Buttons         ###
###########################
        self.select = Button(master, text="Select File(s)", fg="black", command=self.select)
        self.select.grid(row=3,column=13,sticky=E+W)
        self.load = Button(master, text="Load Source", command=self.load)
        self.load.grid(row=2, column=13,sticky=E+W)
        self.load = Button(master, text="Fetch Files")#, command=self.load)
        self.load.grid(row=4, column=13,sticky=E+W)

############################
### Source List Box      ###
############################
        #self.scrollbar=Scrollbar(frame,orient=VERTICAL)          #Activate for second Scrollbar on Sources
        self.sources = Listbox(master,selectmode=SINGLE,relief=RAISED,height=10,exportselection=0)
        for item in SOURCES:
            self.sources.insert(END, item)
        #self.scrollbar.config(command=self.sources.yview)
        #self.scrollbar.pack(side=LEFT,fill=Y)
        self.sources.grid(row=2,column=0,rowspan=5,columnspan=3)

############################
### Avail. File List Box ###
############################
        self.fileScrollBar=Scrollbar(master,orient=VERTICAL)
        self.files = Listbox(master,yscrollcommand=self.fileScrollBar.set,selectmode=MULTIPLE,
                              relief=RAISED,height=10,width=40)
        for item in ['To view files','load a source']:
            self.files.insert(END, item)
        self.fileScrollBar.config(command=self.files.yview)
        self.fileScrollBar.grid(column=9,row=2,rowspan=5,sticky=N+S)
        self.files.grid(column=4,row=2,rowspan=5,columnspan=5)
        self.current = None


########################################
#######                       ##########
####### FUNCTION DEFENITIONS  ##########
#######                       ##########
########################################

#--------------------
# Loads Files into GUI
#--------------------
    def load(self):
        global FILES
        source = self.sources.curselection()
        lengthCWD = len(CWD)
        lengthCWD = lengthCWD-3
        location = str(CWD[0:lengthCWD]) + str(SOURCES[int(source[0])])+'/FILES'
        input = open(location,'r')
        FILES = input.readlines()
        lengthFiles = len(FILES)
        filesOut = []
        global segyFiles, nativeFiles
        segyFiles = []
        nativeFiles = []
        for item in FILES[0:lengthFiles]:
            length = len(item)
            letterCount=0
            for letter in item[0:length]:
                letterCount= letterCount+1
                if letter is ':':
                    fileName = item[letterCount+3:length-1]
                    filesOut.append(fileName)
        for item in filesOut:
            length = len(item)
            extStart = length - 4
            extStart2 = length -2
            extStart3 = length - 1
            extensionOne = item[extStart:length]
            extensionTwo = item[extStart2:length]
            extensionThree = item[extStart3:length]
            if extensionOne == 'segy' or extensionOne == 'SGY' or extensionOne == 'SEGY' or extensionOne == 'sgy':
                segyFiles.append(item)
            if extensionTwo == 'HH' or extensionTwo == 'hh' or extensionThree == 'H' or extensionThree =='h':
                nativeFiles.append(item)
        FILES=segyFiles + nativeFiles
        self.files.delete(0, END)         # clear Existing entries
        for item in FILES:
            self.files.insert(END,item)
        self.fileScrollBar.config(command=self.files.yview)
        self.fileScrollBar.grid(column=9,row=2,rowspan=5,sticky=N+S)
        self.files.grid(column=4,row=2,rowspan=5,columnspan=5)
        self.current = None

#--------------------
# Selects Files to Fetch
#--------------------
    def select(self):
        selection = self.files.curselection()
        global SELECTFILES
        SELECTFILES=[]
        for item in selection:
            selectFiles =  FILES[int(item)]
            SELECTFILES.append(selectFiles)
#        print SELECTFILES
        fileOutputWindow = Toplevel()
        displayRow=1
        
        title=Label(fileOutputWindow,text='Selected File List')
        title.grid(row=0,column=0)        
        
        for file in SELECTFILES: 
            Label(fileOutputWindow,text=file,width=50,fg="red",font=("Helvetica",10,"italic")).grid(row=displayRow)
            displayRow=displayRow+1

      #  Label(master, text="Sources",fg="blue",font=("Helvetica", 14)).grid(row=1,columnspan=3) 


root = Tk()
title=Label(root,text='Madagascar',height=2,borderwidth=1,font=("Helvetica", 16,"bold"))
title.grid(row=0,column=0,columnspan=3)
app = Madagascar(root)
root.mainloop()
