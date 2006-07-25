from Tkinter import *
import os
global SOURCES
global CWD
global FILES

CWD = os.getcwd()
SOURCES=["amoco", "bpait", "marmousi", "marmousi2","sigsbee","pluto","other"]
class Madagascar:
    
    def __init__(self, master):
        frame = Frame(master)
        frame.pack()

###########################
###    Buttons          ###
###########################
        self.select = Button(frame, text="Select Files", fg="black", command=self.select)
        self.select.pack(side=RIGHT)
        self.load = Button(frame, text="Load", command=self.load)
        self.load.pack(side=RIGHT)

############################
###  Canvas Settings     ###
############################
#        self.canvas=Canvas(frame)
#        self.canvas.create_text(20,10,text="String",fill="red")
#        self.canvas.pack(side=BOTTOM, padx=2, pady=2)             


############################
### Source List Box      ###
############################
        #self.scrollbar=Scrollbar(frame,orient=VERTICAL)          #Activate for second Scrollbar on Sources
        self.sources = Listbox(frame,selectmode=SINGLE,relief=RAISED,height=10,exportselection=0)
        for item in SOURCES:
            self.sources.insert(END, item)
        #self.scrollbar.config(command=self.sources.yview)
        #self.scrollbar.pack(side=LEFT,fill=Y)
        self.sources.pack(side=LEFT)
        self.sourceLabel = Label(frame,text="Source")
        self.sourceLabel.pack()

############################
### Avail. File List Box ###
############################
        self.fileScrollBar=Scrollbar(frame,orient=VERTICAL)
        self.files = Listbox(frame,yscrollcommand=self.fileScrollBar.set,selectmode=MULTIPLE,
                              relief=RAISED,height=10)
        for item in ['0','1']:
            self.files.insert(END, item)
        self.fileScrollBar.config(command=self.files.yview)
        self.fileScrollBar.pack(side=RIGHT,fill=Y)
        self.files.pack(side=LEFT)
        self.current = None
        
############################
###  Canvas Settings     ###
############################
#        self.canvas=Canvas(frame)
#        self.canvas.create_text(20,10,text="String",fill="red")
#        self.canvas.pack(side=RIGHT)       

########################################
#######                       ##########
####### FUNCTION DEFENITIONS  ##########
#######                       ##########
########################################

    def yview(self, *args):
        apply(self.b1.yview, args)
        apply(self.b2.yview, args)

#--------------------
# Loads Files into GUI
#--------------------
    def load(self):
        source = self.sources.curselection()
	lengthCWD = len(CWD)
	lengthCWD = lengthCWD-3
        location = str(CWD[0:lengthCWD]) + str(SOURCES[int(source[0])])+'/FILES'
        input = open(location,'r')
        FILES = input.readlines()
        lengthFiles = len(FILES)
        filesOut = []
        for item in FILES[0:lengthFiles]:
            length = len(item)
            letterCount=0
            for letter in item[0:length]:
                letterCount= letterCount+1
                if letter is ':':
                    fileName = item[letterCount+3:length-1]
                    filesOut.append(fileName)
        FILES=filesOut
        self.files.delete(0, END)         # clear Existing entries
        for item in FILES:
            self.files.insert(END,item)
        self.fileScrollBar.config(command=self.files.yview)
        self.fileScrollBar.pack(side=RIGHT,fill=Y)
        self.files.pack(side=LEFT)
        self.current = None
        return FILES

#--------------------
# Selects Files to Fetch
#--------------------
    def select(self):
        selection = self.files.curselection()
        for item in selection:
            print FILES
    
#############################################################################
#############################################################################
#############################################################################



################################
#  Define root frame           #
################################
root = Tk()
title=Label(root,text='Madagascar')
title.pack()
################################
# Call Madagascar              #
################################
app = Madagascar(root)
################################
#      Run loop                #
################################
root.mainloop()
