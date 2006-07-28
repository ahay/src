from Tkinter import *
import os
#global SOURCES
#global CWD
#global FILES

CWD = os.getcwd()
SOURCES=["amoco", "bpait", "marmousi", "marmousi2","sigsbee","pluto","wggom"]
#SELECTFILES = ['No files selected','Please Select']                #Global SELECTFILES
class Madagascar:

    def __init__(self, master):
    
###########################
###     Labels          ###
###########################
        Label(master, text="Sources",fg="blue",font=("Helvetica", 14)).grid(row=1,columnspan=3) 
        Label(master, text="Files",fg="blue",font=("Helvetica", 14)).grid(row=1,column=4,columnspan=5)
        Label(master, text="Header Info",fg="blue",font=("Helvetica", 14)).grid(row=20,column=0,columnspan=3)
        Label(master, text="Images and Plots",fg="blue",font=("Helvetica", 14)).grid(row=30,column=0,columnspan=3)

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
        
        self.unit1 = Entry(master,width=5)
        self.unit2 = Entry(master,width=5)
        self.unit3 = Entry(master,width=5)
        
        self.unit1.grid(row=21, column=9)
        self.unit2.grid(row=22, column=9)
        self.unit3.grid(row=23, column=9)

###########################
###     Buttons         ###
###########################
        self.select = Button(master, text="Select File(s)", fg="black", command=self.select)
        self.select.grid(row=3,column=13,sticky=E+W)
        
        self.load = Button(master, text="Load Source", command=self.load)
        self.load.grid(row=2, column=13,sticky=E+W)
        
        self.fetch = Button(master, text="Fetch Files", command=self.fetch)    
        self.fetch.grid(row=4, column=13,sticky=E+W)
        
        self.convert = Button(master, text="Convert Files", command=self.convert)
        self.convert.grid(row=5,column=13,sticky=E+W)

        self.view = Button(master, text="View File Info", command=self.view)
        self.view.grid(row=6,column=13,sticky=E+W)
        
        self.header = Button(master, text="Update Header",command=self.header)
        self.header.grid(row=21, column=13,sticky=E+W)
    
        self.archive = Button(master, text="Archive", command=self.archive)
        self.archive.grid(row=41, column=13, sticky=E+W)
        
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
        if source == 'other':
            importFileWindow = Toplevel()
        lengthCWD = len(CWD)
        lengthCWD = lengthCWD-3
        location = str(CWD[0:lengthCWD]) + str(SOURCES[int(source[0])])+'/FILES'
        global LOCATION                   # Used by function fetch
        LOCATION=str(SOURCES[int(source[0])])
        if LOCATION == "marmousi":
            LOCATION = "marm"
        if LOCATION == "amoco":
            LOCATION = "Amoco"
        if LOCATION == "marmousi2":
            LOCATION = "marm2"
        if LOCATION == "wggom":
            Location = "gom"
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
            extStart2 = length -3
            extStart3 = length -2
            extStart4 = length - 1
            extensionOne = item[extStart:length]
            extensionTwo = item[extStart2:length]
            extensionThree = item[extStart3:length]
            extensionFour = item[extStart4:length]
            if extensionOne == 'segy' or extensionOne == 'SEGY' or extensionTwo == 'SGY' or extensionTwo == 'sgy':
                segyFiles.append(item)
            if extensionThree == 'HH' or extensionThree == 'hh' or extensionFour == 'H' or extensionFour =='h':
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


#--------------------
# Updated Header Info    
#--------------------
    def header(self):
        n1=self.n1.get()
        n2=self.n2.get()
        n3=self.n3.get()
        N=[n1,n2,n3]

        o1=self.o1.get()
        o2=self.o2.get()
        o3=self.o3.get()
        o2=self.o2.get()
        O=[o1,o2,o3]

        d1=self.d1.get()
        d2=self.d2.get()
        d3=self.d3.get()
        D=[d1,d2,d3]

        label1=self.label1.get()
        label2=self.label2.get()
        label3=self.label3.get()
        Labels=[label1,label2,label3]

        unit1=self.unit1.get()
        unit2=self.unit2.get()
        unit3=self.unit3.get()
        Unit=[unit1,unit2,unit3]

        headerOutputWindow = Toplevel()
        title=Label(headerOutputWindow,text='Header Info')
        title.grid(row=0,column=0,columnspan=3)

        ###   N Data Collection
        nRow=1
        counter=0
        for item in ['1','2','3']:
            labelText='n'+item
            Label(headerOutputWindow,text=labelText,fg="black",font=("Helvetica",10,"italic")).grid(row=nRow)
            Label(headerOutputWindow,text=N[counter],fg='red',font=("Helvetica",10,"italic")).grid(row=nRow,column=1)
            nRow = nRow+1
            counter=counter+1
        ###   O Data Collection 
        oRow=1
        counter=0
        for item in ['1','2','3']:
            labelText='o'+item
            Label(headerOutputWindow,text=labelText,fg="black",font=("Helvetica",10,"italic")).grid(row=oRow,column=2)
            Label(headerOutputWindow,text=O[counter],fg='red',font=("Helvetica",10,"italic")).grid(row=oRow,column=3)
            oRow = oRow+1
            counter=counter+1
        ###   D Data Collection 
        dRow=1
        counter=0
        for item in ['1','2','3']:
            labelText='d'+item
            Label(headerOutputWindow,text=labelText,fg="black",font=("Helvetica",10,"italic")).grid(row=dRow,column=4)
            Label(headerOutputWindow,text=D[counter],fg='red',font=("Helvetica",10,"italic")).grid(row=dRow,column=5)
            dRow = dRow+1
            counter=counter+1
        ###   Label Data Collection 
        labelRow=1
        counter=0
        for item in ['1','2','3']:
            labelText='Label '+item
            Label(headerOutputWindow,text=labelText,fg="black",font=("Helvetica",10,"italic")).grid(row=labelRow,column=6)
            Label(headerOutputWindow,text=Labels[counter],fg='red',font=("Helvetica",10,"italic")).grid(row=labelRow,column=7)
            labelRow = labelRow+1
            counter=counter+1
        ###   Unit Data Collection 
        unitRow=1
        counter=0
        for item in ['1','2','3']:
            labelText='Unit'+item
            Label(headerOutputWindow,text=labelText,fg="black",font=("Helvetica",10,"italic")).grid(row=unitRow,column=8)
            Label(headerOutputWindow,text=Unit[counter],fg='red',font=("Helvetica",10,"italic")).grid(row=unitRow,column=9)
            unitRow = unitRow+1
            counter=counter+1

        ###  RSFPut function 
        newHeaderInfo = ''    
        rsfFiles = RSFfiles
        counter = 0
        for file in rsfFiles:
            if o1 is not "":
                newHeaderInfo = " o1="+ o1 
            if o2 is not "":
                newHeaderInfo = newHeaderInfo + " o2="+ o2 
            if o3 is not "":
                newHeaderInfo = newHeaderInfo + " o3="+ o3 
            if d1 is not "":
                newHeaderInfo = newHeaderInfo + " d1="+ d1
            if d2 is not "":
                newHeaderInfo = newHeaderInfo + " d2="+ d2
            if d3 is not "":
                newHeaderInfo = newHeaderInfo + " d3="+ d3
            if n1 is not "":
                newHeaderInfo = newHeaderInfo + " n1="+ n1
            if n2 is not "":
                newHeaderInfo = newHeaderInfo + " n2="+ n2
            if n3 is not "":
                newHeaderInfo = newHeaderInfo + " n3="+ n3
            input = "Flow(\'" +"update_" + file + "\',\'" + file + "\',\'put " + newHeaderInfo + "\')"
            command=LOCATION + "/SConstruct"
            SConstruct=open(command,'a')
            SConstruct.write(input)
            SConstruct.write('\n')
            SConstruct.close()
            newFile = "update_"+file
            RSFfiles.append(newFile)
        command = "cd " + LOCATION + "\n  pwd \n  scons "
        os.system(command)
        print "rsf files"
        print RSFfiles
#--------------------
# Fetch Files   
#--------------------
    def fetch(self):
    #    location = os.system("pwd")
        command = "mkdir " + LOCATION
        os.system(command)
        command2=LOCATION + '/SConstruct' 
        print command2
        SConstruct=open(command2,'w')
        SConstruct.write("from rsfproj import *")
        SConstruct.write("\n")
        for file in SELECTFILES:
            input = "Fetch(\"" + file + "\"," + "\"" + LOCATION + "\")"
            SConstruct.write(input)
            SConstruct.write("\n") 
        print file
        SConstruct.close()
        command = "cd " + LOCATION + " \n  pwd \n  scons "
        os.system(command)
        #os.system("pwd")

    def view(self):
        fileInfo = Toplevel()
        fileNameRow=0
        columnNum = 0
        for file in RSFfiles:
            Label(fileInfo,text=file,fg="black",font=("Helvetica",12)).grid(row=fileNameRow,column=columnNum)
            command = "cd " + LOCATION + "\nsfin " + file + ".rsf > info_" + file
            makeFile = os.system(command)
            command2=LOCATION + "/info_" + file
            readFile = open(command2,'r')
            info = readFile.readlines()
            readFile.close()
            fileInfoRow=1
            for item in info:
                Label(fileInfo,text=item,fg="red",font=("Helvetica",10,"italic")).grid(row=fileInfoRow,column=columnNum)
                fileInfoRow = fileInfoRow + 1
            columnNum=columnNum+1
    
    def convert(self):
        global RSFfiles
        RSFfiles=[]
        command2=LOCATION + '/SConstruct' 
        SConstruct=open(command2,'a')
        SConstruct.write("#Convert Files to RSF\n")
        for file in SELECTFILES:
            for file2 in segyFiles:
                if file is file2:
                    fileOut = ''
                    fileLength =len(file)
                    for letter in file[0:fileLength]:
                        if letter is not '.':
                            fileOut = fileOut + letter
                        if letter == '.':
                            break
                    RSFfiles.append(fileOut)
                    command = "Flow(\'"+ fileOut + "\',\'" + file + "\', 'segyread tape=$SOURCE',stdin=0)"
                    SConstruct.write(command)
                    SConstruct.write('\n')
            for file3 in nativeFiles:
                if file is file3:
                    print file
        command3="cd "+LOCATION+"\n" + "scons"
        os.system(command3)

    def archive(self):
        pass

#########################################################################################################
#########################################################################################################
#########################################################################################################

## Declare Root Frame
## Run App Madagascar on root frame
## Create Title Bar and insert Logo

root = Tk()
title=Label(root,text='Madagascar',height=2,borderwidth=1,font=("Helvetica", 16,"bold"))
title.grid(row=0,column=3,columnspan=3)
#logo=PhotoImage(file="Madagascar2.xbm")
#logo.grid(row=0,column=0,columnspan=3)
app = Madagascar(root)
root.mainloop()
