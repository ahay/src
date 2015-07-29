from Tkinter import *
import os
CWD = os.getcwd()
SOURCES=["amoco", "bpait", "marmousi", "marmousi2","sigsbee","pluto","wggom"]


class Madagascar:

    def __init__(self, master):
  
############################ 
# ######################## #
# #      LABELS          # #
# ######################## #
############################

        Label(master, text="Sources",fg="blue",font=("Helvetica", 14)).grid(row=1,columnspan=3) 
        Label(master, text="Files",fg="blue",font=("Helvetica", 14)).grid(row=1,column=4,columnspan=5)
        Label(master, text="Header Info",fg="blue",font=("Helvetica", 14)).grid(row=20,column=0,columnspan=3)
#        Label(master, text="Processing Operations",fg='blue',font=('Helvetica',14)).grid(row=28,column=0,columnspan=4)
        Label(master, text="Images and Plots",fg="blue",font=("Helvetica", 14)).grid(row=30,column=0,columnspan=3)

###############################
## SF PUT FUNCTION OPTIONS   ##   (HEADER INFO)
###############################
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

###############################
## SFWINDOW FUNCTION OPTIONS ##  (WINDOW)
###############################        
        Label(master, text="Window Data",fg="blue",font=("Helvetica",12,'bold')).grid(row=31,column=4,columnspan=3) 

        Label(master, text="Sampling",fg="black",font=("Helvetica",10,"underline",'bold')).grid(row=32,column=1,columnspan=1)
        Label(master, text="Start",fg="black",font=("Helvetica",10,"underline",'bold')).grid(row=32,column=3,columnspan=1)
        Label(master, text="Jump",fg="black",font=("Helvetica",10,"underline",'bold')).grid(row=32,column=5,columnspan=1) 
        Label(master, text="Min",fg="black",font=("Helvetica",10,"underline",'bold')).grid(row=32,column=6,columnspan=1)
        Label(master, text="Max",fg="black",font=("Helvetica",10,"underline",'bold')).grid(row=32,column=7,columnspan=1)
        Label(master, text="Size",fg="black",font=("Helvetica",10,"underline",'bold')).grid(row=32,column=8,columnspan=1)

        Label(master, text="time:").grid(row=33,column=0,sticky=W)
        Label(master, text="offset:").grid(row=34,column=0,sticky=W)
        Label(master, text="shot:").grid(row=35,column=0,sticky=W)            

###############################
## SF GREY FUNCTION OPTIONS  ##   (RASTOR PLOT)
###############################      
        Label(master, text="Raster Plot",fg="blue",font=("Helvetica",12,'bold')).grid(row=39,column=4,columnspan=3)
        Label(master, text="Color Scheme").grid(row=40,column=0,sticky=W)
        Label(master, text="Gainpanel").grid(row=40,column=2)
        Label(master, text="Clip (%)").grid(row=40,column=4)
        
         
############################ 
# ######################## #
# #    ENTRY CELLS       # #
# ######################## #
############################

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
        

##########################
###   WINDOW ENTRY     ### 
########################## 
        self.samp1 = Entry(master,width=5)
        self.samp2 = Entry(master,width=5)
        self.samp3 = Entry(master,width=5)
        
        self.samp1.grid(row=33, column=1)
        self.samp2.grid(row=34, column=1)
        self.samp3.grid(row=35, column=1)

        self.start1 = Entry(master,width=5)
        self.start2 = Entry(master,width=5)
        self.start3 = Entry(master,width=5)
        
        self.start1.grid(row=33, column=3)
        self.start2.grid(row=34, column=3)
        self.start3.grid(row=35, column=3)
        
        self.jump1 = Entry(master,width=5)
        self.jump2 = Entry(master,width=5)
        self.jump3 = Entry(master,width=5)
        
        self.jump1.grid(row=33, column=5)
        self.jump2.grid(row=34, column=5)
        self.jump3.grid(row=35, column=5)
        
        self.min1 = Entry(master,width=5)
        self.min2 = Entry(master,width=5)
        self.min3 = Entry(master,width=5)
        
        self.min1.grid(row=33, column=6)
        self.min2.grid(row=34, column=6)
        self.min3.grid(row=35, column=6)
   
        self.max1 = Entry(master,width=5)
        self.max2 = Entry(master,width=5)
        self.max3 = Entry(master,width=5)
        
        self.max1.grid(row=33, column=7)
        self.max2.grid(row=34, column=7)
        self.max3.grid(row=35, column=7)

        self.size1 = Entry(master,width=5)
        self.size2 = Entry(master,width=5)
        self.size3 = Entry(master,width=5)
        
        self.size1.grid(row=33, column=8)
        self.size2.grid(row=34, column=8)
        self.size3.grid(row=35, column=8)

        self.clip = Entry(master,width=3)
        
        self.clip.grid(row=40, column=5)

############################ 
# ######################## #
# #      BUTTONS         # #
# ######################## #
############################

        self.select = Button(master, text="Select File(s)", fg="black", command=self.select)
        self.select.grid(row=3,column=13,sticky=E+W)
        
        self.load = Button(master, text="Load Source", command=self.load)
        self.load.grid(row=2, column=13,sticky=E+W)
        
        self.fetch = Button(master, text="Fetch Files", command=self.fetch)    
        self.fetch.grid(row=4, column=13,sticky=E+W)
        
        self.convert = Button(master, text="Convert Files", command=self.convert)
        self.convert.grid(row=5,column=13,sticky=E+W)
        
        self.convert = Button(master, text="Concatinate Files", command=self.loadCat)
        self.convert.grid(row=23,column=13,sticky=E+W)

        self.view = Button(master, text="View File Info", command=self.view)
        self.view.grid(row=21,column=13,sticky=E+W)
        
        self.header = Button(master, text="Update Header",command=self.header)
        self.header.grid(row=22, column=13,sticky=E+W)
    
        self.archive = Button(master, text="Archive", command=self.archive)
        self.archive.grid(row=51, column=13, sticky=E+W)
    
        self.archive = Button(master, text="Plot", command=self.plot)
        self.archive.grid(row=41, column=13, sticky=E+W)
       
 
############################ 
# ######################## #
# #     LIST BOXES       # #
# ######################## #
############################

############################
### Source List Box      ###
############################
        self.sources = Listbox(master,selectmode=SINGLE,relief=RAISED,height=10,exportselection=0)
        for item in SOURCES:
            self.sources.insert(END, item)
        self.sources.grid(row=2,column=0,rowspan=5,columnspan=3)

############################
### Avail. File List Box ###
############################
        self.fileScrollBar=Scrollbar(master,orient=VERTICAL)
        self.files = Listbox(master,yscrollcommand=self.fileScrollBar.set,selectmode=MULTIPLE,
                              relief=RAISED,height=10,width=40,exportselection=0)
        for item in ['To view files','load a source']:
            self.files.insert(END, item)
      #  self.fileScrollBar.config(command=self.files.yview)
        self.fileScrollBar.grid(column=9,row=2,rowspan=5,sticky=N+S)
        self.files.grid(column=4,row=2,rowspan=5,columnspan=5)
        self.current = None

#########################
### Color List Box    ###
#########################
        global colorOpts
        colorOpts = ['i','I','J','K','F','R','W','G','T' ]
    #    self.colorScrollBar=Scrollbar(master,orient=VERTICAL)
        self.color = Listbox(master,selectmode=SINGLE,relief=RAISED,height=3,width=2,exportselection=0)
        for item in colorOpts:
            self.color.insert(END, item)
    #    self.colorScrollBar.config(command=self.color.yview)
    #    self.colorScrollBar.grid(column=2,row=40,rowspan=5,sticky=N+S)
        self.color.grid(column=1,row=40,rowspan=1,columnspan=1)
        self.current = None

#########################
### Gain Panel List   ###
#########################
        global gainpanelOpts
        gainpanelOpts = ['a','e']
        self.gainpanel = Listbox(master,selectmode=SINGLE,relief=RAISED,height=2,width=2,exportselection=0)
        for item in gainpanelOpts:
            self.gainpanel.insert(END, item)
        self.gainpanel.grid(column=3,row=40,rowspan=1,columnspan=1)
        self.gainpanel.current = None


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
            LOCATION = "gom"
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

        ###  RSFPut function 
        command=LOCATION + '/SConstruct'
        input=open(command,'r')
        OldSConstruct = input.readlines()
        input.close()
        SConstruct = open(command,'w')
        sectionHead = '# Update Header \n'
        for item in OldSConstruct:
            if item != sectionHead:
                SConstruct.write(item)
            if item == sectionHead:
                break
        SConstruct.write(sectionHead)
        counter = 0
        newFiles = []
        for file in RSFfiles:
            newHeaderInfo = ''
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
            if label1 is not "":
                newHeaderInfo = newHeaderInfo + " label1=" + label1
            if label2 is not "":
                newHeaderInfo = newHeaderInfo + " label2=" + label2
            if label3 is not "":
                newHeaderInfo = newHeaderInfo + " label3=" + label3
            if unit1 is not "":
                newHeaderInfo = newHeaderInfo + " unit1=" + unit1
            if unit2 is not "":
                newHeaderInfo = newHeaderInfo + " unit2=" + unit2
            if unit3 is not "":
                newHeaderInfo = newHeaderInfo + " unit3=" + unit3
            input = "Flow(\'" +"update_" + file + "\',\'" + file + "\',\'put " + newHeaderInfo + "\')"
            command=LOCATION + "/SConstruct"
            SConstruct=open(command,'a')
            SConstruct.write(input)
            SConstruct.write('\n')
            SConstruct.close()
            newFile = "update_"+file
            newFiles.append(newFile)
        command = "cd " + LOCATION + "\n  pwd \n  scons &"
        os.system(command)
        print 'newFiles'
        print newFiles
        length = len(RSFfiles)
        for file in RSFfiles:
            trash = RSFfiles.pop()
        if length > 1:
            trash = RSFfiles.pop()
        print 'midRSF'
        print RSFfiles
        for file in newFiles:
            RSFfiles.append(file)
        print 'RSF'
        print RSFfiles
#--------------------
# Fetch Files   
#--------------------
    def fetch(self):
    #    location = os.system("pwd")
        command = "mkdir " + LOCATION
        os.system(command)
        command2=LOCATION + '/SConstruct' 
        SConstruct=open(command2,'w')
        SConstruct.write("from rsf.proj import *")
        SConstruct.write("\n")
        sectionHeader = '# Fetch Files from repository \n'
        SConstruct.write(sectionHeader)
        for file in SELECTFILES:
            input = "Fetch(\"" + file + "\"," + "\"" + LOCATION + "\")"
            SConstruct.write(input)
            SConstruct.write("\n") 
        SConstruct.close()
        command = "cd " + LOCATION + " \n  pwd \n  scons &"
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
        SConstruct.write("# Convert Files to RSF\n")
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
                    fileOut = ''
                    fileLength = len(file)
                    for letter in file[0:fileLength]:
                        if letter is not '.':
                            fileOut = fileOut + letter
                        if letter == '.':
                            break
                    RSFfiles.append(fileOut)
                    command = "Flow(\'"+ fileOut + "\',\'" + file + "\', 'dd form=native')"
                    SConstruct.write(command)
                    SConstruct.write('\n')
                    SConstruct.write('# Update Header \n')
        SConstruct.close()
        command3="cd "+LOCATION+"\n" + "scons &"
        os.system(command3)

    def plot(self):
        command=LOCATION + '/SConstruct'
        input=open(command,'r')
        OldSConstruct = input.readlines()
        input.close()
        SConstruct = open(command,'w')
        sectionHeader = '# Plotting Section\n'
        for item in OldSConstruct:
            if item != sectionHeader:
                SConstruct.write(item)
            if str(item) == sectionHeader:
                break 
        SConstruct.write(sectionHeader)
 
##
## Get Values to Plug into SFwindow and SFgrey
##
        samp1=self.samp1.get()
        samp2=self.samp2.get()
        samp3=self.samp3.get()

        start1=self.start1.get()
        start2=self.start2.get()
        start3=self.start3.get()
        
        jump1=self.jump1.get()
        jump2=self.jump2.get()
        jump3=self.jump3.get()

        min1=self.min1.get()
        min2=self.min2.get()
        min3=self.min3.get()

        max1=self.max1.get()
        max2=self.max2.get()
        max3=self.max3.get()

        size1=self.size1.get()
        size2=self.size2.get()
        size3=self.size3.get()

        clip=self.clip.get()

        if self.color.curselection() is not ():
            colorTuple = self.color.curselection()
            selection = int(colorTuple[0])
            color = colorOpts[selection]
        else:
            color = 'I'

        if self.gainpanel.curselection() is not ():
            gainTuple = self.gainpanel.curselection()
            choice = int(gainTuple[0])
            gainpanel = gainpanelOpts[choice]
        else:
            gainpanel = 'a'

        for file in RSFfiles:
            windowInfo = ''
            if samp1 is not '':
                windowInfo = windowInfo + ' d1='+ samp1
            if samp2 is not '':
                windowInfo = windowInfo + ' d2='+ samp2
            if samp3 is not '':
                windowInfo = windowInfo + ' d3='+ samp3
            if start1 is not '':
                windowInfo = windowInfo + ' f1='+ start1
            if start2 is not '':
                windowInfo = windowInfo + ' f2='+ start2
            if start3 is not '':
                windowInfo = windowInfo + ' f3='+ start3
            if jump1 is not '':
                windowInfo = windowInfo + ' j1='+ jump1
            if jump2 is not '':
                windowInfo = windowInfo + ' j2='+ jump2
            if jump3 is not '':
                windowInfo = windowInfo + ' j3='+ jump3
            if min1 is not '':
                windowInfo = windowInfo + ' min1='+ min1
            if min2 is not '':
                windowInfo = windowInfo + ' min2='+ min2
            if min3 is not '':
                windowInfo = windowInfo + ' min3='+ min3
            if max1 is not '':
                windowInfo = windowInfo + ' max1='+ max1
            if max2 is not '':
                windowInfo = windowInfo + ' max2='+ max2
            if max3 is not '':
                windowInfo = windowInfo + ' max3='+ max3
            if size1 is not '':
                windowInfo = windowInfo + ' size1='+ size1
            if size2 is not '':
                windowInfo = windowInfo + ' size2='+ size2
            if size3 is not '':
                windowInfo = windowInfo + ' size3='+ size3

            greyInfo=''
            if clip is not '':
                greyInfo = greyInfo + ' pclip=' + clip
            if color is not '':
                greyInfo = greyInfo + ' color=' + color
            if gainpanel is not '':
                greyInfo = greyInfo + ' gainpanel=' + gainpanel            
            
            windowCommand = '\'window $SOURCE ' + windowInfo            
            command = "Result(\'"+ file +"\',"+ windowCommand + ' | grey ' + greyInfo  + '\')'
            
            SConstruct.write(command)
            SConstruct.write('\n')
        
        SConstruct.write('End()')
        SConstruct.close()    
        command3="cd "+LOCATION+"\n" + "scons view &"
        os.system(command3)

    def loadCat(self): 
        Label(root, text="Select Files To Combine").grid(row=25,column=0,columnspan=2,pady=2) 
        Label(root, text="Cat Axis").grid(row=25,column=7,columnspan=1,pady=2)  
        self.catFiles = Listbox(root,selectmode=MULTIPLE,relief=RAISED,height=3,width=35,exportselection=0)
        for item in RSFfiles:
            self.catFiles.insert(END, item)
        self.catFiles.grid(row=25,column=2,rowspan=3,columnspan=5)
        self.catAxis = Entry(root,width=3)
        self.catAxis.grid(row=25, column=8,pady=2) 
        self.concatinate = Button(root, text="Combine", command=self.cat)
        self.concatinate.grid(row=25,column=13,sticky=E+W)

    def cat(self):
        choices = self.catFiles.curselection()
        axisNum=self.catAxis.get()
        CATfiles = []
        for item in choices:
            number = int(item)
            selectFiles =  RSFfiles[number]
            CATfiles.append(selectFiles)
        command=LOCATION + '/SConstruct'
        input=open(command,'r')
        OldSConstruct = input.readlines()
        input.close()
        SConstruct = open(command,'w')
        plotHeader = '# Plotting Section\n'
        sectionHeader = '# Concatinate Files \n'
        for item in OldSConstruct:
            if item != plotHeader:
                SConstruct.write(item)
            if str(item) == plotHeader:
                break 
        SConstruct.write(sectionHeader)    
        catRules ='\'cat ${SOURCES[0:2]} axis=' + axisNum + '\''
        command2 = 'Flow(\'catFile_'+ LOCATION + '\',' + str(CATfiles) + ',' + catRules +  ',stdin=0) \n'
        SConstruct.write(command2)
        command3="cd "+LOCATION+"\n" + "scons &"
        os.system(command3) 
        for file in RSFfiles:
            trash = RSFfiles.pop()  
        trash = RSFfiles.pop()            # Remove last entry 
        newFile='catFile_' + LOCATION 
        RSFfiles.append(newFile)
        print RSFfiles
        
    def archive(self):
        command = 'cd ' + LOCATION + '\n' '\\rm update* info*'
        print command
        os.system(command)
       

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
