from Tkinter import *
from Flow import Flow
import Util,tkMessageBox

class Browser(Frame):
    '''
    The Browser is a Frame that contains all the searching, buttons and other tools 
    used to quickly browse and search for Madagascar programs.
    
    If it is attached to a Canvas, then it can AddFlows directly to that
    Canvas.
    
    Otherwise, it can standAlone, meaning that it will exclude certain features, 
    such as Buttons to add items to a Canvas.
    
    To run in standalone mode directly call this program from the command line.
    
    '''

    def __init__(self,master,programs,programNames,canvas=None,standAlone=False):

        Frame.__init__(self,master)

        self.parent = master 
        self.programs = programs
        self.programNames = programNames

        self.standAlone = standAlone # Am I running as my own program?
        self.canvas     = canvas

        menubar = Menu(self)

        menubar.add_command(label="Search by functionality",command=self.__searchDescriptions)
        master.config(menu=menubar)

        programFrame = Frame(self,bd=3)
        programFrame.pack(side=RIGHT,fill=Y,expand=1)

        lbscrolly = Scrollbar(programFrame,orient=VERTICAL)
        lbscrolly.pack(side=RIGHT,fill=Y)

        searchFrame = Frame(programFrame)
        searchFrame.pack(side=TOP)
        self.programEntry = Entry(searchFrame)
        self.programEntry.bind('<Return>',self.__searchProgramList)

        entryButton = Button(searchFrame,text="Find program")
        entryButton.bind('<Button-1>',self.__searchProgramList)

        self.programEntry.pack(side=LEFT)
        entryButton.pack(side=RIGHT)

        self.programList = Listbox(programFrame,height=18,
            yscrollcommand=lbscrolly.set)

        lbscrolly['command'] = self.programList.yview

        self.programList.bind('<Double-Button-1>',self.__updateDocumentation)

        if not standAlone:
            buttonFrame = Frame(programFrame)
            buttonFrame.pack(side=BOTTOM,fill=X,expand=1)

            addButton = Button(buttonFrame,text="Add as Flow")
            plotButton = Button(buttonFrame,text="Add as Plot")
            resultButton = Button(buttonFrame,text="Add as Result")

            addButton.bind('<Button-1>',
                lambda event, arg=Flow.FLOW: self.createFlow(event,arg) )
            plotButton.bind('<Button-1>',
                lambda event, arg=Flow.PLOT: self.createFlow(event,arg) )
            resultButton.bind('<Button-1>',
                lambda event, arg=Flow.RESULT: self.createFlow(event,arg) )
            addButton.pack(side=LEFT)
            plotButton.pack(side=LEFT)
            resultButton.pack(side=RIGHT)

        self.programList.pack(side=LEFT,fill=BOTH,expand=1)
        
        docFrame = Frame(self)
        docFrame.pack(side=LEFT,fill=BOTH,expand=1)

        dscrolly = Scrollbar(docFrame,orient=VERTICAL)
        dscrolly.pack(side=RIGHT,fill=Y)
        self.documentation = Text(docFrame,width=80,state=DISABLED,yscrollcommand=dscrolly.set)
        dscrolly['command'] = self.documentation.yview
        self.documentation.pack(side=LEFT,fill=Y,expand=1)

        for name in self.programNames:
            self.programList.insert(END,name)

    def createFlow(self,event,ftype):
        ''' 
        Create and attach a Flow to a Canvas, if we have one.
        Otherwise, you have called this method erroneously.
        '''
        name = self.__getListItem()
        if not Util.checkNameType(name,ftype):
            tkMessageBox.showwarning("Bad flow type:", "%s may only be a Plot or a Result." % name)
        else:
            if self.canvas:
                self.canvas.addFlow(self.programs[name],ftype)
            else:
                raise Exception("No canvas to add Flows to")

    def __updateDocumentation(self,event):
        '''
        Update the Text widget to display the selected
        programs self-documentation.
        '''
        name = self.__getListItem()
        selfdoc = self.programs[name].selfdoc()
        self.documentation.configure(state=NORMAL)
        self.documentation.delete(1.0,END)
        self.documentation.insert(END,selfdoc)
        self.documentation.configure(state=DISABLED)

    def __getListItem(self):
        '''
        Get the name of the program that the user has
        selected.
        '''
        item = int(self.programList.curselection()[0])
        name = self.programNames[item]
        return name

    def __changeListSelection(self,text):
        '''
        Change the List Selection to match the closest item with 
        the text in it.
        
        This is a partial match move.
        '''
        closest = None
        closestIndex = -1
        for name in self.programNames:
            if text in name:
                closest = name
                break

        if closest:
            closestIndex = self.programNames.index(name)
            self.programList.selection_clear(0,END)
            self.programList.selection_set(closestIndex)
            self.programList.see(closestIndex)
            self.__updateDocumentation(None)
            return True
        else:
            return False

    def __searchProgramList(self,event):
        '''
        Get the partial name of the program that the user is looking
        for, and then move the listbox to the closest matching entry.
        '''
        text = self.programEntry.get()

        if not self.__changeListSelection(text):
            tkMessageBox.showwarning(
            "Bad program name",
            "No program in Madagascar with the name: %s" % text)

    def __searchDescriptions(self):
        '''
        Create a window to allow the user to look for items with
        words that match the description.
        
        This is basically, keyword search.
        '''
        window = Toplevel()

        window.title("Search by functionality")
        window.resizable(False,False)
        
        frame = Frame(window) 
        frame.pack(side=TOP,expand=1,fill=X)

        entry = Entry(frame)
        entry.pack(side=LEFT,expand=1,fill=X)


        button = Button(frame,text="Search")
        button.pack(side=RIGHT)

        tframe = Frame(window)
        tframe.pack(side=BOTTOM)

        scrollbar = Scrollbar(tframe,orient=VERTICAL)
        scrollbar.pack(side=RIGHT,fill=Y)

        listbox = Listbox(tframe,width=40,height=20,yscrollcommand=scrollbar.set)

        scrollbar['command'] = listbox.yview
        listbox.pack(expand=1,fill=Y,side=TOP)

        button2 = Button(tframe,text="Select")
        button2.pack(side=BOTTOM,expand=1,fill=X)

        def __getAndSearch(event):

            text = entry.get().lower()

            hits = []
            hitDescs = []
            for name in self.programNames:
                desc = self.programs[name].desc.lower()

                if text in desc:
                    hits.append(name)
                    hitDescs.append(desc)
            listbox.delete(0,END)
            for i in range(len(hits)):
                listbox.insert(END,'%s - %s' % (hits[i],hitDescs[i]))

        entry.bind('<Return>',__getAndSearch)

        button.bind('<Button-1>',__getAndSearch)
        
        def __seeMoreInfo(event):
            
            item = listbox.get(ACTIVE)
            
            name = item.split('-')[0].replace(' ','')

            self.__changeListSelection(name)
            
        listbox.bind('<Double-Button-1>',__seeMoreInfo) 
        button2.bind('<Button-1>',__seeMoreInfo)

def launch():
    ''' Create a standalone Browser '''
    programs, names = Util.getPrograms()
    tk = Tk()
    b = Browser(tk,programs,names,standAlone=True)
    b.pack(fill=BOTH,expand=1)
    tk.title('Browse Madagascar programs')
    tk.resizable(False,False)
    tk.mainloop()    
    
if __name__ == '__main__':
    launch()



