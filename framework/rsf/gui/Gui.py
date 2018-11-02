from Tkinter import *
from Flow import Flow
from Flow import LinkedFlow
from Sandbox import Sandbox
from Browser import Browser

import tkFileDialog, tkMessageBox, Util, os, random,string, subprocess, pickle, signal

class Gui(Tk):
    '''
    The Gui class is the main class for using TkMadagascar.

    It creates and maintains the relationships between all necessary
    widgets and subwindows.

    The Gui is primarily responsible for file operations:
    --- Loading state information
    --- Saving state information
    --- Exporting to SConstruct

    '''

    def __init__(self):
        Tk.__init__(self)
        self.title('tkMadagascar')

        # Are we showing the browser?
        # showBrowser is the state variable that tells us
        # whether we are showing or not.
        self.showBrowser = IntVar()
        self.showBrowser.set(0)

        self.showLog     = IntVar()
        self.showLog.set(0)

        # The name of a temporary SConstruct file that we are
        # using for running these scripts.
        self.scons = None
        # The process instance that is running, so that we can poll it.
        # Are we running?  We are running if the process is alive.
        self.process = None
        # Name of the log file, and the position in the logfile.
        self.logfile   = None
        self.logpos    = 0L
        self.logWindow = None
        self.log       = None

        # Get the program list once for the Browser window...
        # The browser could do this automatically, but this way the
        # browser could browse any programs that satisfy the formatting
        # requirements.
        self.programs, self.programNames = Util.getPrograms()

        # Of course, create the screen widgets.
        self.__createWidgets()

    def __createWidgets(self):
        '''
        Create all widgets attached to this GUI.
        '''
        self.__createCanvas()
        self.__createMenu()

    def exportFlows(self,name=None):
        '''
        Save all flows on our canvas to a file.

        If name is not provided, then open a File Chooser to get a name
        and save it there.
        '''
        # Get all items on this canvas
        items = self.canvas.getItems()

        if not name:
            # Get a file name using a file chooser
            formats = [('SConstruct','SConstruct*')] # Restrict to SConstruct names
            name = tkFileDialog.asksaveasfilename(title="Export to SConstruct",filetypes=formats)

        # If we found a filename.
        if len(name) > 0:
            try:
                file = open(name,'w')
                file.write('from rsf.proj import *\n\n')

                for item in items:
                    # Write each item out if it was properly configured.
                    # If not, then display a dialog warning the user that
                    # the item will be ignored.
                    command = item.getCommand()
                    if command:
                        file.write(command+'\n')
                    else:
                        if isinstance(item,LinkedFlow):
                            tkMessageBox.showwarning("LinkedFlow not configured",
                            "LinkedFlow %d was not completely configured.  It will not be exported." % (item.id))
                        else:
                            tkMessageBox.showwarning("Flow not configured",
                                "%s %d was not configured.  It will not be exported." % (item.ftype,item.id))
                        continue
                file.write('End()')
                file.close()
            except Exception, e:
                tkMessageBox.showwarning("Exception while exporting",e)

    def save(self):
        '''
        Save the state of all items attached to the current Sandbox (canvas).

        This is distinctly different from exportFlows, in that we save the items
        using Python's pickle package.

        You CANNOT execute objects saved using this method on the command line.

        Objects are saved to *.tkm file format (tkmadagascar) format.
        '''
        formats = [('tkMadagascar object format','*.tkm')]
        # Get a name to save all objects to
        name = tkFileDialog.asksaveasfilename(title="Save state to:",
                                              filetypes=formats)
        # If we got a name
        if len(name) > 0:
            # Then get all items attached to this canvas.
            try:
                items = self.canvas.getItems()
                nitem = len(items)

                file = open(name,'w')
                file.write('nitem=%d\n' % nitem) # How many items are pickled?
                for item in items:
                    pickle.dump(item,file) # Write each item out using pickle
                file.close()

                tkMessageBox.showinfo("Saved state:","Succesfully saved state information to: %s" % name)

            except Exception, e:
                tkMessageBox.showwarning("Failed to save state",e)

    def load(self):
        '''
        Load the state of the saved objects from the selected file.

        Add these objects to a Sandbox.  To prevent race conditions, you may
        not load objects once new objects have been added to a Sandbox.

        Valid file suffixes are *.tkm.
        '''
        # Get a file name to load
        formats = [('tkMadagascar object format','*.tkm')]
        name = tkFileDialog.askopenfilename(title="Load state from:",
                                            filetypes=formats)
        if len(name) > 0:
            try:
                file = open(name,'r')
                # How many items are there to read?
                nitem = int(file.readline().split('=')[1])

                flows = []
                # Load all flows first
                # if we fail here, then nothing is added
                # to the canvas.  We have default state.
                for i in range(nitem):
                    flow = pickle.load(file)
                    flows.append(flow)

                # Now try and add the flows to the canvas.
                for flow in flows:
		            if isinstance(flow,LinkedFlow):
		                self.canvas.loadLinkedFlow(flow)
		            else:
		                self.canvas.loadFlow(flow)
                tkMessageBox.showinfo("Load state:","Successfully restored state, you may proceed.")
            except Exception,e :
                tkMessageBox.showwarning("Failed to load state", e)
                self.canvas.reset()

    def getSConstruct(self):
        '''
        Get a unique name for a temporary SConstruct in the current
        directory.
        '''
        if not self.scons:
            drct = os.getcwd()
            # Grab 10 random characters from all available letters in a string
            # and combine them together without space
            rstr = ''.join([random.choice(string.ascii_letters) for i in range(10)])
            name = 'SConstruct.'+rstr
            # The path is the absolute path name.
            path = os.path.join(drct,name)
            # Save this information.
            log = 'SConstruct.log.'+rstr
            self.scons     = name
            self.logfile   = log
        return self.scons

    def deleteSConstruct(self):
        '''
        Delete the temporary SConstruct that we were using for running
        in the current directory.
        '''
        if self.scons:
            try:
                os.remove(self.scons)
                os.remove(self.logfile)
                tkMessageBox.showinfo("Deleted temporary files:","Removed: %s and %s" % (self.scons,self.logfile))
                self.scons = None
                self.logfile   = None
                self.logpos = 0L
            except Exception, e:
                tkMessageBox.showwarning("Error deleting temporary files:",e)

    def waitForProcess(self):
        '''
        Wait for the current process to finish, while updating the log file.
        '''
        if self.process:
            rc = self.process.poll()

            if self.log:
                file = open(self.logfile,'r')
                file.seek(self.logpos)
                contents = file.read()
                # Remember where we were, so that we only update the changes.
                self.logpos = file.tell()
                self.log.insert(END,contents)
                self.log.yview_moveto(0.99)
                file.close()

            if not rc and rc != 0: # Process still running
                self.after(500,self.waitForProcess )

            elif rc == 0:
                self.log.insert(END,"\nSUCCESSFUL RUN\n")
                self.process = None

            elif rc > 0:
                self.log.insert(END,"\nRUN TERMINATED DUE TO ERRORS\n")
                self.process = None
            elif rc < 0:
                self.log.insert(END,"\nRUN KILLED BY SIGNAL %d\n" % rc)
                self.process = None

    def killProcess(self):

        if self.process:
            if self.log:
                self.log.insert(END,"****\nKILLING PROCESS: %d\n****\n" % self.process.pid)

            rc = self.process.poll()
            if rc == None:
                self.process.terminate()# Sigint?

            rc = self.process.poll()
            if rc == None:
                self.process.kill()# Sigint?
            self.after(500,self.waitForProcess)

    def __destroyLogWindow(self):
        Toplevel.destroy(self.logWindow)
        self.logWindow = None
        self.log       = None
        self.showLog.set(0)
        self.viewmenu.entryconfigure(1,state=NORMAL)


    def __showLogWindow(self):
        if not self.logWindow:
            self.showLog.set(1)
            self.viewmenu.entryconfigure(1,state=DISABLED)
            self.logWindow = Toplevel()
            self.logWindow.title("Log")
            scroll = Scrollbar(self.logWindow)
            self.log = Text(self.logWindow,yscrollcommand=scroll.set)
            menubar = Menu(self.logWindow,tearoff=0)
            menubar.add_command(label="Kill process",command=self.killProcess)

            self.logWindow.config(menu=menubar)
            scroll['command'] = self.log.yview
            scroll.pack(fill=Y,side=RIGHT)
            self.log.pack(expand=1, fill=Y,side=LEFT)
            self.logWindow.protocol("WM_DELETE_WINDOW",self.__destroyLogWindow)

    def runCommand(self,command,delete=False):
        '''
        Execute the given command on the command line using the subprocess
        module.

        This command will save all flows attached to the current Sandbox
        to a temporary SConstruct, and then execute the command on that
        SConstruct.

        If delete, then remove the temporary SConstruct file.
        '''
        if not self.process: # DO NOT RUN if another process is running.

            self.__showLogWindow()

            name = self.getSConstruct()

            command += ' -f %s' % name

            self.exportFlows(name)
            self.log.insert(END,
'''
Saving flows to: %s
Will execute: %s
''' % (name, command))

            logfile = open(self.logfile,'a')
            self.process = subprocess.Popen(command, shell=True,stdout=logfile,stderr=subprocess.STDOUT)

            self.log.insert(END,"PROCESS ID: %d\n" % self.process.pid)
            self.after(500,self.waitForProcess)
            logfile.close()

    def __createMenu(self):
        '''
        Create the menubar for this GUI
        '''

        menubar = Menu(self)

        self.filemenu = Menu(menubar,tearoff=0)
        self.filemenu.add_command(label="New",command=self.canvas.reset)
        self.filemenu.add_separator()
        self.filemenu.add_command(label="Load state",command=self.load)
        self.filemenu.add_command(label="Save state",command=self.save)
        self.filemenu.add_separator()
        self.filemenu.add_command(label="Export to SConstruct",
            command=self.exportFlows)
        self.filemenu.add_separator()
        self.filemenu.add_command(label="Quit",command=self.destroy)

        menubar.add_cascade(label="File",menu=self.filemenu)

        self.viewmenu = Menu(menubar,tearoff=0)
        self.viewmenu.add_checkbutton(label="Program browser",
                variable=self.showBrowser,
                command=self.__showProgramBrowser)
        self.viewmenu.add_checkbutton(label="Log",
                variable=self.showLog,
                command=self.__showLogWindow)
        menubar.add_cascade(label="View",menu=self.viewmenu)

        self.runmenu = Menu(menubar,tearoff=0)

        self.runmenu.add_command(label="scons",
            command=lambda arg='scons': self.runCommand(arg))

        self.runmenu.add_command(label="scons view",
            command=lambda arg='scons view': self.runCommand(arg))

        self.runmenu.add_command(label="scons lock",
            command=lambda arg='scons lock': self.runCommand(arg))

        self.runmenu.add_command(label="scons -c",
            command=lambda arg='scons -c': self.runCommand(arg,delete=True))

        self.runmenu.add_command(label="scons -n",
            command=lambda arg='scons -n': self.runCommand(arg))

        menubar.add_cascade(label="Run...",menu=self.runmenu)
        self.config(menu=menubar)


    def __createCanvas(self):
        '''
        Create and save the Sandbox associated with this GUI.
        '''

        self.canvasFrame = Frame(self)
        self.canvasFrame.pack(fill=BOTH,expand=1)

        self.canvas = Sandbox(self.canvasFrame)

        self.csby = Scrollbar(self.canvasFrame,orient=VERTICAL,command=self.canvas.yview)
        self.csby.pack(side=RIGHT,fill=Y)

        self.canvas.pack(fill=BOTH,expand=1)
        self.csbx = Scrollbar(self.canvasFrame,orient=HORIZONTAL,command=self.canvas.xview)
        self.csbx.pack(side=BOTTOM,fill=X)

        self.canvas['xscrollcommand'] = self.csbx.set
        self.canvas['yscrollcommand'] = self.csby.set


    def __showProgramBrowser(self):
        '''
        Show the Program Browser in a separate window.
        '''
        if self.showBrowser.get() == 1:
            self.viewmenu.entryconfigure(0,state=DISABLED)
            top = Toplevel()
            top.title('Program browser - tkMadagascar')
            top.resizable(False,False)
            def closeBrowser():
                self.showBrowser.set(0)
                self.viewmenu.entryconfigure(0,state=NORMAL)
                top.destroy()

            top.protocol("WM_DELETE_WINDOW",closeBrowser)
            browserFrame = Browser(top,self.programs,self.programNames,self.canvas)
            browserFrame.pack()

    def destroy(self):
        if self.process:
            rc = self.process.poll()
            if rc == None:
                tkMessageBox.showwarning("Process %d still active:" % self.process.pid, "tkMadagascar still has a child process active.  You may not quit until this process is either terminated or finishes.")
        else:
            if tkMessageBox.askyesno("Close tkMadagascar?", "Are you sure you want to quit?"):
                self.deleteSConstruct()
                Tk.destroy(self)

def launch():
    gui = Gui()
    gui.mainloop()

if __name__=='__main__': # Execute the main program
    launch()
