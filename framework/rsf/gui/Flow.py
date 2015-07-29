from Tkinter import *
import tkMessageBox
import Program

class Flow():
    '''
    A Flow is the workhorse of Madagascar and SCons.  

    Basically, a Flow contains information on how to build an RSF file from other 
    RSF files using the given commands with certain parameters.

    Flows are uninitialized when created, but once set, they "know" how to build
    targets, and are ready to go into an SConstruct.

    Flows can create windows to get parameter values.

    Flows can read and write themselves to files for saving, and can be exported
    to SConstructs for execution.

    Flows also can draw themselves onto a tkCanvas and have fields that help
    with keeping track of themselves on tkCanvases.  

    There are three (3) types of Flows:
    -> Flow
    -> Plot
    -> Result
    The types are functionally identical.  They are separated only to
    allow us to give graphical reinforcement that they are different from 
    one another.
    '''

    # These are "enumerated" types of Flows
    FLOW = 'Flow'
    PLOT = 'Plot'
    RESULT = 'Result'

    # Some additional class wide parameters
    WIDTH    = 100
    HEIGHT   = 50
    BACKGROUND_COLOR       = 'white'

    def __init__(self,id,program,ftype='Flow'):
        '''
        Create a new Flow.
        
         - program is an instance of Program, which tells this flow, 
           which parameters it has and it's program name.
         - id is a number that uniquely identifies this flow.
         - ftype is the type of the Flow (see enumerated types in Flow class)
        '''
        
        # Save identifying variables
        self.ftype = ftype
        self.tag  = 'flow-%d' % id 
        self.id   = id
        self.name = program.name

        # Input and Output represent other Flows that we are piping to or from
        self.input  = None
        self.output = None
        self.linkedFlow = None

        # Copy the programs parameter dictionary, making 
        # deep copies of each parameter, to ensure that 
        # we can have multiple copies of the same program
        # that have independent parameters.
        self.pars = {}
        for key,par in program.pars.items():
            self.pars[key] = par.copy()
            
        keys = self.pars.keys()
        keys.sort()
        keys.remove('output')
        keys.remove('input')
        keys.remove('custom')
        # Put input, output and custom at specific locations in the list
        keys.insert(0,'input')
        keys.insert(0,'output')
        keys.insert(len(keys),'custom')

        self.keys = keys
        # Set tells us when this Flow has been configured by the user.  
        # Until Set, this Flow may not return a command string, and subsequently
        # will not be output to an SConstruct.
        self.set = False

        # Rectangle and text are unique object IDs for representing this Flow on 
        # a Canvas,
        self.rectangle = None
        self.text      = None
        self.outtext   = None
        # Rectangles and text objects that display the unique ID assigned to this flow
        self.rectid    = None
        self.textid    = None

        self.coords    = None

        # Set some Flow, Plot, Result specific colors.
        if self.ftype == Flow.FLOW:
            self.OUTLINE_COLOR='black'
            self.LINKED_COLOR = 'grey'
        elif self.ftype == Flow.PLOT:
            self.OUTLINE_COLOR= 'purple'
            self.LINKED_COLOR = 'violet'
        else:  # Result
            self.OUTLINE_COLOR= 'blue'
            self.LINKED_COLOR = 'sky blue'

        self.canvas = None # Which canvas am I on?
        self.selfdoc = program.selfdoc()

    def __getstate__(self):
        ''' 
        Save this Flow using pickle.
        '''
        if self.canvas:
            self.coords = self.canvas.coords(self.rectangle)
        else:
            self.coords = None

        state = (self.name,self.pars,
            self.keys,self.ftype,
            self.coords,self.OUTLINE_COLOR,self.LINKED_COLOR, self.set,self.selfdoc)

        return state

    def __setstate__(self,state):
        '''
        Restore this Flow using pickle.
        '''
        (self.name,
        self.pars,self.keys,
        self.ftype,
        self.coords,self.OUTLINE_COLOR,self.LINKED_COLOR,self.set,self.selfdoc) = state
        # Set these to defaults
        # If we are a part of a LinkedFlow
        # then the LinkedFlow will re-establish the proper links
        self.output = None
        self.input = None
        self.linkedFlow = None
        self.outtext = None
        self.tag = None
        self.id  = -1
        
    def setTag(self,id):
        self.tag = 'flow-%d' % id
        self.id  = id

    def getCommand(self):
        '''
        Return the string that would be put in an SConstruct for this Flow.
        This string must be properly formatted for execution in an SConstruct.

        If this Flow is piped on either input or output, this method returns the
        command segment corresponding to this Flow.  It must be combined with the
        command segment from its linked flows to make a full command.

        If this Flow has not been configured, then return None indicating that this
        Flow is not to be used.
        '''

        if not self.set:
            return None
        
        command = ''
        command += self.name + ' '

        for par in self.keys: # For all parameters
            if par == 'output' or par == 'input' or par == 'custom':
                # Skip these keywords
                continue
            else:
                # Get the value, and make sure it was set.
                # If it wasn't set, then assume we are using defaults.
                value = self.pars[par].value
                if (value != None) and (value != ''):
                    if '#' in self.pars[par].name: 
                        command += value+' '
                    else:
                        command += '%s=%s ' % (par,value)

        # Custom is just a string of custom commands, no key=val notation
        custom = self.pars['custom'].value
        if custom and custom != '':
            command += custom
        
        # Now we have to do some tricky logic.
        # Are we being Piped?
        # If so, then we only want to return this Flows contribution to the
        # overall command.  
        # If not, then we need to return the full command as would go in
        # an SConstruct.
        header = ""
        if self.output: # I am piped to another program
            command += ' | \n        '
        elif self.name in Program.plotCommands:
            command = command.strip()+"')"
        else: # I am the last in the chain
            command += "\n        ''')"
            
            
        if not self.input:
            if self.output: 
                return command 
                
        # I am not piping from someone else, so go ahead and add
        # the rest of the features
            output   = self.pars['output'].value
            if output and output != '':
                header   = self.ftype+"('%s'," % output
            else:
                self.unconfigure()
                return None

            input = self.pars['input'].value
              
            if input and input != '':
                if self.name in Program.plotCommands:
                    command = header+"'%s','" % input + command
                else:
                    command = header+"'%s',\n        '''\n        " % input + command
            elif self.ftype == Flow.PLOT or self.ftype == Flow.RESULT:
                command = header+"\n        '''\n        " + command
            else:
                command = header+"None,\n        '''\n        " + command
       
        return command

    def delete(self):
        '''
        Remove this Flow from the canvas that it is drawn on.
        '''
        if self.canvas:
            self.canvas.delete(self.rectangle)
            self.canvas.delete(self.text)
            self.canvas.delete(self.rectid)
            self.canvas.delete(self.textid)
            if self.outtext:
                self.canvas.delete(self.outtext)
                self.outtext = None
            self.rectangle = None
            self.text      = None
            self.textid      = None
            self.rectid      = None

    def draw(self,canvas,x,y):
        '''
        Draw this flow on the given canvas, with the upper-left corner of the Flow
        at the (x,y) coordinates (given in pixels).

        If you want to change how Flows look, then this is the place to modify.
        '''
        assert isinstance(canvas,Canvas)
        self.canvas = canvas # I am drawn on this canvas
        # Create background rectangle
        self.rectangle = canvas.create_rectangle(x,y,
            x+Flow.WIDTH,
            y+Flow.HEIGHT,
            fill=Flow.BACKGROUND_COLOR,
            outline=self.OUTLINE_COLOR,
            width=3,
            tags = (self.tag,self.tag+'-rect'))

        # Create text (name)
        self.text = canvas.create_text(x+Flow.WIDTH/2,y+Flow.HEIGHT/5,
            text=self.name,
            tags=(self.tag,self.tag+'-text'))
        # Create flow tag text
        self.rectid = canvas.create_rectangle(x,y,x+15,y+15,outline='red',
            width=2,tags=(self.tag,self.tag+'-rectid'))
        self.textid = canvas.create_text(x+7,y+7,text='%d' % self.id,
            tags=(self.tag,self.tag+'-textid'))

    def getOutputLinkCoordinates(self):
        '''
        If we are piping our Output to another Flow, then we need to put that Flow
        in a particular place relative to this Flow (stair-stepping effect).  

        This method returns the coordinates where we should move the Output Flow to.
        '''
        if self.canvas:
            coords = self.canvas.coords(self.rectangle)

            x0 = coords[0] + Flow.WIDTH/4
            y0 = coords[1] + Flow.HEIGHT
            x1 = coords[2] + Flow.WIDTH + Flow.WIDTH/4 # this is the stair-step
            y1 = coords[3] + 2*Flow.HEIGHT

            return x0,y0,x1,y1

    def move(self,dx,dy):
        '''
        Move this Flow by dx and dy relative to its current location on its canvas.
        '''
        if self.canvas: 
            self.canvas.move(self.tag,dx,dy)
            self.canvas.lift(self.tag)

    def moveToCoordinates(self,coords):
        '''
        Move this flow to the following coordinates.

        coords is (x0,y0,x1,y1) the same format as is returned by canvas.coords()
        '''
        if self.canvas: # If I am drawn
            x0,y0,x1,y1 = self.canvas.coords(self.rectangle)
            n0,m0,n1,m1 = coords
            self.canvas.move(self.tag,n0-x0,m0-y0)
            self.canvas.lift(self.tag)
    
    def addLinkTag(self,newtag):
        '''
        Add a tag corresponding to the LinkedFlow that this Flow is part of
        to all canvas items in this Flow.
        '''
        if self.canvas: # I must be drawn somewhere before I can add tags
            tags = self.canvas.gettags(self.rectangle)
            newtags = list(tags)
            newtags.insert(0,newtag)
            self.canvas.itemconfig(self.rectangle, 
                tags=tuple(newtags), outline=self.LINKED_COLOR)
                
            tags = self.canvas.gettags(self.text)
            newtags = list(tags)
            newtags.insert(0,newtag)
            self.canvas.itemconfig(self.text,      tags=tuple(newtags))
            
            if self.output:
                if self.outtext != None:
                    self.canvas.delete(self.outtext)
                else:
                    self.canvas
            elif self.outtext:
                self.canvas.itemconfig(self.outtext,tags=tuple(newtags))
            
        else:
            raise Exception('Flow has not been painted to a canvas')

    def link(self,lf=None,input=None,output=None):
        '''
        Link this Flow to a LinkedFlow. 
        OR Link this Flow to another Flow's OUTPUT
        OR Link this Flow to another Flow's INPUT
        OR a combination of the above.
        '''

       
        if input: # Link the Output from input to this Flow's INPUT
            assert isinstance(input,Flow) or isinstance(input,LinkedFlow)
            self.input = input
            self.input.configure()
        if output: # Link the OUTPUT from this Flow to the input to another Flow
            assert isinstance(output,Flow) or isinstance(output,LinkedFlow)
            self.output = output
            # Move the output to the right location
            output.moveToCoordinates(self.getOutputLinkCoordinates())
        if lf: # Link this Flow to a LinkedFlow
            self.addLinkTag(lf.tag)
            self.linkedFlow = lf
        
    def unlink(self):
        '''
        Unlink this flow from any flows that it is Linked to as well as
        any LinkedFlows that it is apart of.
        '''
        if self.canvas: # Remove link tags, and reset colors
            tags = self.canvas.gettags(self.rectangle)
            newtags = list(tags)
            newtags.pop(0)
            self.canvas.itemconfig(self.rectangle,
                tags=tuple(newtags),
                outline=self.OUTLINE_COLOR)

            tags = self.canvas.gettags(self.text)
            newtags = list(tags)
            newtags.pop(0)
            self.canvas.itemconfig(self.text,tags=tuple(newtags))

        self.input  = None
        self.output = None


    def getParameters(self):
        '''
        Open a window that allows us to enter values for the Parameters
        for this Flow.

        When the window is closed, this Flow is considered to be configured.

        Corner is a tuple of (x,y) of the upper left corners of this window.
        If corner is passed, then the window will be placed at those coordinates.

        This method can be directly called.  No need for a Canvas.
        
        Presently, there is no Parameter checking to ensure that we have all
        parameters required for programs (because Madagascar documentation does not
        require these values to be set for program parameters).  If we want
        to implement parameter checking in the future we would do so in here.
        '''
        # Create a window
        window = Toplevel()
        window.grid()
        window.resizable(False,False)
        window.title("Edit parameters:")

        # Create a frame to go in this Window
        allFrame = Frame(window)
        allFrame.grid()
        
        Entries = [] # List of Entries for all parameters
        pars = self.keys 
        
        EntryFrame = Frame(allFrame)
        EntryFrame.grid(row=0,column=0)
        for ip in range(len(pars)): # For each parameter
            par = pars[ip]
            value = self.pars[par].value
            if not value: value = ''
            description = str(par) # Get tooltip description

            # Make a label with the parameters' name
            name  = Label(EntryFrame,text='%s=' % par,justify=RIGHT)

            # Create an Entry for this parameter next to its label
            # Special handling for Input and Output, incase we are piping
            if (par == 'output' and self.output) or (par == 'input' and self.input): # we are piped
                entry = Entry(EntryFrame)
                entry.insert(0,"-----PIPED-----")
                entry.configure(state=DISABLED)
            else:
                entry = Entry(EntryFrame)
                entry.insert(0,value)

            Entries.append(entry) # Save for later
            name.grid  (column=0,row=ip,sticky=E+W)
            entry.grid (column=1,row=ip,sticky=E+W)

        def __setParameters(event,preview=False,field=None):
            '''
            Get the values for the parameters from our window.
            '''

            for j in range(len(pars)):
                value = Entries[j].get()
                try:
                    self.pars[pars[j]].set(value)
                except Exception, e:
                    tkMessageBox.showwarning("Bad parameter value:",e)
                    
            if preview:
                wasSet = True

                if not self.set: 
                    self.set = True
                    wasSet = False
                
                if self.linkedFlow:
                    temp = self.linkedFlow.getCommand()
                else:
                    temp = self.getCommand()

                self.set = wasSet
                if temp: 
                    text.delete(1.0,END)
                    text.insert(END,temp)
            else:
                self.configure()
                window.destroy()
        
       
        
        tscroll = Scrollbar(allFrame,orient=VERTICAL)
        text = Text(allFrame,width=60,yscrollcommand=tscroll.set,wrap=WORD)
        tscroll['command'] = text.yview
        text.insert(END,"Click Preview to see how this Flow would be in an SConstruct.")
        text.grid(row=0,column=1,rowspan=len(pars),sticky=N+S)
        tscroll.grid(row=0,column=2,sticky=N+S,rowspan=len(pars))
        
        def __viewSelfDoc(event):
            text.delete(1.0,END)
            text.insert(END,self.selfdoc)

        b = Button(allFrame,text="Accept")
        b.bind('<Button-1>',__setParameters)
        b.grid(row=len(pars),column=0,sticky=E+W)

        bFrame = Frame(allFrame)
        bFrame.grid(row=len(pars),column=1)
        b2 = Button(bFrame,text="Preview")
        b2.bind('<Button-1>',lambda event, preview=True, field=text: __setParameters(event,preview,field) )
        b2.grid(row=0, column=0)
        
        b2 = Button(bFrame,text="View Documentation")
        b2.bind('<Button-1>', __viewSelfDoc)
        b2.grid(row=0, column=1)
        
    def unconfigure(self):
        '''
        This Flow was not properly configured, because it has an invalid output value,
        but it was marked as being configured.
        
        Unmark its configured status, and display a warning that it was not configured.
        '''
        self.set = False
        if self.canvas: 
            if self.outtext:
                self.canvas.delete(self.outtext)
                self.outtext = None
            self.canvas.itemconfigure(self.rectangle,fill=self.BACKGROUND_COLOR)    

    def drawOutputText(self):
        if self.outtext: 
            self.canvas.delete(self.outtext)
        x,y,x1,y1 = self.canvas.coords(self.rectangle)
        self.outtext = self.canvas.create_text(x+Flow.WIDTH/2,
                            y+2*Flow.HEIGHT/3,
                            text=self.pars['output'].value,
                            tags=self.canvas.gettags(self.text),
                            fill='black', font=("arial","12","bold"))
    def configure(self):
        '''
        Set this Flow as having been configured.
        '''
        if self.canvas: 
            if self.pars['output'].value:
                self.drawOutputText()
            else:
                self.unconfigure() # I was not properly configured.
                return None
            self.canvas.itemconfigure(self.tag+'-rect',fill='green')
        self.set = True



class LinkedFlow():
    ''' 
    A LinkedFlow is NOT a Flow.  It exists independently of a Flow.

    Rather it is only a container for multiple Flows.  It does not directly
    draw itself on any Canvas. 
    '''
    def __init__(self,id):
        ''' 
        Create a LinkedFlow with the given tag that uniquely identifies it
        from all other items.

        This inherits the ftype from the first Flow added to it.
        '''
        self.flows = []
        self.tag = 'link-%d' % id
        self.id  = id 
        self.ftype = None

    def __getstate__(self):
        '''
        Save information about this LinkedFlow using pickle.
        Saves the tag, id, and flow type and all tags of flows attached to 
        this linked flow.
        
        The order of the flows is preserved when pickled.
        '''
        
        return (self.ftype,self.flows)

    def __setstate__(self,state):
        '''
        Restore this LinkedFlow using pickle.
        
        Gets the tag, id, flow type and tags of flows that were associated
        with this LinkedFlow.
        
        Sets the flows attached to this to [], the list must be repopulated
        manually because we do not want independent copies of Flows floating
        in this LinkedFlow.
        '''
        self.ftype,self.flows = state
        self.tag = None
        self.id  = -1
    
    def setTag(self,id):
        self.tag = 'link-%d' % id
        self.id  = id

    def link(self,flow):
        '''
        Link a Flow to this LinkedFlow.  Flows that are added are automatically
        attached at the END of this LinkedFlow, meaning that you may not insert
        new Flows into a LinkedFlow.

        OR add a LinkedFlow to this LinkedFlow.  Same rules apply as above.
        '''

        if self.ftype and flow.ftype != self.ftype: 
            tkMessageBox.showwarning("Bad link type","You may not link together flows of different types (e.g. Flow with Result).  In this case you tried to link %s to %s!" % (self.ftype,flow.ftype))
        elif isinstance(flow,LinkedFlow): # add a LinkedFlow to this Flow
            for iflow in flow.flows:
                iflow.unlink()
                self.addFlow(iflow)
        else: # Add a Flow
            # Get the length, because that will be the index of the Flow that
            # we are attaching to this.
            index = len(self.flows)
            if index > 0: 
            # If there are already Flows, then setup the INPUT and OUTPUT
            # pipes for this flow and the flow before it.
                lastFlow = self.flows[index-1]
                lastFlow.link(output=flow)
                flow.link(lf=self,input=lastFlow)
            else:
            # Otherwise, just add this flow.
                self.ftype = flow.ftype
                flow.link(lf=self)
            self.flows.append(flow)

    def delete(self):
        '''
        Delete all Flows in this LinkedFlow.  Called when this LinkedFlow is destroyed.
        '''
        for flow in self.flows:
            flow.delete()

    def move(self,dx,dy):
        '''
        Move this LinkedFlow by dx,dy relative to its current location.
        '''
        for flow in self.flows:
            flow.move(dx,dy)

    def unlink(self,flow,tag):
        ''' 
        Unlink the flow from the list of flows attached to this LinkedFlow,
        if there are more than two flows in this LinkedFlow prior to detaching
        figure out how to divide the remaining flows into two LinkedFlows
        that wrap the removed flow index.

        Pass a tag for a new Linked Flow, if necessary.
        
        Returns: a list of unlinked Flows and/or LinkedFlows that include all
        Flows that are split as a result of this action.
        '''

        index = self.flows.index(flow)
        # Splice the list of flows into two sets, before and after
        before = self.flows[:index]
        after  = self.flows[index+1:]

        newFlows = []

        if len(before) > 1: # keep this LinkedFlow
            self.flows = []
            for iflow in before:
                self.flows.append(iflow)
            newFlows.append(self) 
        else:
            if len(before) == 1: # Get the first flow
                tflow = before[0]
                tflow.unlink()
                newFlows.append(tflow)

        if len(after) > 1: # make a new LinkedFlow
            lf = LinkedFlow(tag)
            for iflow in after:
                iflow.unlink()
                lf.addFlow(iflow)
            newFlows.append(lf)
        else:
            if len(after) == 1: # Get the last flow
                tflow = after[0]
                tflow.unlink()
                newFlows.append(tflow)
        flow.unlink()
        newFlows.append(flow)
        return newFlows

    def getCommand(self):
        '''
        Get the command string that this LinkedFlow would produce, which
        is the concatenation of the command strings for each flow
        attached to this already.
        
        We have to do some nasty logic to get the formatting correct
        for the inputs and outputs.
        '''
        
        # The first and last flows tell us what the input and output 
        # targets and sources are for this LinkedFlow.
        last = self.flows[len(self.flows)-1]
        first = self.flows[0]

        input = first.pars['input'].value
        output = last.pars['output'].value
        
        # Setup the output command
        command = self.ftype+"("
        if output and output != '':
            command += "'%s'," % output
        else:
            # If we have no output target, then we did not properly
            # configure the last flow in this LinkedFlow
            # bail on getting the command
            last.unconfigure()
            return None
            
        # Handle some nastiness here, depending on which input we have    
        if input and input != '': # String input
            command += "'%s',\n        '''\n        " % input
        elif self.ftype == Flow.PLOT or self.ftype == Flow.RESULT: 
        # If our Input is None and we are a Plot or a Result
        # then we don't want to have a None input.
            command += "\n        '''\n        "
        else:
        # Otherwise, we should have a None input.
            command += "None,\n        '''\n        "
        
        # Check to make sure that all flows
        # are configured.  If not, bail.
        set = True
        # For all other flows
        for flow in self.flows:
            if not flow.set:
                return None
            else:
                command += flow.getCommand()
              
        return command
        
        
        
        
