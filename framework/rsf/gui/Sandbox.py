from Tkinter import *
from Flow import Flow
from Flow import LinkedFlow

import tkMessageBox

class Sandbox(Canvas):
    '''
    The Sandbox class contains the Canvas where all of the graphical action
    happens.
    
    Sandboxes are responsible for keeping track of Flows attached to them,
    and allow users to decide how to interact with Flows.
    
    A Sandbox must be attached to a Toplevel or root widget.
    '''  

    HEIGHT = 600
    WIDTH  = 800
    SCROLL_HEIGHT = 2400
    SCROLL_WIDTH  = 3200
    def __init__(self, master):
        Canvas.__init__(self, master,
            width=Sandbox.WIDTH,
            height=Sandbox.HEIGHT,
            scrollregion=(0,0,
                Sandbox.SCROLL_WIDTH,
                Sandbox.SCROLL_HEIGHT)
        )   

        # self.flows lists all flows, plots and results on the canvas
        # self.items lists all independent elements on the canvas
        #     (i.e. flows,plots,results and LinkedFlows, LinkedPlots and
        #           LinkedResults)
        self.flows = {}
        self.items = {}

        self.dragging = None
        self.drag_coords = (0,0)

        self.linkid = 0 # number of linked flows,plots or results
        self.flowid = 0 # number of flows, plots or results
        self.__setMouseBindings()

    def __setMouseBindings(self):
        '''
        Set canvas mouse bindings.
        '''
        self.bind('<B1-Motion>',self.__dragItem)
        self.bind('<ButtonRelease-1>',self.__endDrag)
        self.bind('<ButtonPress-1>',self.__startDrag)
        self.bind('<Button-3>',self.__openPopupMenu)

    def getItems(self):
        '''
        Return a list of only separate items attached to this canvas.
        
        i.e. LinkedFlows and Flows, Plots and Results that are not
        part of a LinkedFlow.
        '''
        items = []
        for tag,item in self.items.items():
            items.append(item)
        return items

    def getFlowID(self):
        '''
        Get a unique ID to identify a flow on this canvas.
        '''
        oid = self.flowid
        self.flowid += 1
        return oid

    def getLinkID(self):
        '''
        Get a unique ID to identify a linked flow on this canvas.
        '''
        oid = self.linkid
        self.linkid += 1
        return oid

    def __findItemAtMouse(self,x,y,index=-1):
        '''
        Get the LinkedFlow or Flow at the location (x,y) in this
        canvas.
        
        Return a list of Flows at this location, in top-bottom order.
        
        If index, return the item corresponding to that index from the list.
        Example:  
        -> If there are two Flows then index=0 returns the top-most Flow.
        ->        OR                   index=1 returns the bottom   Flow.
        
        For LinkedFlows: index=1 returns the Flow under the cursor instead of 
        the LinkedFlow.
        '''       
        ids = self.find_overlapping(x,y,x+1,y+1)
        unique = None
        if len(ids) > 0:
            unique = []
            idlist = list(ids)
            idlist.reverse() # Get top-most first
            for item in idlist:
                tags = self.gettags(item)
                for tag in tags:
                    try:
                        flow = self.items[tag]
                        if flow not in unique:
                            unique.append(flow)
                    except:
                        pass
                    try:
                        flow = self.flows[tag]
                        if flow not in unique:
                            unique.append(flow)
                    except:
                        pass
        if index >= 0 and unique:
            return unique[index]
        else:
            return unique

    def __openPopupMenu(self,event):
        '''
        Open a popup menu at the Flow or LinkedFlow under
        the cursor.
        '''
        
        x = event.x
        y = event.y

        flow = self.__findItemAtMouse(x,y,0)

        def __openDialog():
            ''' 
            Open the Edit Parameters window for the Flow at
            this location.
            '''
            # Get the first flow that is at this location
            flows = self.__findItemAtMouse(x,y) 
            if isinstance(flows[0],LinkedFlow):
                try:
                    flow = flows[1]
                except:
                    pass
            else:
                try:
                    flow = flows[0]
                except:
                    pass

            if flow: flow.getParameters()

        def __linkFlows(event, iflow):
            ''' 
            Link iflow to the flow or LinkedFlow that the user selects.
            
            If nothing is selected, then do nothing, but show a messagebox.
            '''
            x = event.x
            y = event.y

            flow = self.__findItemAtMouse(x,y,0)

            try:
                if isinstance(flow,LinkedFlow):
                    flow.link(iflow)
                    self.items.pop(iflow.tag)
                else:
                    if flow != None and flow != iflow:
                        self.addLinkedFlow(flow,iflow)
                    else:
                        raise Exception("Invalid flow selected to link to.")
            except Exception, e:
                if flow:
                    tkMessageBox.showwarning("Bad attempt at linking",e)
            finally:
                self.unbind("<Button-1>")
                self.__setMouseBindings()
                self.configure(cursor="arrow")

        def __linkFlow():
            ''' 
            Link the flow that the popup is over to another flow that
            the user selects.
            
            Change the mouse cursor to indicate that the user is selecting
            an item.
            '''
            try:
                self.bind('<Button-1>', 
                    lambda event, arg=flow: __linkFlows(event,arg) )
                self.configure(cursor='crosshair')
            except:
                print 'bad attempt at linking'

        def __unlinkFlow():
            '''
            Find the LinkedFlow and the Flow that are under the cursor
            and remove that Flow from its LinkedFlow.
            '''
            
            flows = self.__findItemAtMouse(x,y)
            try:
                linked = flows[0]; flow = flows[1]
                self.items.pop(linked.tag)
                
                split = linked.unlink(flow,self.getLinkID())

                for flow in split:
                    self.items[flow.tag] = flow
            except Exception, e:
                tkMessageBox.showwarning("Failed to unlink:",e)

        def __deleteFlow():
            '''
            Delete the Flow or LinkedFlow under the cursor.
            '''
            try:
                if isinstance(flow,LinkedFlow):
                    if tkMessageBox.askyesno("Delete Linked Flow?",
                        "Delete all flows attached to this Linked Flow?  (To delete a single flow, unlink it first)"):
                        flow.delete()
                        self.items.pop(flow.tag)
                else:
                    if tkMessageBox.askyesno("Delete Flow?", 
                        "Are you sure you want to delete this %s?" % flow.ftype):
                        
                        flow.delete()
                        self.items.pop(flow.tag)
                        self.flows.pop(flow.tag)
            except Exception, e:
                tkMessageBox.showwarning("Failed to delete flow:", e)

        if flow: # we found an item
            # create a menu
            popup = Menu(self, tearoff=0)
            popup.add_command(label="Edit Parameters",command=__openDialog) 
            popup.add_command(label="Link",command=__linkFlow)
            popup.add_command(label="Unlink",command=__unlinkFlow)
            popup.add_separator()
            popup.add_command(label="Delete",command=__deleteFlow)
            popup.add_separator()
            popup.add_command(label="Close this menu...")

            # Snagged off the internet
            try:
                popup.tk_popup(event.x_root, event.y_root, 0)
            finally:
                # make sure to release the grab (Tk 8.0a1 only)
                popup.grab_release()

    def __addFlowAtMouse(self,event,flow):
        '''
        Let the user specify where to place the newly created
        Flow by clicking somewhere on this canvas. 
        
        Afterwards, fix the MouseBindings and cursor.
        '''
        x = event.x
        y = event.y
        flow.draw(self,x,y)
        self.unbind('<Button-1>')
        self.configure(cursor="arrow")
        self.__setMouseBindings()
        self.items[flow.tag] = flow
        self.flows[flow.tag] = flow

    def addFlow(self,program,ftype):
        '''
        Add a Flow on this canvas, but let the user choose
        where to place it.
        '''
        flow = Flow(self.getFlowID(),program,ftype)
        self.configure(cursor="crosshair")
        self.bind('<Button-1>',
            lambda event, arg=flow: self.__addFlowAtMouse(event,arg) )

    def loadFlow(self,flow):
        '''
        Restore a Flow from a pickled state.
        '''
        # Set the unique flow count to be higher than the
        # maximum value from the pickled flows.
        flow.setTag(self.getFlowID())

        # Draw this flow where it was.
        if flow.coords:
            flow.draw(self,flow.coords[0],flow.coords[1])

        # If it was configured, then configure it.
	    if flow.set:
	        flow.configure()
        
        # Keep track of this flow.
        self.flows[flow.tag] = flow
        self.items[flow.tag] = flow
        
    def addLinkedFlow(self,flow1,flow2):
        '''
        Add a LinkedFlow by attaching two Flows to a new LinkedFlow.
        '''
        lf = LinkedFlow(self.getLinkID())
        lf.link(flow1)
        lf.link(flow2)

        self.items[lf.tag] = lf
        self.items.pop(flow1.tag)
        self.items.pop(flow2.tag)

    def loadLinkedFlow(self,lf):
        '''
        Restore a LinkedFlow from a pickled state.  
        You MUST restore all Flows that this LinkedFlow has attached
        to this Canvas before you may call this method.
        '''
        # Which flows are attached?
        lf.setTag(self.getLinkID())
        
        tflows = []
        for flow in lf.flows:
            self.loadFlow(flow)
            self.items.pop(flow.tag)
            tflows.append(flow)
        
        lf.flows = []
        for flow in tflows:
            lf.link(flow)
        
        self.items[lf.tag] = lf
        
    def reset(self):
        '''
        Reset this canvas.  Delete all objects attached to it.
        '''
        for tag,item in self.items.items():
            item.delete()
            if tag in self.flows.keys():
                self.flows.pop(tag)
        
        for tag,item in self.flows.items():
            item.delete()
        
        self.items = {}
        self.flows = {}
        self.linkid = 0
        self.flowid = 0
        

    def __startDrag(self,event):
        '''
        Start dragging an object if there is one under the cursor.
        '''
        x = event.x; y = event.y
        flow = self.__findItemAtMouse(x,y,0)
        if flow: self.dragging = flow
        self.drag_coords = (x,y)

    def __dragItem(self,event):
        '''
        Move the object that the mouse selected.
        '''
        x = event.x
        y = event.y
        if self.dragging:
            dx = x-self.drag_coords[0]
            dy = y-self.drag_coords[1]
            self.dragging.move(dx,dy)
        self.drag_coords = (x,y)

    def __endDrag(self,event):
        '''
        Quit dragging the object.
        '''
        self.dragging = None 
