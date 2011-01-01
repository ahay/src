#!/usr/bin/env python
''' 
Open tkMadagascar GUI
'''
import sys
try:
    from Tkinter import *
except Exception, e:
    print 'Could not find Tkinter installation: ', e
    sys.exit(1)
    
try:
    import rsf.gui.Gui as gui
    gui.launch()
except Exception, e:
    print 'Failed to start tkMadagascar: ', e
