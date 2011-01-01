#!/usr/bin/env python
''' 
Open tkMadagascar program browser.
'''
import sys

try:
    from Tkinter import *
except Exception, e:
    print 'Could not find Tkinter installation: ', e
    sys.exit(1)
 
try:
    import rsf.gui.Browser as browser 
    browser.launch()
except Exception, e:
    print 'Failed to start Browser: ', e
        

