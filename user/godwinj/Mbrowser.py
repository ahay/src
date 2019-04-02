#!/usr/bin/env python
'''
Open tkMadagascar program browser.
'''
from __future__ import print_function
import sys

try:
    from Tkinter import *
except Exception as e:
    print('Could not find Tkinter installation: ', e)
    sys.exit(1)

try:
    import rsf.gui.Browser as browser
    browser.launch()
except Exception as e:
    print('Failed to start Browser: ', e)


