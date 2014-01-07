import os

_methods = {
    'gazdag': 'Phase-Shift Migration',
    'kirch': 'Kirchhoff Time Migration',
    'lowrank': 'Two-Step Lowrank',
    'oway': 'One-Way Wave Equation',
    'stolt': 'Stolt Migration',
    'vc': 'Velocity Continuation',
    'fowler': 'Fowler',
    'fakirmig': 'First-Arrival Kirchhoff',
    'ffd': 'FFD',
    'lrostep': 'One-Step Lowrank',
    }

def method():
    thisdir = os.path.basename(os.getcwd())
    if thisdir in _methods.keys():
        return _methods[thisdir]
    else:
        return ''
    
