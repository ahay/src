import sys, os

import rsf

import sfplot
import sflibplot
import sfpens
import sfsumain
import sfsuplot
import sfgeneric
import sfmain
import sfseismic
import sfbash
import sfbrowaeys
import sfcuda
import sfeffsilva
import sffomels
import sfgchliu
import sfgee
import sfgodwinj
import sfhess
import sfivlad
import sfjeff
import sfjennings
import sfjyan
import sfkourkina
import sfllisiw
import sfmccowan
import sfpetsc
import sfpsava
import sfroman
import sfslim
import sfsongxl
import sftrip
import sfyliu

import vpplot


def selfdoc():
    'Display man page'
    prognm = os.path.basename(sys.argv[0])
    if prognm[0] == 'M' and prognm[-3:] == '.py':
        # User testing code in his local directory
        prognm = 'sf' + prognm[1:-3]
        msg = 'Self-doc may be out of synch, do "scons install" in RSFSRC'
        sys.stderr.write(msg+'\n')

    prog = rsf.doc.progs.get(prognm)
    if prog != None:
        prog.document()
    else:
        sys.stderr.write('No installed man page for ' + prognm+'\n')
