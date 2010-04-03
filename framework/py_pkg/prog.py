import sys, os

import rsf

import sfplot
import sflibplot
import sfpens
import sfsumain
import sfsuplot
import sfgeneric
import sfseismic
import sfmain
import sfwemva
import sfmccowan
import sfslim
import sfpsava
import sfjeff
import sfroman
import sfkourkina
import sfgodwinj
import sfbash
import sfjyan
import sfsolv
import sfpetsc
import sftaup
import sfsongxl
import sfstatoilutil
import sfgchliu
import sfslsmig
import sfgee
import sfjennings
import sftrip
import sfbrowaeys
import sffomels
import sfeffsilva
import sfyliu
import sfllisiw
import sfivlad

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
