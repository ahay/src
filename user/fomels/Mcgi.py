#!/usr/bin/env python
'A generic CGI script'

##   Copyright (C) 2012 University of Texas at Austin
##  
##   This program is free software; you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation; either version 2 of the License, or
##   (at your option) any later version.
##  
##   This program is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##  
##   You should have received a copy of the GNU General Public License
##   along with this program; if not, write to the Free Software
##   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
import os, sys

import cgi
import cgitb; cgitb.enable()  # for troubleshooting

import shutil, time

class Lock(object):
    'lock to avoid concurrent access'
    def __init__(self):
        self.locked = False
        self.sysout = sys.stdout.fileno()
        self.stdout = os.dup(self.sysout)
        self.devnull = open(os.devnull,'w',0)
 
    def lock(self):
        assert not self.locked

        while 1:
            try:
                os.mkdir('.lock')
                self.locked = True
                break
            except:
                time.sleep(1)

        # silence stdout
        os.dup2(self.devnull.fileno(),self.sysout)

    def unlock(self):
        assert self.locked
        self.locked = False
        os.rmdir('.lock')

        # unsilence stdout
        os.dup2(self.stdout,self.sysout)

    # auto-unlock when lock object is deleted
    def __del__(self):
        if self.locked:
            self.unlock()
        self.devnull.close()

def picture(directory,figure,par):
    png = figure+".png"
    fig = os.path.join("Fig",figure+".vpl")

    lock = Lock()

    if os.chdir(directory):
        print "Content-type: text/html\n"
        print "<html><body>Wrong directory \"%s\".</body></html>" % directory
    else:
        lock.lock()
        
        fail = \
            os.system("source env.sh && scons %s %s" % (par,fig)) or \
            os.system("source env.sh && vpconvert pen=gd fat=3 serifs=n bgcolor=b %s %s" % (fig,png))

        lock.unlock()

        if fail:
            print "Content-type: text/html\n"
            print "<html><body>Madagascar failure.</body></html>"
        else:    
            print "Content-type: image/png\n"
            shutil.copyfileobj(open(png,"rb"), sys.stdout)

if __name__ == "__main__":
    form = cgi.FieldStorage()

    d = form.getvalue("dir")
    f = form.getvalue("fig")

    if not d or not f:
        print "Content-type: text/html\n"
        print "<html><body>Need dir= and fig=</body></html>"
        sys.exit(1)

    par = ''
    for key in form.keys():
        if key != "dir" and key != "fig":
            par += ' %s=%s' % (key,form.getvalue(key))

    picture(d,f,par)
