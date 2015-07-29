#!/usr/bin/env python
'''
Finds if programs in the current directory
are used in any official examples.
'''
import glob, os

book = os.path.join(os.environ.get('RSFSRC'),'book')

for main in glob.glob('M*.c'):
    name = main[1:-2]
    print name + ':\n'
    grep = os.popen('grep %s %s/*/*/*/SConstruct %s/packages/*.py' %
                    (name,book,book), 'r').read()
    if not grep:
        #        os.system('svn delete --force %s' % main)
        print "NO"
    else:
        print grep

for sub in glob.glob('[a-z]*.c'):
    head = sub[:-1]+'h'
    print sub + ':\n'
    grep = os.popen('grep %s *.c' % head, 'r').read()
    if not grep:
        os.system('svn delete %s' % sub)
    else:
        print grep

