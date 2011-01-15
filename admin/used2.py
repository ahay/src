#!/usr/bin/env python
'''
Finds if programs in the current directory
are used in any official examples.
'''
import glob, os
import rsf.doc
import rsf.prog

book = os.path.join(os.environ.get('RSFSRC'),'book')

names = []
for main in glob.glob('M*.c'):
    name = main[1:-2]
    uses = rsf.doc.progs['sf'+name].uses
    if uses:
        print name + ': ' + str(uses)
    else:
        names.append(main)
print 'svn delete ' + ' '.join(names)

for sub in glob.glob('[a-z]*.c'):
    head = sub[:-1]+'h'
    print sub + ':\n'
    grep = os.popen('grep %s *.c' % head, 'r').read()
    if not grep:
        os.system('svn delete %s' % sub)
    else:
        print grep

