######################################################################
# SELF-DOCUMENTATION
######################################################################
import rsfdoc

def selfdoc(target=None,source=None,env=None):
    src = str(source[0])
    doc = open(str(target[0]),"w")
    rsfdoc.getprog(src,doc)
    doc.close()

Doc = Builder (action = Action(selfdoc),src_suffix='.c',suffix='.doc')
#######################################################################

env = Environment()
env.Append(BUILDERS = {'Doc' : Doc})

Export('env')
SConscript(dirs=['seis/rsf','seis/main','seis/proc','seis/imag'],
           name='SConstruct')
SConscript(dirs=['vplot/lib','vplot/main'],
           name='SConstruct')

import os

doc = """echo 'import sys, os' > $TARGET
echo "sys.path.append(os.environ.get('RSFROOT'))" >> $TARGET
echo "import rsfdoc" >> $TARGET
echo " " >> $TARGET
echo "rsfprog={}" >> $TARGET
cat $SOURCES >> $TARGET"""

mains = ['seis/main','seis/proc','seis/imag','vplot/main']

env.Command('rsfprogs.py',map(lambda x: os.path.join(x,'RSFdoc'),mains),doc)
Clean('rsfprogs.py','rsfprogs.pyc')
Default('rsfprogs.py')

