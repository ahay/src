SConscript(dirs=['seis/rsf','seis/main','seis/proc','seis/imag'],
           name='SConstruct')
SConscript(dirs=['vplot/lib','vplot/main'],
           name='SConstruct')

import os

env = Environment()

doc = """echo 'import sys, os' > $TARGET
echo "sys.path.append(os.environ.get('RSFROOT'))" >> $TARGET
echo "import rsfdoc" >> $TARGET
echo " " >> $TARGET
echo "rsfprog={}" >> $TARGET
cat $SOURCES >> $TARGET"""

mains = ['seis/main','seis/proc','seis/imag','vplot/main']

env.Command('rsfprogs.py',map(lambda x: os.path.join(x,'RSFdoc'),mains),doc)

Clean('rsfprogs.py','rsfprogs.pyc')
