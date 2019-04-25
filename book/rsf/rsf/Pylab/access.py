from pylab import *
switch_backend('agg')

labels = '''
google.com
direct
sep.stanford.edu
yahoo.com
seg.org
opendtect.org
other
'''.split()

#cwp.mines.edu

values = [
15827,
6610,
679,
639,
574,
301,
187+1808]

figure(1, figsize=(8,8))
pie(values, labels=labels)
title('Source Referral')

savefig('junk_py.eps')
