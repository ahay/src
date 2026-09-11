# Example flow demonstrating the rsf.proj DSL.
#
# Run with:
#   cp example-flow.SConstruct /tmp/writing-rsf-flows-demo/SConstruct
#   cd /tmp/writing-rsf-flows-demo && scons

from rsf.proj import *

# A Flow(target, source, command) produces <target>.rsf by running <command>.
# "None" as the source means the flow has no RSF input.
Flow('spike', None, 'spike n1=1000 k1=300')

# Multi-stage flow via pipe. The "sf" prefix is implicit inside a Flow string.
Flow('filter', 'spike', 'bandpass fhi=2 phase=y')

# A parameter can be interpolated at SConstruct time.
fhi = 4
Flow('filter2', 'spike', 'bandpass fhi=%g phase=y' % fhi)

# A Plot wraps a vplot program. It produces a .vpl file but is not a "final" result.
Plot('filter', 'wiggle clip=0.02 title="Filtered spike"')
Plot('filter2', 'wiggle clip=0.02 title="Filter (fhi=%g)"' % fhi)

# Result == Plot, but marks the output as a deliverable (goes into Fig/ and
# is shown by `scons view`).
Result('filter', 'wiggle clip=0.02 title="Filtered spike"')

End()
