# GP-Localize

Running algo:

run.sh

Revising the following configuration:

'dat' is for dataset, refer to #traj.

traj0: test on IEQ data, traj name: working

traj1: test on wireless (single, multiple fields), traj name: work

traj2: test on real-world pioneer data, traj name: real

'''
envwork: IEQ data
wirework: WSS data
real: real-world pioneer data
'''

field info:

fieldnum_numbering

field number starts from 0

sh genCase.sh [dataset id] [traj name] [field list] [algo list] [no. of steps] [output file]
