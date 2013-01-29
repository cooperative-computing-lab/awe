
import work_queue as WQ

import os

WQ.set_debug_flag('all')

Q = WQ.WorkQueue(name = 'test-badi',
                 port = WQ.WORK_QUEUE_RANDOM_PORT,
                 exclusive = True,
                 catalog = True,
                 shutdown = False)

script = '/tmp/test.sh'
with open(script, 'w') as fd:
    fd.write("""\
#!/usr/bin/env bash
echo "running"
echo "sleeping"
sleep 1
echo "done"
""")

os.system('chmod a+x %s' % script)

print 'Submitting'

for i in xrange(100):
    print i
    T = WQ.Task('./test.sh %s' % i)
    T.specify_file(script, 'test.sh')
    T.specify_tag(str(i))
    Q.submit(T)

print 'Waiting'
while not Q.empty():
    R = Q.wait(5)
    if R:
        print R.tag
        print R.output
        print R.return_status
        print R.result
        print
