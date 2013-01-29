
from awe.workqueue import TagSet

import work_queue as WQ
import random
import datetime
import time

from collections import defaultdict

random.seed(42)

DEBUG = False
WAIT = 1
GPU_PROB = 0.5
NTASKS = 500
MAXREPS = float('inf')

LOGFD = open('cancel.log', 'w')
ETIMESFD = open('exectimes.log', 'w')
LOGHEADER = False

TAGS = TagSet(maxreps=MAXREPS)


if DEBUG:
    WQ.set_debug_flag('all')


def now():
    t = datetime.datetime.now()
    return time.mktime(t.timetuple()) + t.microsecond/10.**6


def task(i):
    p = random.random()
    if p < GPU_PROB:
        bounds = (1,2)
    else:
        bounds = (20,90)
    t = WQ.Task('echo %s;sleep %sm' % (i, random.randint(*bounds)))
    t.specify_tag(str(i))
    return t

def mktasks(count):
    for i in xrange(count):
        yield task(i)

def total_tasks(Q): return Q.stats.tasks_running + Q.stats.tasks_waiting
def total_workers(Q): return Q.stats.workers_busy + Q.stats.workers_cancelling + Q.stats.workers_ready

def can_duplicate(Q, tags):
    return len(tags) > 0 and tags.can_duplicate() and total_tasks(Q) < total_workers(Q)

def submit(t, Q):
    Q.submit(t)
    print 'Adding tag', t.tag, 'to TAGS'
    TAGS.add(t.tag)

def status(Q): return Q.stats.tasks_running, Q.stats.tasks_waiting, total_tasks(Q), total_workers(Q), len(TAGS)

def log(Q):
    global LOGHEADER
    t = now()
    if not LOGHEADER:
        LOGFD.write('# time tasks_running tasks_waiting total_tasks workers_busy workers_init workers_ready total_workers\n')
        LOGHEADER = True

    status = '%s %s %s %s %s %s %s %s' % \
        (t,
         Q.stats.tasks_running, Q.stats.tasks_waiting, total_tasks(Q),
         Q.stats.workers_busy, Q.stats.workers_init, Q.stats.workers_ready, total_workers(Q))
    # print status
    LOGFD.write(status + '\n')
    LOGFD.flush()

def log_result(t):
    ETIMESFD.write('%d\n' % t.cmd_execution_time)
    ETIMESFD.flush()


def main():
    Q = WQ.WorkQueue(name='test-cancel', catalog=True, exclusive=True)
    print 'Port:', Q.port

    log(Q)
    for t in mktasks(NTASKS):
        log(Q)
        submit(t, Q)
        log(Q)

    while total_tasks(Q) > 0:
        log(Q)
        t = Q.wait(WAIT)
        log(Q)
        if t and t.result == 0:
            TAGS.discard(t.tag)
            print 'Got', t.tag, t.cmd_execution_time / 10.**6
            log_result(t)
            while Q.cancel_by_tasktag(t.tag):
                print now(), 'canceled', t.tag, status(Q), TAGS
                log(Q)

        while can_duplicate(Q, TAGS):
            tag = TAGS.select()
            if tag is None: break
            print now(), 'duplicating', tag
            log(Q)
            submit(task(tag), Q)
            log(Q)
            TAGS.clean()

        print status(Q), TAGS

if __name__ == '__main__': main()
