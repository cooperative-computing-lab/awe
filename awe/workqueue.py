
import awe

import work_queue as WQ

import os, tempfile


### A process can only support a single WorkQueue instance
_AWE_WORK_QUEUE = None


### workaround for now.
##+ These are the names of the input/output filess to be materialized on the worker

WORKER_PDB_NAME     = 'structure.pdb'
WORKER_WEIGHTS_NAME = 'weights.dat'
WORKER_COLOR_NAME   = 'color.dat'
WORKER_CELL_NAME    = 'cell.dat'
WORKER_RESULTS_NAME = 'results-%s.tar'

RESULT_POSITIONS    = 'structure2.pdb'
RESULT_WEIGHTS      = 'weights2.dat'
RESULT_COLOR        = 'color2.dat'
RESULT_CELL         = 'cell2.dat'



class WorkQueueException       (Exception): pass
class WorkQueueWorkerException (Exception): pass

class Config(object):
    """
    Class for configuring a WorkQueue instance
    """

    def __init__(self):

        self.name      = 'awe'
        self.port      = WQ.WORK_QUEUE_RANDOM_PORT
        self.schedule  = WQ.WORK_QUEUE_SCHEDULE_TIME
        self.exclusive = True
        self.catalog   = True
        self.debug     = 'all'
        self.shutdown  = False
        self.fastabort = -1

        self.waittime  = 10 # in seconds


        self._executable = None
        self._cache = set()

    executable = property(lambda self: self._executable)
    getcache   = property(lambda self: self._cache)

    def execute(self, path):
        self._executable = path
        self.cache(path)

    def cache(self, *files):
        for path in files:
            self._cache.add(path)

    def _mk_wq(self):
        global _AWE_WORK_QUEUE
        if _AWE_WORK_QUEUE is not None:
            ### warn
            awe.log('WARNING: using previously created WorkQueue instance')
            return _AWE_WORK_QUEUE
        else:
            WQ.set_debug_flag(self.debug)
            wq = WQ.WorkQueue(name      = self.name,
                              port      = self.port,
                              shutdown  = self.shutdown,
                              catalog   = self.catalog,
                              exclusive = self.exclusive)
            wq.specify_algorithm(self.schedule)

            typ = type(self.fastabort)
            if typ is float or typ is int:
                wq.activate_fast_abort(self.fastabort)

            _AWE_WORK_QUEUE = wq

            return wq


class WorkQueue(object):
    def __init__(self, cfg):

        self.cfg    = cfg
        self.wq     = self.cfg._mk_wq()

        self.stats  = awe.stats.WQStats()

        self.tmpdir = tempfile.mkdtemp(prefix='awe-tmp.')


    empty = property(lambda self: self.wq.empty())

    def __del__(self):
        import os
        os.rmdir(self.tmpdir)


    @awe.trace()
    def update_wq_stats(self):
        self.stats.wq(self.wq)

    @awe.trace()
    def update_task_stats(self, task):
        self.stats.task(task)

    @awe.trace()
    def new_task(self, params):
        cmd = self.cfg.executable
        task = WQ.Task('./' + cmd)

        ### executable
        task.specify_file(self.cfg.executable)

        ### cached files
        for path in self.cfg.getcache:
            task.specify_file(path)

        ### convert the walker parameters for WQWorker
        task.specify_buffer(params['weight'] , WORKER_WEIGHTS_NAME , cache=False)
        task.specify_buffer(params['color']  , WORKER_COLOR_NAME   , cache=False)
        task.specify_buffer(params['cell']   , WORKER_CELL_NAME    , cache=False)
        task.specify_buffer(params['pdb']    , WORKER_PDB_NAME     , cache=False)
        task.specify_tag   (params['id'])

        ### result file:
        result = os.path.join(self.tmpdir, WORKER_RESULTS_NAME % task.tag)
        task.specify_output_file(result, cache=False)


        return task

    @awe.trace()
    def submit(self, task):
        return self.wq.submit(task)

    @awe.trace()
    def wait(self, *args, **kws):
        return self.wq.wait(*args, **kws)


    @awe.trace()
    def _load_result_file(task):

        path = task.output_files[0]
        with tarfile.open(path) as tar:

            pdbstring    = tar.getmember(RESULT_POSITIONS ).tobuf()
            weightstring = tar.getmember(RESULT_WEIGHTS   ).tobuf()
            colorstring  = tar.getmember(RESULT_COLOR     ).tobuf()
            cellstring   = tar.getmember(RESULT_CELL      ).tobuf()

            ss           = awe.io.StringStream(pdbstring)

            walker       = awe.aweclasses.Walker(
                coords   = mdtools.prody.parsePDBStream(ss).getCoords(),
                weights  = float(weightstring),
                color    = awe.aweclasses.Color(colorstring),
                cell     = int(cellstring),
                wid      = int(task.tag)
                )

        os.unlink(path)
        return walker


    @awe.trace()
    def recv(self):

        while True:

            awe.log('DEBUG waiting for task')
            task = self.wait(self.cfg.waittime)
            self.update_wq_stats()

            if task:

                if not task.return_status == 0:
                    raise WorkQueueWorkerException, \
                        task.output + '\n\nTask %s failed with %d' % (task.tag, task.return_status)

                self.update_task_stats(task)
                return self._load_result_file(task)
