
import awe

import work_queue as WQ


### A process can only support a single WorkQueue instance
_AWE_WORK_QUEUE = None

class WorkQueueException       (Exception): pass
class WorkQueueWorkerException (Exception): pass

class Config(object):
    """
    Class for configuring a WorkQueue instance
    """

    def __init__(self):

        self.name      = 'awe'
        self.port      = WQ.WORK_QUEUE_RANDOM_PORT
        self.schedule  = WQ.WORK_QUEUE_SCHEDULE_FCFS
        self.exclusive = True
        self.catalogue = True
        self.debug     = 'all'
        self.shutdown  = False
        self.fastabort = 0

        self.waittime  = 10 # in seconds


        self._executable = None
        self._cache = list()

    executable = property(lambda self: self._executable)
    cache      = property(lambda self: self._cache)

    def execute(self, path):
        self._executable = path

    def cache(self, *files):
        for path in files:
            self._cache.append(path)

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
                              catalogue = self.catalogue,
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

    def new_task(self, params):
        cmd = self.cfg.executable
        task = WQ.Task('./' + cmd)

        ### executable
        awe.log('Executable: %s' % self.cfg.executable)
        task.specify_file(self.cfg.executable)

        ### cached files
        for path in self.cfg.cache:
            awe.log('Caching %s' % path)
            task.specify_file(path)

        ### result file:
        result = os.path.join(self.tmpdir, awe.io.WORKER_RESULTS_NAME % task.tag)
        awe.log('Execting result: %s' % result)
        task.specify_output_file(result, cache=False)

        ### convert the walker parameters for WQWorker
        task.specify_buffer(params['weight'] , awe.io.WORKER_WEIGHTS_NAME , cache=False)
        task.specify_buffer(params['color']  , awe.io.WORKER_COLOR_NAME   , cache=False)
        task.specify_buffer(params['cell']   , awe.io.WORKER_CELL_NAME    , cache=False)
        task.specify_buffer(params['pdb']    , awe.io.WORKER_PDB_NAME     , cache=False)
        task.specify_tag   (params['id'])


        return task

    def submit(self, task):
        return self.wq.submit(task)

    def wait(self, *args, **kws):
        return self.wq.wait(*args, **kws)


    def _load_result_file(task):

        path = task.output_files[0]
        with tarfile.open(path) as tar:

            pdbstring    = tar.getmember(awe.io.RESULT_POSITIONS ).tobuf()
            weightstring = tar.getmember(awe.io.RESULT_WEIGHTS   ).tobuf()
            colorstring  = tar.getmember(awe.io.RESULT_COLOR     ).tobuf()
            cellstring   = tar.getmember(awe.io.RESULT_CELL      ).tobuf()

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


    def recv(self):

        while True:
            task = self.wait(self.cfg.waittime)
            self.update_wq_stats()

            if task:

                if not task.return_status == 0:
                    raise WorkQueueWorkerException, \
                        task.output + '\n\nTask %s failed with %d' % (task.tag, task.return_status)

                return self._load_result_file(task)
