

class typecheck(object):

    def __init__(self, *args, **kws):
        self.args = args
        self.kws  = kws

    def typecheck(self, value, expected, name=''):
        typ = type(value)
        if typ is not expected: raise TypeError, '%s expected %s but got %s' % (name or 'param', expected, typ)

    def __call__(self, fn):
        def wrapped(*args, **kws):
            for v, t in zip(args, self.args):
                self.typecheck(v, t)
            for n, e in self.kws.iteritems():
                if n in kws:
                    self.typecheck(kws[n], e, n)

            return fn(*args,**kws)
        return wrapped
