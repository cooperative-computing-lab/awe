
import traceback

class TypeException (Exception): pass

class _typecheck(object):

    def __init__(self, method=True):
        self.method = method

    def check(self, *args, **kws):
        if self.method:
            self.args = [None] + list(args)
        else:
            self.args = args
        self.kws  = kws

        return self

    def typecheck(self, value, expected, name='', arg=-1):
        typ = type(value)
        if typ is not expected:
             raise TypeException, '%s expected: %s, but got: %s' % (name or 'param %s' % arg, expected, typ)

    def __call__(self, fn):
        def wrapped(*args, **kws):
 
            try:
                i = -1
                for v, t in zip(args, self.args):
                    i += 1
                    if self.method: continue
                    self.typecheck(v, t, arg=i)
                del i

                for n, e in self.kws.iteritems():
                    if n in kws:
                        self.typecheck(kws[n], e, name=n)
            except TypeException, ex:
                stack = traceback.extract_stack()
                stack = ''.join(traceback.format_list(stack[:-1]))
                stack = '\n\t'.join(('\t' + stack).split('\n'))
                raise TypeException, '%s:\n\n%s\n\t%s' % (fn, stack, ex)

            return fn(*args,**kws)
        return wrapped

def typecheck(*args, **kws):
    tc = _typecheck(method=True)
    tc.check(*args, **kws)
    return tc


def typecheckfn(*args, **kws):
    tc = _typecheck(method=False)
    tc.check(*args, **kws)
    return tc
