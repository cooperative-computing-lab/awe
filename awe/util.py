"""
This file is part of AWE
Copyright (C) 2012- University of Notre Dame
This software is distributed under the GNU General Public License.
See the file COPYING for details.
"""


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

        wrapped.func_name = fn.func_name
        wrapped.func_doc = fn.func_doc
        return wrapped

class returns(object):
    def __init__(self, expected):
        self._expected = expected

    @property
    def expected(self): return self._expected

    def typecheck(self, value):
        typ = type(value)
        return typ is self.expected

    def __call__(self, fn):
        def wrapped(*args, **kws):
            result = fn(*args, **kws)
            if self.typecheck(result):
                return result
            else:
                stack = traceback.extract_stack()
                stack = ''.join(traceback.format_list(stack[:-1]))
                stack = '\n\t'.join(('\t' + stack).split('\n'))
                raise TypeError, 'Result of %s(*%s, **%s) should be %s but is %s' % \
                    (fn, args, kws, self.expected, type(result))

        wrapped.func_name = fn.func_name
        wrapped.func_doc = '%s -> %s\n\n%s' % (fn.func_name, self.expected, fn.func_doc or '')
        return wrapped


def typecheck(*args, **kws):
    tc = _typecheck(method=True)
    tc.check(*args, **kws)
    return tc


def typecheckfn(*args, **kws):
    tc = _typecheck(method=False)
    tc.check(*args, **kws)
    return tc


def deprecated(fn):
    def wrapped(*args, **kws):
        print 'WARNING: call to deprecated function: %s' % fn.func_name
        return fn(*args, **kws)
    wrapped.func_name = fn.func_name
    wrapped.func_doc  = fn.func_doc
    return wrapped




def checkpicklable(d):
    for v in d.itervalues():
        try:
            slots = v.__slots__
            hasslots = True
        except AttributeError:
            hasslots = False
        try:
            getstate = v.__getstate__
            hasslots = True
        except AttributeError:
            hasgetstate = False

        if hasslots and not hasgetstate:
            print type(v), 'has slots but not __getstate__'

        try:
            d2 = d.__dict__
            checkpicklable(d2)
        except AttributeError: pass
