# -*- mode: Python; indent-tabs-mode: nil -*-  #
"""
This file is part of AWE
Copyright (C) 2012- University of Notre Dame
This software is distributed under the GNU General Public License.
See the file COPYING for details.
"""


import traceback
import os

class TypeException (Exception): pass

class _typecheck(object):

    """
    awe.util._typecheck

    Utility for checking that the types of parameters passed to functions are
    of the correct type.

    Fields:
        method - the function/method to check
        args   - arguments to the function/method
        kws    - keyword arguments to the function/method
    """

    def __init__(self, method=True):

        """
        awe.util._typecheck.__init__

        Initialize a new _typecheck instance.

        Parameters:
            method - the function/method to verify parameters for 
        
        Returns:
            None 
        """

        self.method = method

    def check(self, *args, **kws):

        """
        awe.util._typecheck.__init__

        Assign the arguments and keyword arguments to their instance fields,
        including accounting for the fact that Python sends a reference to the
        object instance when calling a method.

        Parameters:
            args - arguments to the function/method to check
            kws  - keyword arguments to the function/method to check

        Returns:
            The calling _typecheck instance
        """

        if self.method:
            # If calling a method, prepend None to account for the fact that
            # an object reference is passed along with the argument list.
            self.args = [None] + list(args)
        else:
            self.args = args
        self.kws  = kws

        return self

    def typecheck(self, value, expected, name='', arg=-1):

        """
        awe.util._typecheck.typecheck

        Check that a value's type matches the expected type for that value.

        Parameters: 
            value    - the value to check
            expected - the expected type of the value
            name     - the name of the variable/parameter containing the value
            arg      - the position of the parameter in an argument list
        """

        typ = type(value)
        if typ is not expected:
             raise TypeException, '%s expected: %s, but got: %s' % (name or 'param %s' % arg, expected, typ)

    def __call__(self, fn):
        
        """
        awe.util._typecheck.__call__

        The action taken when a _typecheck instance is called. This allows it
        to wrap functions and have them called like normal while performing the
        _typecheck functionality behind-the-scenes.

        Parameters:
            fn - the function to check and run

        Returns:
            The function wrapped in a function to do typechecking
        """

        def wrapped(*args, **kws):

            try:    # Iterate over the argument list and check type
                i = -1
                for v, t in zip(args, self.args):
                    # if it is a method, skip this part
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

        # Ensure that wrapping the function does not overwrite its
        # visible identity (i.e. it should display the same name and docstring)
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


def abspath(p): return os.path.abspath(os.path.expanduser(os.path.expandvars(p)))

def makedirs_parent(p):
    path = abspath(p)
    d    = os.path.dirname(path)
    if not os.path.exists(d):
        os.makedirs(d)
