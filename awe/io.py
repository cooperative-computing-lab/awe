# -*- mode: Python; indent-tabs-mode: nil -*-  #
"""
This file is part of AWE
Copyright (C) 2012- University of Notre Dame
This software is distributed under the GNU General Public License.
See the file COPYING for details.
"""


TRACE = False

class trace(object):
    def __init__(self, values=False):
        self.print_values = values

    def prettyargs(self, *args, **kws):
        if self.print_values:
            pretty = lambda a: repr(a)
        else:
            pretty = lambda a: str(type(a))

        pargs = map(pretty, args)
        pkws  = [ '%s=%s' % (k, pretty(v)) for k, v in kws.items() ]

        return ', '.join(pargs + pkws)

    def __call__(self, fn):
        def wrapped(*args, **kws):
            global TRACE
            if TRACE:
                print 'TRACE calling %s(%s)' % (fn.func_name, self.prettyargs(*args, **kws))
            return fn(*args, **kws)
        return wrapped



class StringStream(object):

    """
    write strings to a representation in memory
    """

    def __init__(self, s=None):
        """
        Create the stream.
        Parameters:
          s: if *s* is a string, split the string by '\n', if a list, then use it
        """

        self.reset()

        if type(s) is str:
            self._buffer = s.split('\n')
        elif type(s) is list:
            self._buffer = s

    def write(self, s):
        self._read = False
        self._buffer.append(s)

    def reset(self):
        self._buffer = list()
        self._str    = str()
        self._read   = False

    def read(self):
        if not self._read:
            self._str  = ''.join(self._buffer)
            self._read = True

        return self._str

    def readlines(self):
        return self._buffer


def log(string):
    print string
