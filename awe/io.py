# -*- mode: Python; indent-tabs-mode: nil -*-  #
"""
This file is part of AWE
Copyright (C) 2012- University of Notre Dame
This software is distributed under the GNU General Public License.
See the file COPYING for details.
"""


TRACE = False

class trace(object):

    """
    awe.io.trace

    A utility for printing function arguments and keyword arguments for
    development and debugging purposes.

    Fields:
        print_values - flag for printing argument values or types

    Methods:
        prettyargs - transform function arguments into string representations
                     of their values or types
    """

    def __init__(self, values=False):
        
        """
        awe.io.trace.__init__

        Initialize a new instance of the trace class.

        Parameters:
            values - 

        """
        self.print_values = values

    def prettyargs(self, *args, **kws):
        
        """
        awe.io.trace.prettyargs

        Transform a list of arguments and keyword arguments into a string
        representations (either of their types or values).

        Parameters:
            args - a list of function arguments
            kws  - a dictionary of function keyword arguments

        Returns:
            A comma-separated string listing the args and kws
        """

        if self.print_values:
            # Outputs the argument value
            pretty = lambda a: repr(a)
        else:
            # Outputs the argument type
            pretty = lambda a: str(type(a))

        # Generate the string representations
        pargs = map(pretty, args)
        pkws  = [ '%s=%s' % (k, pretty(v)) for k, v in kws.items() ]

        # Make the arguments a comma-delimited string
        return ', '.join(pargs + pkws)

    def __call__(self, fn):
        
        """
        awe.io.trace.__call__

        Wrap a function to provide trace functionality to its arguments.

        Parameters:
            fn - a function/method to wrap

        Returns:
            The supplied function/method wrapped with trace functionality
        """

        def wrapped(*args, **kws):
            global TRACE
            if TRACE:
                print 'TRACE calling %s(%s)' % (fn.func_name, self.prettyargs(*args, **kws))
            return fn(*args, **kws)
        return wrapped



class StringStream(object):

    """
    awe.io.StringStream

    Store values in a buffer for changing to and from string values. This
    allows an internal representation of arbitrary data as strings in memory.

    Fields:
        _buffer - a list of strings
        _read   - whether the StringStream is in read or write mode
        _str    - a string representation of the buffer

    Methods:
        write     - add an item to the buffer
        reset     - clear the buffer and string representation, and enter
                    write mode
        read      - transform the buffer into a string representation
        readlines - return the buffer
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

        """
        awe.io.StringStream.write

        Write to the buffer (add a string to it).

        Parameters:
            s - a string to add to the buffer

        Returns:
            None
        """

        self._read = False
        self._buffer.append(s)

    def reset(self):

        """
        awe.io.StringStream.clear

        Clear the buffer and string representation, and set to write mode.

        Parameters:
            None

        Returns:
            None
        """

        self._buffer = list()
        self._str    = str()
        self._read   = False

    def read(self):
        
        """
        awe.io.StringStream.read

        Read from the buffer (turn it into a string representation of its
        contents with no separators).
        
        Parameters:
            None

        Returns:
            A string representation of the contents of the buffer.
        """

        if not self._read:
            self._str  = ''.join(self._buffer)
            self._read = True

        return self._str

    def readlines(self):

        """
        awe.io.StringStream.readlines

        Return a list of lines in the string (i.e., the buffer).

        Parameters:
            None

        Returns:
            The buffer (a list)
        """

        return self._buffer


def log(string):
    """
    awe.io.log

    Alias for a print statement. Likely to be used in the event that logging
    requires some other action in the future.

    Parameters:
        string - a value to be printed (any primitive or object with __str__)

    Returns:
        None
    """

    print string
