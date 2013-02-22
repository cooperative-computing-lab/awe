"""
This file is part of AWE
Copyright (C) 2012- University of Notre Dame
This software is distributed under the GNU General Public License.
See the file COPYING for details.
"""


from distutils.core import setup

import os.path


def get_executables():
    import glob
    return glob.glob(os.path.join('executables', '*'))



setup( author       = "Badi' Abdul-Wahid",
       author_email = 'abdulwahidc@gmail.com',
       url          = 'https://bitbucket.org/badi/awe',
       name         = 'awe',
       packages     = ['awe', 'trax'],
       scripts      = get_executables(),
       platforms    = ['Linux', 'Mac OS X'],
       description  = 'Adaptive Weighted Ensemble method using WorkQueue'
       )
