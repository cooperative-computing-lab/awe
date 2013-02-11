#!/usr/bin/env python

"""
This file is part of trax
Copyright (C) 2012- University of Notre Dame
This software is distributed under the GNU General Public License.
See the file COPYING for details.
"""


from setuptools import setup


setup( author       = "Badi' Abdul-Wahid",
       author_email = 'abdulwahidc@gmail.com',
       url          = 'https://github.com/badi/trax',
       name         = 'trax',
       version      = '1.0',
       packages     = ['trax'],
       test_suite   = 'tests',
       platforms    = ['Linux', 'Mac OS X'],
       description  = 'Transactional logging library'
       )
