# -*- mode: Python; indent-tabs-mode: nil -*-  #
"""
This file is part of AWE
Copyright (C) 2012- University of Notre Dame
This software is distributed under the GNU General Public License.
See the file COPYING for details.
"""

from .io_tools import trace, log
from .util import typecheck, typecheckfn

from . import stats
from . import aweclasses
from . import workqueue
from . import io_tools
from . import resample
from . import structures

from .aweclasses import Walker, AWE, Cell, System, SinkStates
from .workqueue import Config
from .stats import time
from .structures import PDB

from . import voronoi

import prody

### compatibility with different versions of prody
try:
    prody.setVerbosity('error')
except AttributeError:
    prody.confProDy(verbosity='error')

stats.time.start()
