
import awe


class Color(object):
    def __init__(self, name):
        self._name = name
    def __str__(self):
        return self._name
    def __repr__(self):
        return 'Color(%s)' % self._name



class Walker(object):

    def __init__(self, coords=None, weight=0., color=None, cell=None, wid=-1):

        assert type(pdb)    is np.array
        assert type(weight) is float
        assert type(color)  is Color
        assert type(cell)   is int
        assert type(wid)    is int

        self.coords = coords
        self.weight = weight
        self.color  = color
        self.cell   = cell
        self.id     = wid

    natoms = property(lambda self: len(self._coords))


class WalkerGroup(object):

    def __init__(self, count=None, natoms=None, topology=None, walkers=None):

        assert type(count) is int
        assert type(natoms) is int
        assert type(topology) is mdtools.prody.AtomGroup
        # TODO: assert type(walkers) is

        self._count    = count
        self._natoms   = natoms

        self.topology  = topology
        self.positions = -1 * np.ones((count, natoms, 3))
        self.weights   =      np.ones(count)
        self.colors    =      np.empty(count, dtype=str)
        self.cells     =      np.zeros(count, dtype=int)

        self._ix       = 0

        if walkers:
            map(self.add, walkers)


    def topology(self, pdb):

        assert type(pdb) is mdtools.prody.AtomGroup
        self.topology = pdb

    def add(self, walker, ix=None):

        assert type(walker) is Walker

        if ix is None: i = self._ix
        else:          i = ix

        self[i] = walker

        if ix is None:
            self._ix          += 1

    def __setitem__(self, i, walker):

        assert walker.natoms == self._natoms
        assert i              < self._count

        self.positions[ix] = walker.coords
        self.weights  [ix] = walker.weight
        self.colors   [ix] = str(walker.color)
        self.cells    [ix] = walker.cell


    def __getitem__(self, k):
        w = Walker(coords = self.positions    [k],
                   weight = self.weights      [k],
                   color  = Color(self.colors [k]),
                   cell   = self.cells        [k])
        return w

    def get_pdb(self, k):
        pdb = self.topology.copy()
        pdb.setCoords(self.positions[k])
        return pdb

    def get_task_params(self, k):

        ss  = awe.io.StringStream()
        pdb = mdtools.prody.writePDBStream(ss, self.get_pdb(k))

        return {'pdb'    : pdb                    ,
                'weight' : str( self.weights [k]) ,
                'color'  :      self.colors  [k]  ,
                'cell'   : str( self.cells   [k]) }
