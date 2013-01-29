

import cPickle as pickle

def checkpoint():
    with open('checkpoint') as fd:
        obj = pickle.load(fd)
        return obj
