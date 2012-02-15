


class IResampler(object):

    """
    Interface that all resampling methods should implement
    """

    def __call__(self, walkers):
        """
        Parameters:
          *walkers* : an object of type awe.aweclasses.WalkerGroup

        Returns:
          a new awe.aweclasses.WalkerGroup instance
        """

        raise NotImplementedError



class DoNothing(IResampler):

    """
    the identity function
    """

    def __call__(self walkers):
        return walkers
