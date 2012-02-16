
import awe


class IResampler(object):

    """
    Interface that all resampling methods should implement
    """

    def resample(self, walkers):
        """
        Parameters:
          *walkers* : an object of type awe.aweclasses.WalkerGroup

        Returns:
          a new awe.aweclasses.WalkerGroup instance
        """

        raise NotImplementedError

    @awe.trace()
    def __call__(self, walkers):
        assert type(walkers) is awe.aweclasses.WalkerGroup
        ws2 = self.resample(walkers)
        assert type(ws2) is awe.aweclasses.WalkerGroup
        return ws2



class DoNothing(IResampler):

    """
    the identity function
    """

    @awe.trace()
    def resample(self, walkers):
        return walkers
