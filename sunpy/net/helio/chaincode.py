import sunpy.net.helio
from sunpy.util.decorators import deprecated

__all__ = ['Chaincode']


@deprecated("5.0", alternative="sunpy.net.helio.Chaincode")
class Chaincode(sunpy.net.helio.Chaincode):
    pass
