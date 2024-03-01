from . import goes
from . import adapt

__all__ = ['goes', 'adapt']

# Define a custom __dir__ to restrict tab-completion to __all__
def __dir__():
    return __all__