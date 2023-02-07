import sunpy.database
from sunpy.util.decorators import deprecated

__all__ = ['Database']


@deprecated("5.0", alternative="sunpy.database.Database")
class Database(sunpy.database.Database):
    pass
