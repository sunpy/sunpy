""" Converting Astropy Table Object into Pandas Dataframe """
# This module was developed with funding from
# Google Summer of Code 2015
# author - Ankit Kumar  <ankitkmr.iitk@gmail.com>

from pandas import DataFrame

from astropy.table import Table, Column
from astropy.utils import OrderedDict
from astropy.table.column import MaskedColumn

def _to_pandas(astropyTable):
    """
    Return a :class:`pandas.DataFrame` instance

    Returns
    -------
    dataframe : :class:`pandas.DataFrame`
        A pandas :class:`pandas.DataFrame` instance

    Raises
    ------
    ImportError
        If pandas is not installed
    ValueError
        If the Table contains mixin or multi-dimensional columns
    """
    from pandas import DataFrame

    if astropyTable.has_mixin_columns:
        raise ValueError("Cannot convert a table with mixin columns to a pandas DataFrame")

    if any(getattr(col, 'ndim', 1) > 1 for col in astropyTable.columns.values()):
        raise ValueError("Cannot convert a table with multi-dimensional columns to a pandas DataFrame")

    out = OrderedDict()

    for name, column in astropyTable.columns.items():
        if isinstance(column, MaskedColumn):
            if column.dtype.kind in ['i', 'u']:
                out[name] = column.astype(float).filled(np.nan)
            elif column.dtype.kind in ['f', 'c']:
                out[name] = column.filled(np.nan)
            else:
                out[name] = column.astype(np.object).filled(np.nan)
        else:
            out[name] = column

        if out[name].dtype.byteorder not in ('=', '|'):
            out[name] = out[name].byteswap().newbyteorder()

    return DataFrame(out)

