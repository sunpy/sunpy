import re
import string
import importlib
from abc import ABC, abstractmethod
from textwrap import dedent
from functools import wraps
from collections.abc import Sequence

import numpy as np

import astropy.units as u
from astropy.table import Column, QTable, Row, Table, TableAttribute

from sunpy.util.util import get_width

__all__ = ['QueryResponseColumn', 'BaseQueryResponse',
           'QueryResponseRow', 'QueryResponseTable', 'BaseClient',
           'convert_row_to_table']


class BaseQueryResponse(Sequence):
    """
    An Abstract Base Class for results returned from BaseClient.

    Notes
    -----
    * A QueryResponse object must be able to be instantiated with only one
      iterable argument. (i.e. the ``__init__`` must only have one required
      argument).
    * The `client` property must be settable.
    * The base class does not prescribe how you store the results from your
      client, only that it must be possible to represent them as an astropy
      table in the ``build_table`` method.
    * ``__getitem__`` **must** return an instance of the type it was called on.
      I.e. it must always return an object of ``type(self)``.

    """

    @abstractmethod
    def build_table(self):
        """
        Return an `astropy.table.Table` representation of the query response.
        """

    @property
    @abstractmethod
    def client(self):
        """
        An instance of `BaseClient` used to generate the results.

        Generally this is used to fetch the results later.

        .. note::

            In general, this doesn't have to be the same instance of
            ``BaseClient``, this is left to the client developer. If there is a
            significant connection overhead in creating an instance of a client
            you might want it to be the same instance as used for the search.
        """

    @client.setter
    @abstractmethod
    def client(self, value):
        pass

    @property
    @abstractmethod
    def blocks(self):
        """
        A `collections.abc.Sequence` object which contains the records
        contained within the Query Response.
        """

    def response_block_properties(self):
        """
        Returns a set of class attributes on all the response blocks.

        Returns
        -------
        s : `set`
            List of strings, containing attribute names in the response blocks.
        """
        return set()

    def __str__(self):
        """Print out human-readable summary of records retrieved"""
        return '\n'.join(self.build_table().pformat(show_dtype=False))

    def __repr__(self):
        """Print out human-readable summary of records retrieved"""
        return object.__repr__(self) + "\n" + str(self)

    def _repr_html_(self):
        return self.build_table()._repr_html_()

    def show(self, *cols):
        """
        Returns response tables with desired columns for the Query.

        Parameters
        ----------
        \\*cols : `tuple`
            Name of columns to be shown.

        Returns
        -------
        `astropy.table.Table`
            A table showing values for specified columns.
        """
        table = self.build_table()
        if len(cols) == 0:
            return table
        tablecols = table.columns
        valid_cols = [col for col in cols if col in tablecols]
        return table[valid_cols]


class QueryResponseRow(Row):
    """
    A row subclass which knows about the client of the parent table.
    """

    def as_table(self):
        """
        Return this Row as a length one Table
        """
        return self.table[self.index:self.index + 1]

    def get(self, key, default=None):
        """
        Extract a value from the row if the key is present otherwise return the value of ``default``
        """
        if key in self.colnames:
            return self[key]
        return default

    @property
    def response_block_map(self):
        """
        A dictionary designed to be used to format a filename.

        This takes all the columns in this Row and lower cases them and
        replaces spaces with underscores. Also removes any characters not
        allowed in Python identifiers.
        """
        def key_clean(key):
            key = re.sub(f'[{re.escape(string.punctuation)}]', '_', key)
            key = key.replace(' ', '_')
            key = ''.join(char for char in key
                          if char.isidentifier() or char.isnumeric())
            return key.lower()

        return {key_clean(key): value for key, value in zip(self.colnames, self)}


class QueryResponseColumn(Column):
    """
    A column subclass which knows about the client of the parent table.
    """

    def as_table(self):
        """
        Return this Row as a length one Table
        """
        return self.parent_table[(self.name,)]


class QueryResponseTable(QTable):
    __doc__ = QTable.__doc__

    Row = QueryResponseRow
    Column = QueryResponseColumn

    client = TableAttribute()
    display_keys = TableAttribute(default=slice(None))
    hide_keys = TableAttribute()

    size_column = None

    def unhide_columns(self):
        """
        Modify this table so that all columns are displayed.
        """
        self.display_keys = slice(None)
        self.hide_keys = None
        return self

    def _reorder_columns(self, first_columns, remove_empty=True):
        """
        Generate a new version of this table with ``first_columns`` at the start.

        Parameters
        ----------
        first_columns : list
            The column names to put at the start of the table.
        remove_empty : bool, optional
            Remove columns where all values are `None`.
            Defaults to ``True``.

        Returns
        -------
        new_table : QueryResponseTable
            A sliced version of this table instance so that the columns are
            reordered.
        """
        all_cols = list(self.colnames)
        first_names = [n for n in first_columns if n in all_cols]
        extra_cols = [col for col in all_cols if col not in first_names]
        all_cols = first_names + extra_cols
        new_table = self[[col for col in all_cols if self[col] is not None]]

        if remove_empty:
            empty_cols = [col.info.name for col in self.itercols()
                          if col.info.dtype.kind == 'O' and all(val is None for val in col)]
            new_table.remove_columns(empty_cols)

        return new_table

    @property
    def _display_table(self):
        """
        Apply the display_keys and hide_keys attributes to the table.

        This removes any keys in hide keys and then slices by any keys in
        display_keys to return the correct table.
        """

        keys = list(self.colnames)
        if self.hide_keys:
            # Index only the keys not in hide keys in order
            [keys.remove(key) for key in self.hide_keys if key in keys]

        if self.display_keys != slice(None):
            keys = [dk for dk in self.display_keys if dk in keys]

        table = self[keys]
        # The slicing operation resets display and hide keys to default, but we
        # have already applied it
        table.unhide_columns()

        return table

    def __str__(self):
        """Print out human-readable summary of records retrieved"""
        return '\n'.join(self._display_table.pformat(show_dtype=False))

    def __repr__(self):
        """Print out human-readable summary of records retrieved"""
        return object.__repr__(self) + "\n" + str(self._display_table)

    def _repr_html_(self):
        return QTable._repr_html_(self._display_table)

    def show(self, *cols):
        """
        Return a table with only ``cols`` present.

        If no ``cols`` are specified, all columns will be shown, including any
        hidden by default.

        This differs slightly from ``QueryResponseTable[cols]`` as it allows
        keys which are not in the table to be requested.
        """
        table = self.copy()
        table.unhide_columns()

        if len(cols) == 0:
            return table

        valid_cols = [col for col in cols if col in table.colnames]
        table = table[valid_cols]

        # The slicing operation resets display and hide keys to default, but we
        # want to bypass it here.
        table.unhide_columns()
        return table

    def path_format_keys(self):
        """
        Returns all the names that can be used to format filenames.

        Each one corresponds to a single column in the table, and the format
        syntax should match the dtype of that column, i.e. for a ``Time``
        object or a ``Quantity``.
        """
        rbp = set(self[0].response_block_map.keys())
        for row in self[1:]:
            rbp.intersection(row.response_block_map.keys())
        return rbp

    def total_size(self):
        """
        Returns the total size of all files in a query.

        Derived classes must set the 'size_column' class attribute to make use
        of this.
        """
        if self.size_column not in self.colnames:
            return np.nan * u.byte
        sizes = self[self.size_column]
        # Strip negative filesizes
        total = np.nansum(sizes[sizes > 0])
        if not (total > 0 * u.byte):
            return np.nan * u.byte
        # Find the first power of 3 below the total filesize
        power = 10**(np.floor(np.log10(total.to_value(u.byte)) // 3) * 3)
        # Create mapping from prefix value to prefix name
        prefix_dict = {p[2]: p[0][0] for p in u.core.si_prefixes}
        prefix_unit = u.Unit(f'{prefix_dict[power]}byte')
        return total.to(prefix_unit).round(3)


BaseQueryResponse.register(QueryResponseTable)


def convert_row_to_table(func):
    """
    A wrapper to convert any `~.QueryResponseRow` objects to `~.QueryResponseTable` objects.
    """
    @wraps(func)
    def wrapper(self, query_results, **kwargs):
        if isinstance(query_results, QueryResponseRow):
            query_results = query_results.as_table()
        return func(self, query_results, **kwargs)

    return wrapper


def _print_client(client, html=False, visible_entries=None):
    """
    Given a BaseClient instance will print out each registered attribute.

    Parameters
    ----------
    client : BaseClient
        The instance class to print for.
    html : bool
        Will return a html table instead.

    Returns
    -------
    `str`
        String with the client.
    """
    width = -1 if html else get_width()
    class_name = f"{client.__module__+'.' or ''}{client.__class__.__name__}"
    attrs = client.register_values()
    lines = []
    t = Table(names=["Attr Type", "Name", "Description"],
              dtype=["U80", "U80", "U80"])
    for client_key in attrs.keys():
        # Work around for * attrs having one length.
        if len(attrs[client_key]) == 1 and attrs[client_key][0] == "*":
            t.add_row((client_key.__name__, "All", "All valid values"))
            continue
        for name, desc in attrs[client_key]:
            t.add_row((client_key.__name__, name, desc))
    lines = [class_name, dedent(client.__doc__.partition("\n\n")[0])]
    if html:
        lines = [f"<p>{line}</p>" for line in lines]
    lines.extend(t.pformat(max_lines=visible_entries, show_dtype=False,
                               max_width=width, align="<", html=html))
    return '\n'.join(lines)


class BaseClient(ABC):
    """
    This defines the Abstract Base Class for each download client.

    The BaseClient has several abstract methods that ensure that any subclass enforces the bare minimum API.
    These are `search`, `fetch` and `_can_handle_query`.
    The last one ensures that each download client can be registered with Fido.

    Most download clients should subclass `~sunpy.net.dataretriever.GenericClient`.
    If the structure of `~sunpy.net.dataretriever.GenericClient`
    is not useful you should use `BaseClient`.
    `~sunpy.net.vso.VSOClient` and `~sunpy.net.jsoc.JSOCClient`
    are examples of download clients that subclass ``BaseClient``.
    """

    _registry = dict()

    def __init_subclass__(cls, *args, **kwargs):
        """
        An __init_subclass__ hook initializes all of the subclasses of a given class.
        So for each subclass, it will call this block of code on import.
        This replicates some metaclass magic without the need to be aware of metaclasses.
        Here we use this to register each subclass in a dict that has the `_can_handle_query` attribute.
        This is then passed into the UnifiedDownloaderFactory so we can register them.
        This means that Fido can use the clients internally.
        """
        super().__init_subclass__(**kwargs)

        # We do not want to register GenericClient since its a dummy client.
        if cls.__name__ in ('GenericClient'):
            return

        cls._registry[cls] = cls._can_handle_query

        if hasattr(cls, "_attrs_module"):
            from sunpy.net import attrs

            name, module = cls._attrs_module()
            module_obj = importlib.import_module(module)

            existing_mod = getattr(attrs, name, None)
            if existing_mod and existing_mod is not module_obj:
                raise NameError(f"{name} has already been registered as an attrs name.")

            setattr(attrs, name, module_obj)

            if name not in attrs.__all__:
                attrs.__all__.append(name)

        # Register client attrs after it has registered its own attrs
        from sunpy.net import attr
        values = cls.register_values()
        # If the client has no support, we won't try to register attrs
        if values:
            attr.Attr.update_values({cls: values})

    def __repr__(self):
        """
        Returns the normal repr plus the pretty client __str__.
        """
        return object.__repr__(self) + "\n" + _print_client(visible_entries=15, client=self)

    def __str__(self):
        """
        This enables the "pretty" printing of BaseClient.
        """
        return _print_client(client=self)

    def _repr_html_(self):
        """
        This enables the "pretty" printing of the BaseClient with html.
        """
        return _print_client(visible_entries=15, client=self, html=True)

    @abstractmethod
    def search(self, *args, **kwargs):
        """
        This enables the user to search for data using the client.

        Must return a subclass of `BaseQueryResponse`.
        """

    @abstractmethod
    def fetch(self, query_results, *, path, downloader, **kwargs):
        """
        This enables the user to fetch the data using the client, after a search.

        Parameters
        ----------
        query_results:
            Results to download.
        path : `str` or `pathlib.Path`, optional
            Path to the download directory
        downloader : `parfive.Downloader`
            The download manager to use.

        Returns
        -------
        `parfive.Results`
            The results object, can be `None` if ``wait`` is `False` and
            ``downloader`` is not None.
        """

    @classmethod
    @abstractmethod
    def _can_handle_query(cls, *query):
        """
        This enables the client to register what kind of searches it can handle, to prevent Fido
        using the incorrect client.
        """

    @property
    def info_url(self):
        """
        This should return a string that is a URL to the data server or
        documentation on the data being served.
        """

    @staticmethod
    def check_attr_types_in_query(query, required_attrs={}, optional_attrs={}):
        """
        Check a query against required and optional attributes.

        Returns `True` if *query* contains all the attrs in *required_attrs*,
        and if *query* contains only attrs in both *required_attrs* and *optional_attrs*.
        """
        query_attrs = {type(x) for x in query}
        all_attrs = required_attrs.union(optional_attrs)
        return required_attrs.issubset(query_attrs) and query_attrs.issubset(all_attrs)

    @classmethod
    def register_values(cls, *query):
        """
        This enables the client to register what kind of Attrs it can use directly.

        Returns
        -------
        `dict`
            A dictionary with key values of Attrs and the values are a tuple of
            ("Attr Type", "Name", "Description").
        """
        return {}
