"""
This file backports the TableAttribute functionality from astropy 4.1 it can be
removed when we depend on astropy 4.1.
"""
from copy import deepcopy

from astropy.table.table import QTable

try:
    from astropy.table.table import TableAttribute

except ImportError:
    class MetaAttribute:
        """
        Descriptor to define custom attribute which gets stored in the object
        ``meta`` dict and can have a defined default.
        This descriptor is intended to provide a convenient way to add attributes
        to a subclass of a complex class such as ``Table`` or ``NDData``.
        This requires that the object has an attribute ``meta`` which is a
        dict-like object.  The value of the MetaAttribute will be stored in a
        new dict meta['__attributes__'] that is created when required.
        Classes that define MetaAttributes are encouraged to support initializing
        the attributes via the class ``__init__``.  For example::
            for attr in list(kwargs):
                descr = getattr(self.__class__, attr, None)
                if isinstance(descr, MetaAttribute):
                    setattr(self, attr, kwargs.pop(attr))
        The name of a ``MetaAttribute`` cannot be the same as any of the following:
        - Keyword argument in the owner class ``__init__``
        - Method or attribute of the "parent class", where the parent class is
        taken to be ``owner.__mro__[1]``.
        :param default: default value
        """
        def __init__(self, default=None):
            self.default = default

        def __get__(self, instance, owner):
            # When called without an instance, return self to allow access
            # to descriptor attributes.
            if instance is None:
                return self

            # Get the __attributes__ dict and create if not there already.
            attributes = instance.meta.setdefault('__attributes__', {})
            try:
                value = attributes[self.name]
            except KeyError:
                if self.default is not None:
                    attributes[self.name] = deepcopy(self.default)
                # Return either specified default or None
                value = attributes.get(self.name)
            return value

        def __set__(self, instance, value):
            # Get the __attributes__ dict and create if not there already.
            attributes = instance.meta.setdefault('__attributes__', {})
            attributes[self.name] = value

        def __set_name__(self, owner, name):
            import inspect
            params = [param.name for param in inspect.signature(owner).parameters.values()
                      if param.kind not in (inspect.Parameter.VAR_KEYWORD,
                                            inspect.Parameter.VAR_POSITIONAL)]

            # Reject names from existing params or best guess at parent class
            if name in params or hasattr(owner.__mro__[1], name):
                raise ValueError(f'{name} not allowed as {self.__class__.__name__}')

            self.name = name

        def __repr__(self):
            return f'<{self.__class__.__name__} name={self.name} default={self.default}>'

    class TableAttribute(MetaAttribute):
        """
        Descriptor to define a custom attribute for a Table subclass.
        The value of the ``TableAttribute`` will be stored in a dict named
        ``__attributes__`` that is stored in the table ``meta``.  The attribute
        can be accessed and set in the usual way, and it can be provided when
        creating the object.
        Defining an attribute by this mechanism ensures that it will persist if
        the table is sliced or serialized, for example as a pickle or ECSV file.
        See the `~astropy.utils.metadata.MetaAttribute` documentation for additional
        details.
        Parameters
        ----------
        default : object
            Default value for attribute
        """

    class QTable(QTable):
        def __init__(self, *args, **kwargs):
            # Handle custom (subclass) table attributes that are stored in meta.
            # These are defined as class attributes using the TableAttribute
            # descriptor.  Any such attributes get removed from kwargs here and
            # stored for use after the table is otherwise initialized. Any values
            # provided via kwargs will have precedence over existing values from
            # meta (e.g. from data as a Table or meta via kwargs).
            meta_table_attrs = {}
            if kwargs:
                for attr in list(kwargs):
                    descr = getattr(self.__class__, attr, None)
                    if isinstance(descr, TableAttribute):
                        meta_table_attrs[attr] = kwargs.pop(attr)

            super().__init__(*args, **kwargs)


__all__ = ['QTable', 'TableAttribute']
