from sunpy.net import Fido
from sunpy.net import attrs as a

result = Fido.search(a.Time('2012/3/4', '2012/3/6'), a.Instrument.lasco)
print(result)





# # Reproduces the bug where Fido.search crashes entirely even
# # if one of the client has an error using a FakeClient

# from astropy.table import Table
# import numpy as np
# arr = np.arange(15).reshape(5, 3)
# t = Table(arr, names=('a', 'b', 'c'), meta={'keywords': {'key1': 'val1'}})
# print(t, type(t))

# print(
# t.columns,   # Dict of table columns (access by column name, index, or slice)
# t.colnames,  # List of column names
# t.meta ,     # Dict of meta-data
# len(t) ,     # Number of table rows

# )

# print(t['a'][0], t[0]['a'])

# arr = np.array([])
# t = Table(arr, names=())
# print(t, type(t))