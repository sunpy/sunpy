"""
This module implements SRS File Reader.
"""
__author__ = "Sudarshan Konge"
__email__ = "sudk1896@gmail.com"

from astropy.table import Table, Column

__all__ = ['read']

def read(filepath):
    """
    Method for reading an SRS File.
    """
    File = open(filepath, 'r')
    lines = list()
    for line in File:
        arr = line.split()
        store = list()
        for i in range(0, min(len(arr), 8)):
            store.append(arr[i])
        lines.append(store)

    #Problem: An array of strings. We need to extract
    #three tables from these I, IA and II.

    table = list()
    indices = list()
    
    for i in range(0, len(lines)):
        if (lines[i][0][0] == 'I'):
            indices.append(i)
    indices.append(len(lines))
    indices.sort()

    for i in range(0, len(indices) - 1):
        columns = lines[indices[i]+1]
        temp_table = Table(names = cols, dtype=['object_']*len(cols))
        for j in range(indices[i]+2, indices[i+1]):
            temp_string = lines[j]
            temp_array = temp_string.split()
            if (len(temp_array) == len(cols)):
                temp.add_row(temp_array)
            else:
                temp.add_row(['None']*len(cols))
        table.append(temp_table)

    #"table" now has three different tables i.e.
    #I, IA and II.
    #Make the table columns unit-aware. First convert string to
    #floats, ones which can be converted that is.
    for item in table:
        for cols in item.columns.values():
            try:
                column = item[cols.name].astype(float)
                item.replace_column(cols.name, column)
            except ValueError:
                pass
        
        
    
