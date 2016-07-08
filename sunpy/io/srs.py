"""
This module implements SRS File Reader.
"""
__author__ = "Sudarshan Konge"
__email__ = "sudk1896@gmail.com"

from astropy.table import Table, Column, vstack
import collections

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

    table = list() #This table holds the three tables which we would later merge.
    indices = list()
    
    for i in range(0, len(lines)):
        if (lines[i][0][0] == 'I'):
            indices.append(i)
    indices.append(len(lines))
    indices.sort()
    
    for i in range(0, len(indices) - 1):
        cols = lines[indices[i]+1]
        temp_table = Table(names = cols, dtype=['object_']*len(cols))
        #If a table is empty, we are not adding it.
        for j in range(indices[i]+2, indices[i+1]):
            temp_string = lines[j]
            temp_array = temp_string
            if (len(temp_array) == len(cols)):
                temp_table.add_row(temp_array)
##            else:
##                temp_table.add_row(['None']*len(cols))
        #Make the table data type aware while
        #you're building it.
        if len(temp_table)>0:
            for cols in temp_table.columns.values():
                try:
                    column = temp_table[cols.name].astype(float)
                    temp_table.replace_column(cols.name, column)
                except ValueError:
                    pass
        table.append(temp_table)

    #"table" now has three different tables i.e.
    #I, IA and II.
    #Make the table columns unit-aware. First convert string to
    #floats, ones which can be converted that is.
    attributes = list() 
    
    for item in table:
        for cols in item.columns.values():
            attributes.append(cols.name)
    attributes = list(set(attributes))
    #"attributes" is for the master table.
    
    #We are adding those columns in the tables
    #that the tables don't have and initializing
    #them with 'None'.
    for item in table:
        for attrs in attributes:
            item_attr = [cols.name for cols in item.columns.values()]
            if attrs not in item_attr:
                new_column = Column(['-']*len(item), name=attrs, dtype='object_')
                item.add_column(new_column)

    #Just add a column for ID
    Map = {0:'I', 1:'IA', 2:'II'}
    #Map is for - > 0th table is I, 1st table is IA, 2nd Table is II.
    for i in range(0, len(table)):
        table[i].add_column(Column(data=[Map[i]]*len(table[i]), name='ID', dtype='object_'))
    
    attributes.insert(0, 'ID')
    master = Table(names=attributes, dtype=['object_']*len(attributes))
    #We will store all the three tables as a single table, basically
    #store all rows of all the three (or less) tables in 'master'

    #Why are we doing This ?
    #We first decide an order of the columns in the master table.
    #This order is arbitrary (choose and fix on any order you like).
    #The columns in the three (or less) tables in 'table', don't follow
    #the same order we fixed on 'master'. We need to make them. Once we
    #do that all that remains is to add all the rows to 'master'
    for items in table:
        #Take care of order of columns.
        dict_of_columns = collections.OrderedDict()
        for columns in items.columns.values():
            dict_of_columns[columns.name] = items[columns.name]
        new_table = Table()
        for cols in attributes:
            new_table.add_column(dict_of_columns[cols])
        for rows in new_table:
            master.add_row(rows)
               
    return master
