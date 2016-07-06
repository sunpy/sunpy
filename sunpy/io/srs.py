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
        string = ""
        for i in range(0, min(len(arr), 8)):
            string += (arr[i] + ' ')
        lines.append(string)

    #Problem: An array of strings. We need to extract
    #three tables from these I, IA and II.

    table = list()
    indices = list()
    n = len(lines)
    for i in range(0, n):
        if (lines[i][0] == 'I'):
            indices.append(i)
    indices.append(len(lines))
    indices.sort()

    for i in range(0, len(indices) - 1):
        string = lines[indices[i]+1]
        columns = string.split()
        
        
    
