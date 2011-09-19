#!/usr/bin/env python
#-*- coding:utf-8 -*-
from __future__ import absolute_import
"""
Benchmark for comparing different methods of creating arrays in Python
"""
__authors__ = ["Keith Hughitt, Steven Christe, Albert Shih"]
__email__ = "keith.hughitt@nasa.gov"

import numpy as np
import sys
import benchmark

def main():
    """Main application"""
    timer = benchmark.BenchmarkTimer()
    
    options = timer.parse_arguments()
    timer.print_header("ARRAY CREATION")
    
    run_tests(timer, options.scale_factor)
    
    timer.print_summary()

def run_tests(timer, scale_factor):
    '''Go through each test and print out the results'''
    #initialize timer
    timer.reset()   
    
    #Test 1 - Distance matrix creation (method: list comprehension)
    # 
    #Creates an square matrix where each value in the matrix represents the
    #distance from that point to the center of the matrix.
    size = 4096 * scale_factor
    a = np.array([[x for x in range(size)] for y in range(size)], dtype='f')
    offset = (size - 1) / 2.
    x = (a - offset).flatten()
    y = (a.T - offset).flatten()
    result = np.array(zip(x, y)).reshape(size, size, 2)
    timer.log("%d x %d distance matrix creation (list comp)" % (size, size))
    
    #Test 2 - Distance matrix creation (method: numpy methods)
    #Note: format of result differs from above. effect on performance?
    tempa = (np.arange(size ** 2) %  size) - (size - 1) / 2.
    tempb = tempa.reshape(size, size).transpose().reshape(size ** 2)
    result = np.array(zip(tempa, tempb)) * 1.
    timer.log("%d x %d distance matrix creation (numpy)" % (size, size))
    
if __name__ == '__main__':
    sys.exit(main())