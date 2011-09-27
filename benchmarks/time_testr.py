#!/usr/bin/env python
#-*- coding:utf-8 -*-
#
# <License info will go here...>
#
# Written: 10-Apr-2011
#
# Author: Steven Christe <steven.d.christe@nasa.gov>
# Author: Keith Hughitt <keith.hughitt@nasa.gov>
#
# This module implements Richard Schwartz's design for a time test.
#
# TODO
# 1] In Test 3, sort works in place so the second sort does not need to do
#    any work. Need to give it more arrays to sort.
#pylint: disable=W0404
'''Richard Schwartz's time test (equivalent of time_testr.pro)

    The tests are 
    Test 1 - Matrix Multiplication Large Arrays (500,500) 10*scale_factor times
    Test 2 - Matrix Multiplication Small Array (50,50) 10000*scale_factor times
    Test 3 - Sorting 1 million elements 10*scale_factor times
    Test 4 - Moving 1 million elements 1000*scale_factor times
    Test 5 - indirect addressing 1 million elements 100*scale_factor times
    Test 6 - shifting 1 million elements 1000*scale_factor times
    Test 7 - cosine 1 million elements 100*scale_factor times
    Test 8 - alog 1 million elements 100*scale_factor
    Test 9 - writing and reading bytarr(1e6) 1000*scale_factor times

'''
import numpy as np
import sys
import benchmark
import os
import sys

def main():
    """Main application"""
    timer = benchmark.BenchmarkTimer()
   
    options = timer.parse_arguments()
    timer.print_header("TIME_TESTR(ichard)")
    
    run_tests(timer, options.scale_factor)
    
    timer.print_summary()

def run_tests(timer, scale_factor):
    '''Go through each test and print out the results'''
    #nofileio = True

    siz = 500
    a = np.arange(siz**2, dtype=np.float32).reshape(siz, siz)
    b = np.arange(siz**2, dtype=np.float32).reshape(siz, siz)
    
    #initialize time
    timer.reset()

    #Test 1 - Matrix Multiplication Large Arrays (500,500) 10*scale_factor times
    nrep = 10 * scale_factor
    for i in xrange(nrep):
        c = np.dot(a, b) #pylint: disable=W0612
    timer.log("Matrix Multiplication Large Arrays (500,500) %d times" % nrep)

    #Test 2 - Matrix Multiplication Small Array (50,50) 10000*scale_factor times
    siz = 50
    a = np.arange(siz**2, dtype=np.float32).reshape(siz, siz)
    b = np.arange(siz**2, dtype=np.float32).reshape(siz, siz)
    nrep = 10000 * scale_factor
    
    timer.reset()
    for i in xrange(nrep):
        c = np.dot(a,b) 
    timer.log("Matrix Multiplication Small Array (50,50) %d times" % nrep)
    
    #Test 3 - Sorting 1 million elements 10*scale_factor times
    nrep = 10 * scale_factor
    a = [np.random.uniform(0, 1, 1e6).astype(np.float32) for i in xrange(nrep)]
    timer.reset()

    for i in xrange(nrep):
        a[i].sort()
    timer.log("Sorting 1 million elements %d times" % nrep)
    
    #Test 4 - Moving 1 million elements 1000*scale_factor times
    a = np.random.uniform(0, 1, 1e6).astype(np.float32)
    c = np.zeros(1e6)
    nrep = 1000 * scale_factor
    timer.reset()
    
    for i in xrange(nrep):
        c = a    
    timer.log("Moving 1 million elements %d times" % nrep)

    # Test 5 - indirect addressing 1 million elements 100*scale_factor times
    # this test does not seem right, should really create a random index
    # which accesses the elements
    a = 1e6 * np.random.uniform(0, 1, 1e6).astype(np.float32)
    b = np.zeros(1e6).astype(np.uint32)
    nrep = 100 * scale_factor
    timer.reset()
    
    for i in xrange(nrep):
        c = a[b]
    timer.log('Indirect addressing 1 million elements %d times' % nrep)
    
    #Test 6 - Shift 512 by 512 byte and store
    a = np.random.randint(1, 1e6, size=1e6).astype(np.uint32)
    c = np.zeros(1e6, dtype=np.uint32)
    nrep = 1000 * scale_factor
    timer.reset()
    
    for i in xrange(nrep):
        c = np.roll(a, 12, axis=0)
    timer.log('Shifting 1 million elements %d times' % nrep)
 
    #Test 7 - Cosine 1 million elements 100*scale_factor times
    nrep = 100 * scale_factor
    a = 1e6 * np.random.uniform(0, 1, 1e6).astype(np.float32)
    c = np.zeros(1e6, dtype=np.float32)
    timer.reset()
    
    for i in xrange(nrep * 10):
        c = np.cos(a)
    timer.log('Cosine 1 million elements %d times' % nrep)

    #Test 8 - Alog 1 million elements 100*scale_factor
    nrep = 100 * scale_factor
    for i in xrange(nrep):
        np.log(a)
    timer.log(' Alog 1 million elements %d times' % nrep)

    #Test 9 - Alog 1 million elements 100*scale_factor
    a = np.arange(1e6).astype(np.uint8)
    b = np.zeros(1e6).astype(np.uint8)
    nrep = 1000 * scale_factor
    timer.reset()

    for i in xrange(nrep):
        np.save('test.npy', a)
        b = np.load('test.npy')
    
    timer.log(' Writing and reading bytarr(1e6) %d times' % nrep)
    os.remove('test.npy')

    timer.print_summary()
    
if __name__ == '__main__':
    sys.exit(main())
