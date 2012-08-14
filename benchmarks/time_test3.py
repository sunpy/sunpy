#!/usr/bin/env python
#-*- coding:utf-8 -*-
#
# <License info will go here...>
#
# Written: 5-Mar-2011
#
# Author: Steven Christe <steven.d.christe@nasa.gov>
# Author: Keith Hughitt <keith.hughitt@nasa.gov>
#
# This module is not yet finished. The code in the comments below 
# are what remains to be implemented.
#
# TODO
# 1] Implement File IO test (Test 24)
# 2] Implement smooth test
# 7] Check implementation of fft test
# 8] Check implementation of smooth test
# 9] Need to optimize shifting code for float arrays. Roll is very inefficient!
#
#pylint: disable=C0103,R0912,R0915,W0404
'''Equivalent of time_test3.pro, a performance test, in IDL

    The tests are 
    Test 1 - Empty For Loop
    Test 2 - Empty procedure
    Test 3 - Add 200000 scalar ints
    Test 4 - Scalar arithmetic loop
    Test 5 - Create a 512x512 array filled with bytes(2)
    Test 6 - Mult 512 by 512 byte by constant and store
    Test 7 - Shift 512 by 512 byte and store
    Test 8 - Add constant to 512x512 byte array
    Test 9 - Add two 512 by 512 byte arrays and store
    Test 10 - Mult 512 by 512 floating by constant
    Test 11 - Shift 512 x 512 array
    Test 12 - Add two 512 by 512 floating images
    Test 13 - Generate random numbers
    Test 14 - Invert random matrix
    Test 15 - LU Decomposition of random matrix
    Test 16 - Transpose byte array with FOR loop
    Test 17 - Transpose byte array, row and column ops
    Test 18 - Transpose byte array, TRANSPOSE function
    Test 20 - Log of numbers, vector op
    Test 20 - Forward and inverse FFT
    Test 21 - Smooth 512 by 512 byte array, 5x5 boxcar
     Test 23 - Smooth 512 by 512 floating array, 5x5 boxcar
    Test 24 - Write and read 512 by 512 byte array
'''
import numpy as np
import scipy.fftpack
import scipy.ndimage
import scipy.linalg
import math
import sys
import os
import benchmark

#from collections import deque

def main():
    """Main application"""
    timer = benchmark.BenchmarkTimer()
    
    options = timer.parse_arguments()
    timer.print_header("TIME_TEST3")
    
    run_tests(timer, options.scale_factor)
    
    timer.print_summary()

def run_tests(timer, scale_factor):
    '''Go through each test and print out the results'''
    nofileio = True

    #initialize time
    timer.reset()   
    
    #Test 1 - Empty For loop
    nrep = 2000000 * scale_factor
    for i in xrange(nrep):
        pass
    timer.log("Empty For loop %d times." % nrep)

    #Test 2 - Empty procedure    
    trash = lambda x: 0
 
    nrep = 100000 * scale_factor
    for i in xrange(nrep):
        trash(i)
    timer.log("Call empty procedure (1 param) %d times." % nrep)
    
    #Test 3 - Add 200000 scalar ints
    nrep = 200000 * scale_factor
    for i in xrange(nrep):
        a = i + 1
    timer.log("Add %d integer scalars and store" % nrep)
    
    #Test 4 - Scalar arithmetic loop
    nrep = 50000 * scale_factor
    for i in xrange(nrep):
        a = i + i - 2
        b = a / 2 + 1
        if b != i:
            print(("You screwed up", i, a, b))
    timer.log("%d scalar loops each of 5 ops, 2 =, 1 if" % nrep)

    # Create a 512x512 array filled with bytes
    a = 2 * np.ones([512, 512], dtype=np.uint8)    
    timer.reset()
    
    #Test 5 - Mult 512 by 512 byte by constant and store
    nrep = 30 * scale_factor
    for i in xrange(nrep):
        b = a * 2
    timer.log('Mult 512 by 512 byte by constant and store, %d times.' % nrep)
    
    #Test 6 - Shift 512 by 512 byte and store
    nrep = 300 * scale_factor
    for i in xrange(nrep):
        c = np.roll(np.roll(b, 10, axis=0), 10, axis=1) #pylint: disable=W0612
    timer.log('Shift 512 by 512 byte and store, %d times.' % nrep)
 
    #Test 7 - Add constant to 512x512 byte array
    nrep = 100 * scale_factor
    for i in xrange(nrep):
        b = a + 3
    timer.log('Add constant to 512x512 byte array, %d times' % nrep)

    #Test 8 - Add two 512 by 512 byte arrays and store
    nrep = 80 * scale_factor
    for i in xrange(nrep):
        b = a + b
    timer.log('Add two 512 by 512 byte arrays and store, %d times' % nrep)
    
    #a = [[random.random(0, 1) for s in xrange(512)] for s in xrange(512)]
    a = np.random.uniform(0, 1, (512, 512)).astype(np.float32)
    
    #using roll is very inefficient for shifting a float array, 
    #may want to use collections instead instead, need to implement this
    #d = deque(a)

    timer.reset()
    
    #Test 9 - Mult 512 by 512 floating by constant
    nrep = 30 * scale_factor
    for i in xrange(nrep):
        b = a * 2
    timer.log('Mult 512 by 512 floating by constant, %d times.' % nrep)

    #Test 10 - Shift 512 x 512 array
    nrep = 60 * scale_factor
    for i in xrange(nrep):
        c = np.roll(np.roll(b, 10, axis=0), 10, axis=1)
    #for i in xrange(nrep): c = d.rotate(
    timer.log('Shift 512 x 512 array, %d times' % nrep)
    
    #Test 11 - Add two 512 by 512 floating images
    nrep = 40 * scale_factor
    for i in xrange(nrep):
        b = a + b
    timer.log('Add two 512 by 512 floating images, %d times.' % nrep)

    timer.reset()

    #Test 12 - Generate random numbers
    nrep = 10 * scale_factor  
    for i in xrange(nrep): 
        a = np.random.uniform(0, 1, 100000)
    timer.log('Generated %d random numbers' % (nrep * 100000))

    siz = int(math.sqrt(scale_factor) * 192)
    a = np.random.uniform(0, 1, (siz, siz)).astype(np.float32)
    timer.reset()

    #Test 13 - Invert random matrix
    b = np.linalg.inv(a)
    timer.log('Invert a %d^2 random matrix' % siz)
 
    timer.reset()
    
    #Test 14 - LU Decomposition of random matrix
    scipy.linalg.lu(a) #pylint: disable=E1101
    timer.log('LU Decomposition of a %d^2 random matrix' % siz)

    siz = int(384 * math.sqrt(scale_factor))

    #this following line is not quite right, yields an array 
    #filled differently than the IDL version
    #
    # khughitt 2011/04/06: Try now. Before a and b both referenced
    # same underlying object when indexing
    a = np.arange(siz**2, dtype=np.uint8).reshape(siz, siz)
    b = np.arange(siz**2, dtype=np.uint8).reshape(siz, siz)

    timer.reset()

    #Test 15 - Transpose byte array with FOR loop
    for i in xrange(siz):
        for j in xrange(siz):
            b[j, i] = a[i, j]
    timer.log('Transpose %d^2 byte, FOR loops' % siz)
    
    #Test 16 - Transpose byte array, row and column ops
    for j in xrange(10):
        for i in xrange(siz):
            b[:][i] = a[i][:].transpose()
    timer.log('Transpose %d^2 byte, row and column ops x 10' % siz)
  
    #Test 17 - Transpose byte array, TRANSPOSE function
    for i in xrange(100):
        b = a.transpose()
    timer.log('Transpose %d^2 byte, TRANSPOSE function x 100' % siz)

    siz = 100000 * scale_factor
    a = np.arange(siz, dtype=np.float32) + 1
    b = np.arange(siz, dtype=np.float32) + 1

    timer.reset()
    
    #Test 18 - Log of numbers, FOR loop
    for i in xrange(siz):
        b[i] = math.log(a[i])
    timer.log('Log of %d numbers, FOR loop' % siz)

    #Test 19 - Log of numbers, vector op
    for i in xrange(10):
        b = np.log(a)

    timer.log('Log of %d numbers, vector ops 10 times' % siz)

    n = 2**(17 * scale_factor)
    a = np.arange(n, dtype=np.float32)
    timer.reset()
    
    #Test 20 - Forward and inverse FFT
    b = scipy.fftpack.fft(a)
    b = scipy.fftpack.ifft(b)
    timer.log('%d point forward plus inverse FFT' % n)
 
    nrep = 10 * scale_factor
    a = np.zeros([512, 512], dtype=np.uint8)
    a[200:250, 200:250] = 10

    timer.reset()
    
    #Test 21 - Smooth 512 by 512 byte array, 5x5 boxcar
    for i in xrange(nrep):
        b = scipy.ndimage.filters.uniform_filter(a, size=(5, 5))
    timer.log('Smooth 512 by 512 byte array, 5x5 boxcar, %d times' % nrep)
 
    #Test 22 - Smooth 512 by 512 floating point array, 5x5 boxcar
    nrep = 5 * scale_factor
    a = np.zeros([512, 512], dtype=np.float32)
    a[200:250, 200:250] = 10.0
    timer.reset()
    #need to check to see if this is the same as an IDL smooth
    for i in xrange(nrep):
        b = scipy.ndimage.filters.uniform_filter(a, size=(5, 5))
    timer.log('Smooth 512 by 512 floating array, 5x5 boxcar, %d times' % nrep)
    
    a = np.arange(512**2, dtype=np.uint8).reshape(512, 512)

    # aa =assoc(1,a)
    timer.reset()
    nrep = 40 * scale_factor

    #Test 23 - Write and read 512 by 512 byte array
#    if (not nofileio):
#        # openw, 1, FILEPATH('test.dat', /TMP), 512, $
#        fp = open('/tmp/test.dat', 'r+b') 
#
#        initial = 512 * nrep
#        for i in xrange(nrep):
#            aa[i] = a
#        for i in xrange(nrep):
#            a = aa[i]
#        timer.log('Write and read 512 by 512 byte array x ' + str(nrep))
#        fp.close()
#    else:
#        print('\t\t\t\tSkipped read/write test')

    # Remove the data file
    #if (not nofileio):
    #    os.remove('/tmp/test.dat')
        
if __name__ == '__main__':
    sys.exit(main())

