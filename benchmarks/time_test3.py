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
# TODO:
# 1] Implement File IO test (Test 24)
# 2] Implement smooth test
# 7] Check implementation of fft test
# 8] Check implementation of smooth test
# 9] Need to optimize shifting code for float arrays. Roll is very inefficient!
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
import math
import random
import sys
import benchmark
from scipy import fftpack
from scipy import linalg
from scipy import ndimage
from optparse import OptionParser
from optparse import IndentedHelpFormatter

#from collections import deque

def main(argv):
    """Main application"""
    options = parse_arguments()
    run_tests(options.scale_factor)

def run_tests(scale_factor):
    '''Go through each test and print out the results'''
    nofileio = True
    
    timer = benchmark.BenchmarkTimer()
    timer.print_header("TIME_TEST3")

    #initialize time
    timer.reset()   
    
    #Test 1 - Empty For loop
    nrep = 2000000 * scale_factor
    for i in xrange(nrep):
        pass
    timer.log("Empty For loop %d times." % nrep)

    #Test 2 - Empty procedure    
    time_test_dummy = lambda x: 0
 
    nrep = 1000000 * scale_factor
    for i in xrange(nrep):
        time_test_dummy(1)
    timer.log("Call empty procedure (1 param) %d times." % nrep)
    
    #Test 3 - Add 200000 scalar ints
    nrep = 2000000 * scale_factor
    for i in xrange(nrep):
        a = i + 1
    timer.log("Add %d integer scalars and store" % nrep)
    
    #Test 4 - Scalar arithmetic loop
    nrep = 50000 * scale_factor
    for i in xrange(nrep):
        a = i + i - 2
        b = a / 2 + 1
        if b != i:
            print "You screwed up", i, a, b
    timer.log("%d scalar loops each of 5 ops, 2 =, 1 if" % nrep)

    #Test 5 - Create a 512x512 array filled with bytes(2)
    a = 2 * np.ones([512, 512], dtype=np.uint8)    
    timer.reset()
    
    #Test 6 - Mult 512 by 512 byte by constant and store
    nrep = 30 * scale_factor
    for i in xrange(nrep):
        b = a * 2
    timer.log('Mult 512 by 512 byte by constant and store, %d times.' % nrep)
    
    #Test 7 - Shift 512 by 512 byte and store
    nrep = 300 * scale_factor
    for i in xrange(nrep):
        c = np.roll(np.roll(b, 10, axis=0), 10, axis=1)
    timer.log('Shift 512 by 512 byte and store, %d times.' % nrep)
 
    #Test 8 - Add constant to 512x512 byte array
    nrep = 100 * scale_factor
    for i in xrange(nrep):
        b = a + 3
    timer.log('Add constant to 512x512 byte array, %d times' % nrep)

    #Test 9 - Add two 512 by 512 byte arrays and store
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
    
    #Test 10 - Mult 512 by 512 floating by constant
    nrep = 30 * scale_factor
    for i in xrange(nrep):
        b = a * 2
    timer.log('Mult 512 by 512 floating by constant, %d times.' % nrep)

    #Test 11 - Shift 512 x 512 array
    nrep = 60 * scale_factor
    for i in xrange(nrep):
        c = np.roll(np.roll(b, 10, axis=0), 10, axis=1)
    #for i in xrange(nrep): c = d.rotate(
    timer.log('Shift 512 x 512 array, %d times' % nrep)
    
    #Test 12 - Add two 512 by 512 floating images
    nrep = 40 * scale_factor
    for i in xrange(nrep):
        b = a + b
    timer.log('Add two 512 by 512 floating images, %d times.' % nrep)

    timer.reset()

    #Test 13 - Generate random numbers
    nrep = 10 * scale_factor  
    for i in xrange(nrep): 
        a = np.random.uniform(0, 1, 100000)
    timer.log('Generated %d random numbers' % (nrep * 100000))

    siz = int(math.sqrt(scale_factor) * 192)
    a = np.random.uniform(0, 1, (siz, siz)).astype(np.float32)
    timer.reset()

    #Test 14 - Invert random matrix
    b = np.linalg.inv(a)
    timer.log('Invert a %d^2 random matrix' % siz)
 
    timer.reset()
    
    #Test 15 - LU Decomposition of random matrix
    linalg.lu(a)
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

    #Test 16 - Transpose byte array with FOR loop
    for i in xrange(siz):
        for j in xrange(siz):
            b[j,i] = a[i,j]
    timer.log('Transpose %d^2 byte, FOR loops' % siz)
    
    #Test 17 - Transpose byte array, row and column ops
    for j in xrange(10):
        for i in xrange(siz):
            b[:][i] = a[i][:].transpose()
    timer.log('Transpose %d^2 byte, row and column ops x 10' % siz)
  
    #Test 18 - Transpose byte array, TRANSPOSE function
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

	#Test 20 - Log of numbers, vector op
    for i in xrange(10):
        b = np.log(a)

    timer.log('Log of %d numbers, vector ops 10 times' % siz)

    n = 2**(17 * scale_factor)
    a = np.arange(n, dtype=np.float32)
    timer.reset()
    
    #Test 21 - Forward and inverse FFT
    b = fftpack.fft(a)
    b = fftpack.ifft(b)
    timer.log('%d point forward plus inverse FFT' % n)
 
    nrep = 10 * scale_factor
    a = np.zeros([512, 512], dtype=np.uint8)
    a[200:250, 200:250] = 10

    timer.reset()
    
    #Test 21 - Smooth 512 by 512 byte array, 5x5 boxcar
    for i in xrange(nrep):
        b = ndimage.filters.median_filter(a, size=(5, 5))
    timer.log('Smooth 512 by 512 byte array, 5x5 boxcar, %d times' % nrep)
 
 	#Test 23 - Smooth 512 by 512 floating point array, 5x5 boxcar
    nrep = 5 * scale_factor
    a = np.zeros([512, 512], dtype=np.float32)
    a[200:250, 200:250] = 10.0
    timer.reset()
    #need to check to see if this is the same as an IDL smooth
    for i in xrange(nrep):
        b = ndimage.filters.median_filter(a, size=(5, 5))
    timer.log('Smooth 512 by 512 floating array, 5x5 boxcar, %d times' % nrep)
    
    a = np.arange(512**2, dtype=np.uint8).reshape(512, 512)

    # aa =assoc(1,a)
    timer.reset()
    nrep = 40 * scale_factor

    #Test 24 - Write and read 512 by 512 byte array
    if (not nofileio):
        # openw, 1, FILEPATH('test.dat', /TMP), 512, $
        fp = open('/tmp/test.dat', 'r+b') 

        initial = 512 * nrep
        for i in xrange(nrep):
            aa[i] = a
        for i in xrange(nrep):
            a = aa[i]
        timer.log('Write and read 512 by 512 byte array x ' + str(nrep))
        fp.close()
    else:
        print('\t\t\t\tSkipped read/write test')

    # Remove the data file
    if (not nofileio):
        os.remove('/tmp/test.dat')
        
    timer.print_summary()
    
def parse_arguments():
    ''' Gets command-line arguments and handles validation '''
    parser = OptionParser("%prog [options]", formatter=IndentedHelpFormatter(4,80))
    parser.add_option("-s", "--scale-factor", dest="scale_factor", type="int",
                      help="factor to scale tests by", metavar="NUM", default=1)

    options, args = parser.parse_args()
    
    options.scale_factor

    return options
    
if __name__ == '__main__':
    sys.exit(main(sys.argv))

