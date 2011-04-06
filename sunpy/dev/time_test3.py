#-*- coding:utf-8 -*-
#
# <License info will go here...>
#
# Written:
# Steven Christe <steven.d.christe@nasa.gov> (5-Mar-2011)
# Keith Hughitt <keith.hughitt@nasa.gov>
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

import time as tick
import numpy as np
import datetime
import platform
import math
import random
from scipy import fftpack
from scipy import linalg
from scipy import ndimage

#from collections import deque

time = 0.0
total_time = 0.0
geom_time = 0.0 
ntest = 0
nofileio = True

def time_test_timer(test_summary=None):
    '''Print out a string with time taken for test'''
    global time, total_time, geom_time, ntest
    
    #Get current time
    t = tick.time()
    ntest = ntest + 1
    tt = t - time
    total_time = total_time + tt
    geom_time = geom_time + math.log(tt)
    
    output = '\t%d\t%f\t%s' % (ntest, tt, test_summary)
    print(output)

    time = tick.time()

def time_test_reset():
    '''Reset the global clock'''
    global time
    time = tick.time()

def time_test3(fact=1):
    '''Go through each test and print out the results'''
    global nofileio
    
    #Print system information   
    print_sysinfo()
        
    #initialize time
    time_test_reset()   
    
    #Test 1 - Empty For loop
    nrep = 2000000 * fact
    for i in xrange(nrep):
        pass
    time_test_timer("Empty For loop %d times." % nrep)

    #Test 2 - Empty procedure    
    time_test_dummy = lambda x: 0
 
    nrep = 1000000 * fact
    for i in xrange(nrep):
        time_test_dummy(1)
    time_test_timer("Call empty procedure (1 param) %d times." % nrep)
    
    #Test 3 - Add 200000 scalar ints
    nrep = 2000000 * fact
    for i in xrange(nrep):
        a = i + 1
    time_test_timer("Add %d integer scalars and store" % nrep)
    
    #Test 4 - Scalar arithmetic loop
    nrep = 50000 * fact
    for i in xrange(nrep):
        a = i + i - 2
        b = a / 2 + 1
        if b != i:
            print "You screwed up", i, a, b
    time_test_timer("%d scalar loops each of 5 ops, 2 =, 1 if" % nrep)

    #Test 5 - Create a 512x512 array filled with bytes(2)
    a = 2 * np.ones([512, 512], dtype=np.uint8)    
    time_test_reset()
    
    #Test 6 - Mult 512 by 512 byte by constant and store
    nrep = 30 * fact
    for i in xrange(nrep):
        b = a * 2
    time_test_timer('Mult 512 by 512 byte by constant and store, %d times.' % nrep)
    
    #Test 7 - Shift 512 by 512 byte and store
    nrep = 300 * fact
    for i in xrange(nrep):
        c = np.roll(np.roll(b, 10, axis=0), 10, axis=1)
    time_test_timer('Shift 512 by 512 byte and store, %d times.' % nrep)
 
    #Test 8 - Add constant to 512x512 byte array
    nrep = 100 * fact
    for i in xrange(nrep):
        b = a + 3
    time_test_timer('Add constant to 512x512 byte array, %d times' % nrep)

    #Test 9 - Add two 512 by 512 byte arrays and store
    nrep = 80 * fact
    for i in xrange(nrep):
        b = a + b
    time_test_timer('Add two 512 by 512 byte arrays and store, %d times' % nrep)
    
    #a = [[random.random(0, 1) for s in xrange(512)] for s in xrange(512)]
    a = np.random.uniform(0, 1, (512, 512)).astype(np.float32)
    
    #using roll is very inefficient for shifting a float array, 
    #may want to use collections instead instead, need to implement this
    #d = deque(a)

    time_test_reset()
    
    #Test 10 - Mult 512 by 512 floating by constant
    nrep = 30 * fact
    for i in xrange(nrep):
        b = a * 2
    time_test_timer('Mult 512 by 512 floating by constant, %d times.' % nrep)

    #Test 11 - Shift 512 x 512 array
    nrep = 60 * fact
    for i in xrange(nrep):
        c = np.roll(np.roll(b, 10, axis=0), 10, axis=1)
    #for i in xrange(nrep): c = d.rotate(
    time_test_timer('Shift 512 x 512 array, %d times' % nrep)
    
    #Test 12 - Add two 512 by 512 floating images
    nrep = 40 * fact
    for i in xrange(nrep):
        b = a + b
    time_test_timer('Add two 512 by 512 floating images, %d times.' % nrep)

    time_test_reset()

    #Test 13 - Generate random numbers
    nrep = 10 * fact  
    for i in xrange(nrep): 
        a = np.random.uniform(0, 1, 100000)
    time_test_timer('Generated %d random numbers' % (nrep * 100000))

    siz = int(math.sqrt(fact) * 192)
    a = np.random.uniform(0, 1, (siz, siz)).astype(np.float32)
    time_test_reset()

    #Test 14 - Invert random matrix
    b = np.linalg.inv(a)
    time_test_timer('Invert a %d^2 random matrix' % siz)
 
    time_test_reset()
    
    #Test 15 - LU Decomposition of random matrix
    linalg.lu(a)
    time_test_timer('LU Decomposition of a %d^2 random matrix' % siz)

    siz = int(384 * math.sqrt(fact))

    #this following line is not quite right, yields an array 
    #filled differently than the IDL version
    #
    # khughitt 2011/04/06: Try now. Before a and b both referenced
    # same underlying object when indexing
    a = np.arange(siz**2, dtype=np.uint8).reshape(siz, siz)
    b = np.arange(siz**2, dtype=np.uint8).reshape(siz, siz)

    time_test_reset()

    #Test 16 - Transpose byte array with FOR loop
    for i in xrange(siz):
        for j in xrange(siz):
            b[j,i] = a[i,j]
    time_test_timer('Transpose %d^2 byte, FOR loops' % siz)
    
    #Test 17 - Transpose byte array, row and column ops
    for j in xrange(10):
        for i in xrange(siz):
            b[:][i] = a[i][:].transpose()
    time_test_timer('Transpose %d^2 byte, row and column ops x 10' % siz)
  
    #Test 18 - Transpose byte array, TRANSPOSE function
    for i in xrange(100):
        b = a.transpose()
    time_test_timer('Transpose %d^2 byte, TRANSPOSE function x 100' % siz)

    siz = 100000 * fact
    a = np.arange(siz, dtype=np.float32) + 1
    b = np.arange(siz, dtype=np.float32) + 1

    time_test_reset()
    
    #Test 18 - Log of numbers, FOR loop
    for i in xrange(siz):
        b[i] = math.log(a[i])
    time_test_timer('Log of %d numbers, FOR loop' % siz)

	#Test 20 - Log of numbers, vector op
    for i in xrange(10):
        b = np.log(a)

    time_test_timer('Log of %d numbers, vector ops 10 times' % siz)

    n = 2**(17 * fact)
    a = np.arange(n, dtype=np.float32)
    time_test_reset()
    
    #Test 21 - Forward and inverse FFT
    b = fftpack.fft(a)
    b = fftpack.ifft(b)
    time_test_timer('%d point forward plus inverse FFT' % n)
 
    nrep = 10 * fact
    a = np.zeros([512, 512], dtype=np.uint8)
    a[200:250, 200:250] = 10

    time_test_reset()
    
    #Test 21 - Smooth 512 by 512 byte array, 5x5 boxcar
    for i in xrange(nrep):
        b = ndimage.filters.median_filter(a, size=(5, 5))
    time_test_timer('Smooth 512 by 512 byte array, 5x5 boxcar, %d times' % nrep)
 
 	#Test 23 - Smooth 512 by 512 floating point array, 5x5 boxcar
    nrep = 5 * fact
    a = np.zeros([512, 512], dtype=np.float32)
    a[200:250, 200:250] = 10.0
    time_test_reset()
    #need to check to see if this is the same as an IDL smooth
    for i in xrange(nrep):
        b = ndimage.filters.median_filter(a, size=(5, 5))
    time_test_timer('Smooth 512 by 512 floating array, 5x5 boxcar, %d times' % nrep)
    
    a = np.arange(512**2, dtype=np.uint8).reshape(512, 512)

    # aa =assoc(1,a)
    time_test_reset()
    nrep = 40 * fact

    #Test 24 - Write and read 512 by 512 byte array
    if (not nofileio):
        # openw, 1, FILEPATH('test.dat', /TMP), 512, $
        fp = open('/tmp/test.dat', 'r+b') 

        initial = 512 * nrep
        for i in xrange(nrep):
            aa[i] = a
        for i in xrange(nrep):
            a = aa[i]
        time_test_timer('Write and read 512 by 512 byte array x ' + str(nrep))
        fp.close()
    else:
        print('\t\t\tSkipped read/write test')
            
    # Print results
    print_summary()             

    # Remove the data file
    if (not nofileio):
        os.remove('/tmp/test.dat')       
        
def print_sysinfo():
    """Prints the output header containing system and time information"""
    header = ("|TIME_TEST3 performance for Python %s (%s)\n"
              "|\tOS_FAMILY=%s, OS=%s, ARCH=%s %s\n"
              "|\t%s") % (
        platform.python_version(), platform.python_build()[0], 
        platform.system(), " ".join(platform.dist()), 
        platform.processor(), platform.machine(), 
        datetime.datetime.today().ctime()
    )

    #Display header information
        print header
        
def print_summary():
    """Prints a summary of the test results"""
    global total_time, geom_time, ntest

    geom_mean = math.exp(geom_time / ntest)
    summary = ("\t%f=Total Time,"
               "\t%f=Geometric mean,"
               "\t%d tests.") % (total_time, geom_mean, ntest)

    print(summary) 

def time_test3_cuda(fact=1):
    """PyCUDA port of time_test3.pro"""
    import pycuda.autoinit
    import pycuda.driver as cuda
    import pycuda.gpuarray as gpuarray
    import pycuda.curandom as curandom
    import scikits.cuda.linalg
    
    #
    # khughitt (2011/04/04): Non-CUDA tests from above will go here...
    #
    # Perhaps each test should have its own function? Overhead added would
    # likely be small compared to time it takes to run the tests themselves.
    #

    # Initialize linear algebra extensions to PyCUDA
    scikits.cuda.linalg.init()
    
    #
    # Begin CUDA tests
    #
    siz = int(384 * math.sqrt(fact))

    # Hmm... scikits.cuda.linalg.transpose doesn't currently support int32
    # May need to find another way to do this
    # a = curandom.rand((siz,siz), dtype=np.int32)
    a = curandom.rand((siz,siz))

    for i in xrange(100):
        b = scikits.cuda.linalg.transpose(a, pycuda.autoinit.device)

    time_test_timer('Transpose %d^2 byte, TRANSPOSE function x 100')
