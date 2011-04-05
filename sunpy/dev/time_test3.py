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
# 1] Implement File IO test (demomode = False means to output to file) (Test 24)
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
import platform
import numpy as np
import math
import numpy.matlib as matrix
from scipy import fftpack
import scipy
import random
import datetime
from scipy import ndimage
from scipy import linalg

#from collections import deque

timer_common = 0.0
time = 0.0
output_file = 0.0
total_time = 0.0
geom_time = 0.0 
ntest = 0
demomode = True
nofileio = True

def time_test_timer(name=None):
    '''Print out a string with time taken for test'''
    global timer_common, time, lunno, total_time, geom_time, ntest, demomode, output_file
    #Get current time
    t = tick.time()
    ntest = ntest + 1
    tt = t - time
    total_time = total_time + tt
    geom_time = geom_time + math.log(tt)
    
    if demomode:
        print '\t%d\t%f\t%s' % (ntest, tt, name)
        #print '\t' + str(ntest) + '\t' + str(float(tt)) + '\t' + name
    else:
        output_file.write('\t' + str(ntest) + '\t' + str(float(tt)), '\t', name) 
    time = tick.time()

def time_test_reset(n = None):
    '''Reset the global clock'''
    global time
    time = tick.time()

def time_test3(fact=1):
    '''Go through each test and print out the results'''
    global timer_common, time, lunno, total_time, geom_time, ntest, demomode, output_file, nofileio

    ntest = 0

    currentTime = datetime.datetime.today()
    
    # Header
    header = ("|TIME_TEST3 performance for Python %s (%s)\n"
              "|\tOS_FAMILY=%s, OS=%s, ARCH=%s %s\n"
              "|\t%s") % (
        platform.python_version(), platform.python_build()[0], 
        platform.system(), " ".join(platform.dist()), 
        platform.processor(), platform.architecture()[0], 
        currentTime.ctime()
    )

    #Display header information
    if demomode:
        print header
    else:
        f.write(header)
        
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
     b = a/2 + 1
     if b != i: print "You screwed up", i, a, b

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
    nrep = 100*fact
    for i in xrange(nrep): b = a + 3
    time_test_timer('Add constant to 512x512 byte array, ' + str(nrep) + ' times')

    #Test 9 - Add two 512 by 512 byte arrays and store
    nrep = 80*fact
    for i in xrange(nrep): b = a + b
    time_test_timer('Add two 512 by 512 byte arrays and store, ' + str(nrep) + ' times')
    
    a = [[random.uniform(0,1) for s in xrange(512)] for s in xrange(512)]
    #using roll is very inefficient for shifting a float array, 
    #may want to use collections instead instead, need to implement this
    #d = deque(a)

    time_test_reset()
    
    #Test 10 - Mult 512 by 512 floating by constant
    nrep = 30*fact
    for i in xrange(nrep): b = a * 2
    time_test_timer('Mult 512 by 512 floating by constant, ' + str(nrep) + ' times.')

    #Test 11 - Shift 512 x 512 array
    nrep = 60*fact
    for i in xrange(nrep): c = np.roll(np.roll(b,10,axis=0),10,axis=1)
    #for i in xrange(nrep): c = d.rotate(
    time_test_timer('Shift 512 x 512 array, ' +str(nrep) + ' times')
    
    #Test 12 - Add two 512 by 512 floating images
    nrep = 40*fact
    for i in xrange(nrep): b = a + b
    time_test_timer('Add two 512 by 512 floating images, ' + str(nrep) + ' times.')

    time_test_reset()

    #Test 13 - Generate random numbers
    nrep = 10*fact  
    for i in xrange(nrep): 
        a = [[random.uniform(0,1) for s in xrange(100000)]]
    time_test_timer('Generate ' + str(nrep*100000) + ' random numbers')

    siz = math.sqrt(fact) * 192
    a = [[random.uniform(0,1) for s in xrange(int(siz))] for s in xrange(int(siz))]

    time_test_reset()

    #Test 14 - Invert random matrix
    a = matrix.matrix(a)
    b = a.I
    time_test_timer('Invert a '+str(int(siz)) +'^2 random matrix')
 
    time_test_reset()
    
    #Test 15 - LU Decomposition of random matrix
    linalg.lu(a)
    time_test_timer('LU Decomposition of a ' + str(siz) +'^2 random matrix')
    siz = int(384 * math.sqrt(fact))
    a = np.arange(siz**2,dtype = np.uint8) 
    #this following line is not quite right, yields an array filled differently than the IDL version
    a = a.reshape((siz,siz))
    b = a
    time_test_reset()
    
    #Test 16 - Transpose byte array with FOR loop
    for i in xrange(siz):
        for j in xrange(siz):
            b[j,i] = a[i,j]
    time_test_timer('Transpose ' + str(siz) +'^2 byte, FOR loops')
    
    #Test 17 - Transpose byte array, row and column ops
    for j in xrange(10):
        for i in xrange(siz):
            b[:][i] = a[i][:].transpose()
    time_test_timer('Transpose '+str(siz)+'^2 byte, row and column ops x 10')
  
    #Test 18 - Transpose byte array, TRANSPOE function
    for i in xrange(100): b = a.transpose()
    time_test_timer('Transpose ' + str(siz)+ '^2 byte, TRANSPOSE function x 100')

    siz = 100000*fact
    a = np.arange(siz) + 1
    c = a
    b = a
    time_test_reset()
    
    #Test 18 - Log of numbers, FOR loop
    for i in xrange(siz): b[i] = math.log(a[i])
    time_test_timer('Log of ' + str(siz) + ' numbers, FOR loop')
	#need to reset a as the elements have been changed, not sure why
    a = np.arange(siz) + 1
	#Test 20 - Log of numbers, vector op
    for i in xrange(10): b = np.log(a)

    time_test_timer('Log of ' + str(siz) + ' numbers, vector ops 10 times')

    n = 2**(17*fact)
    a = np.arange(n)
    time_test_reset()
    
    #Test 21 - Forward and inverse FFT
    b = fftpack.fft(a)
    b = fftpack.ifft(b)
# b = fft(a,1)
# b = fft(b,-1)
    time_test_timer(str(n) + ' point forward plus inverse FFT')
 
    nrep = 10L*fact
    a = np.zeros([512, 512],dtype = np.uint8)
    a[200:250, 200:250] = 10

    time_test_reset()
    
    #Test 21 - Smooth 512 by 512 byte array, 5x5 boxcar
    for i in xrange(nrep): b = ndimage.filters.median_filter(a, size = (5,5))
    time_test_timer('Smooth 512 by 512 byte array, 5x5 boxcar, ' + str(nrep) + ' times')
 
 	#Test 23 - Smooth 512 by 512 floating array, 5x5 boxcar
    nrep = 5*fact
    a = np.zeros([512,512],dtype = np.float32)
    a[200:250,200:250] = 10.0
    time_test_reset()
    #need to check to see if this is the same as an IDL smooth
    for i in xrange(nrep): b = ndimage.filters.median_filter(a, size = (5,5))
    time_test_timer('Smooth 512 by 512 floating array, 5x5 boxcar, ' + str(nrep) + ' times')
    
    a = np.arange(512**2,dtype = np.uint8) 
    a = a.reshape((512, 512))
# aa =assoc(1,a)
    time_test_reset()
    nrep = 40*fact

    #Test 24 - Write and read 512 by 512 byte array
# IF ((NOT demomode) AND (NOT nofileio)) THEN BEGIN
    if (demomode and not nofileio):
        f.open('test.dat', 'b')
#     openw, 1, FILEPATH('test.dat', /TMP), 512, $
        initial = 512*nrep #Must be changed for vax
        for i in xrange(nrep): aa[i] = a
        for i in xrange(nrep): a = aa[i]
#     FOR i=0, nrep-1 DO aa[i] = a
#     FOR i=0, nrep-1 DO a=aa[i]
        time_test_timer('Write and read 512 by 512 byte array x ' + str(nrep))
#     close, 1
    else:
        if (nofileio and not demomode):
            print('                      Skipped read/write test')
        else:
            print('                      Skipped read/write test in demo mode')

    if demomode:
        print('\t' + str(total_time) + '=Total Time, \t' + str(math.exp(geom_time/ntest)) + '=Geometric mean,\t' + str(ntest) + ' tests.')
    else:
        f.write('\t' + str(total_time) + '=Total Time, \t' + str(math.exp(geom_time/ntest)) + '=Geometric mean,\t' + str(ntest) + ' tests.')
# ;  Remove the data file
    if (not demomode and not nofileio):
        os.remove('test.dat')       

    if type(output_file) is not type(1.0): output_file.close()
    
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
    # a = curandom.rand((siz,siz), dtype=numpy.int32)
    a = curandom.rand((siz,siz))

    for i in xrange(100):
        b = scikits.cuda.linalg.transpose(a, pycuda.autoinit.device)

    time_test_timer('Transpose %d^2 byte, TRANSPOSE function x 100')
