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
'''Equivalent of time_test3.pro, self.a performance test, in IDL

    The tests are 
    Test 1 - Empty For Loop
    Test 2 - Empty procedure
    Test 3 - Add 200000 scalar ints
    Test 4 - Scalar arithmetic loop
    Test 5 - Create self.a 512x512 array filled with bytes(2)
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

class BenchmarkSuite(benchmark.BenchmarkSuite):
        #Test 1 - Empty For loop
    def bench_1(self, watch):
        nrep = 2000000 * self.scale
        for i in range(nrep):
            pass
        return "Empty For loop %d times." % nrep
    
    def bench_2(self, watch): 
        trash = lambda x: 0
     
        nrep = 1000000 * self.scale
        for i in range(nrep):
            trash(i)
        return "Call empty procedure (1 param) %d times." % nrep
        
    def bench_3(self, watch):
        nrep = 2000000 * self.scale
        for i in range(nrep):
            self.a = i + 1
        return "Add %d integer scalars and store" % nrep
        
    def bench_4(self, watch):
        nrep = 50000 * self.scale
        for i in range(nrep):
            self.a = i + i - 2
            self.b = self.a / 2 + 1
            if self.b != i:
                print(("You screwed up", i, self.a, self.b))
        return "%d scalar loops each of 5 ops, 2 =, 1 if" % nrep
    
    def bench_5(self, watch):
        with watch.paused():
            # Create self.a 512x512 array filled with bytes
            self.a = 2 * np.ones([512, 512], dtype=np.uint8)    
        
        #Test 5 - Mult 512 by 512 byte by constant and store
        nrep = 30 * self.scale
        for i in range(nrep):
            self.b = self.a * 2
        return 'Mult 512 by 512 byte by constant and store, %d times.' % nrep
    
    def bench_6(self, watch):
        #Test 6 - Shift 512 by 512 byte and store
        nrep = 300 * self.scale
        for i in range(nrep):
            c = np.roll(np.roll(self.b, 10, axis=0), 10, axis=1) #pylint: disable=W0612
        return 'Shift 512 by 512 byte and store, %d times.' % nrep
    
    def bench_7(self, watch):
        #Test 7 - Add constant to 512x512 byte array
        nrep = 100 * self.scale
        for i in range(nrep):
            self.b = self.a + 3
        return 'Add constant to 512x512 byte array, %d times' % nrep
    
    def bench_8(self, watch):
        #Test 8 - Add two 512 by 512 byte arrays and store
        nrep = 80 * self.scale
        for i in range(nrep):
            self.b = self.a + self.b
        return 'Add two 512 by 512 byte arrays and store, %d times' % nrep
        
        #self.a = [[random.random(0, 1) for s in range(512)] for s in range(512)]
        self.a = np.random.uniform(0, 1, (512, 512)).astype(np.float32)
        
        #using roll is very inefficient for shifting self.a float array, 
        #may want to use collections instead instead, need to implement this
        #d = deque(self.a)
    
    def bench_9(self, watch):
        #Test 9 - Mult 512 by 512 floating by constant
        nrep = 30 * self.scale
        for i in range(nrep):
            self.b = self.a * 2
        return 'Mult 512 by 512 floating by constant, %d times.' % nrep
    
    def bench_10(self, watch):
        #Test 10 - Shift 512 x 512 array
        nrep = 60 * self.scale
        for i in range(nrep):
            c = np.roll(np.roll(self.b, 10, axis=0), 10, axis=1)
        #for i in range(nrep): c = d.rotate(
        return 'Shift 512 x 512 array, %d times' % nrep
    
    def bench_11(self, watch):
        #Test 11 - Add two 512 by 512 floating images
        nrep = 40 * self.scale
        for i in range(nrep):
            self.b = self.a + self.b
        return 'Add two 512 by 512 floating images, %d times.' % nrep
    
    def bench_12(self, watch):
        #Test 12 - Generate random numbers
        nrep = 10 * self.scale  
        for i in range(nrep): 
            self.a = np.random.uniform(0, 1, 100000)
        return 'Generated %d random numbers' % (nrep * 100000)
    
    def bench_13(self, watch):
        #Test 13 - Invert random matrix
        with watch.paused():
            self.siz = int(math.sqrt(self.scale) * 192)
            self.a = np.random.uniform(0, 1, (self.siz, self.siz)).astype(np.float32)
        self.b = np.linalg.inv(self.a)
        return 'Invert self.a %d^2 random matrix' % self.siz

    def bench_14(self, watch):
        #Test 14 - LU Decomposition of random matrix
        scipy.linalg.lu(self.a) #pylint: disable=E1101
        return 'LU Decomposition of self.a %d^2 random matrix' % self.siz
    
    def bench_15(self, watch):
        with watch.paused():
            self.siz = int(384 * math.sqrt(self.scale))
        
            #this following line is not quite right, yields an array 
            #filled differently than the IDL version
            #
            # khughitt 2011/04/06: Try now. Before self.a and self.b both referenced
            # same underlying object when indexing
            self.a = np.arange(self.siz**2, dtype=np.uint8).reshape(self.siz, self.siz)
            self.b = np.arange(self.siz**2, dtype=np.uint8).reshape(self.siz, self.siz)
    
        #Test 15 - Transpose byte array with FOR loop
        for i in range(self.siz):
            for j in range(self.siz):
                self.b[j, i] = self.a[i, j]
        return 'Transpose %d^2 byte, FOR loops' % self.siz
    def bench_16(self, watch):
        #Test 16 - Transpose byte array, row and column ops
        for j in range(10):
            for i in range(self.siz):
                self.b[:][i] = self.a[i][:].transpose()
        return 'Transpose %d^2 byte, row and column ops x 10' % self.siz
    def bench_17(self, watch):
        #Test 17 - Transpose byte array, TRANSPOSE function
        for i in range(100):
            self.b = self.a.transpose()
        return 'Transpose %d^2 byte, TRANSPOSE function x 100' % self.siz
    def bench_18(self, watch):
        #Test 18 - Log of numbers, FOR loop
        with watch.paused():
            self.siz = 100000 * self.scale
            self.a = np.arange(self.siz, dtype=np.float32) + 1
            self.b = np.arange(self.siz, dtype=np.float32) + 1
        for i in range(self.siz):
            self.b[i] = math.log(self.a[i])
        return 'Log of %d numbers, FOR loop' % self.siz
    def bench_19(self, watch):
        #Test 19 - Log of numbers, vector op
        for i in range(10):
            self.b = np.log(self.a)
    
        return 'Log of %d numbers, vector ops 10 times' % self.siz
    def bench_20(self, watch):
        #Test 20 - Forward and inverse FFT
        with watch.paused():
            n = 2**(17 * self.scale)
            self.a = np.arange(n, dtype=np.float32)
        self.b = scipy.fftpack.fft(self.a)
        self.b = scipy.fftpack.ifft(self.b)
        return '%d point forward plus inverse FFT' % n
    def bench_21(self, watch):
        #Test 21 - Smooth 512 by 512 byte array, 5x5 boxcar
        with watch.paused():
            nrep = 10 * self.scale
            self.a = np.zeros([512, 512], dtype=np.uint8)
            self.a[200:250, 200:250] = 10
        for i in range(nrep):
            self.b = scipy.ndimage.filters.median_filter(self.a, size=(5, 5))
        return 'Smooth 512 by 512 byte array, 5x5 boxcar, %d times' % nrep
    def bench_22(self, watch):     
        #Test 22 - Smooth 512 by 512 floating point array, 5x5 boxcar
        with watch.paused():
            nrep = 5 * self.scale
            self.a = np.zeros([512, 512], dtype=np.float32)
            self.a[200:250, 200:250] = 10.0
        #need to check to see if this is the same as an IDL smooth
        for i in range(nrep):
            self.b = scipy.ndimage.filters.median_filter(self.a, size=(5, 5))
        return 'Smooth 512 by 512 floating array, 5x5 boxcar, %d times' % nrep
        
        self.a = np.arange(512**2, dtype=np.uint8).reshape(512, 512)
    
        # aa =assoc(1,self.a)
        nrep = 40 * self.scale

if __name__ == '__main__':
    bench = BenchmarkSuite('TIME_TEST3')
    sys.exit(bench.main(sys.argv))
