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

class BenchmarkSuite(benchmark.BenchmarkSuite):
    '''Go through each test and print out the results'''
    #nofileio = True
    def bench_1(self, watch):
        with watch.paused():
            siz = 500
            a = np.arange(siz**2, dtype=np.float32).reshape(siz, siz)
            b = np.arange(siz**2, dtype=np.float32).reshape(siz, siz)

        #Test 1 - Matrix Multiplication Large Arrays (500,500) 10*self.scale times
        nrep = 10 * self.scale
        for i in range(nrep):
            c = np.dot(a, b) #pylint: disable=W0612
        return "Matrix Multiplication Large Arrays (500,500) %d times" % nrep
    def bench_2(self, watch):
        #Test 2 - Matrix Multiplication Small Array (50,50) 10000*self.scale times
        with watch.paused():
            siz = 50
            a = np.arange(siz**2, dtype=np.float32).reshape(siz, siz)
            b = np.arange(siz**2, dtype=np.float32).reshape(siz, siz)
            nrep = 10000 * self.scale
        
        for i in range(nrep):
            c = np.dot(a,b) 
        return "Matrix Multiplication Small Array (50,50) %d times" % nrep
    def bench_3(self, watch):
        #Test 3 - Sorting 1 million elements 10*self.scale times
        with watch.paused():
            nrep = 10 * self.scale
            a = [np.random.uniform(0, 1, 1e6).astype(np.float32) for i in range(nrep)]
    
        for i in range(nrep):
            a[i].sort()
        return "Sorting 1 million elements %d times" % nrep
    def bench_4(self, watch):
        #Test 4 - Moving 1 million elements 1000*self.scale times
        with watch.paused():
            a = np.random.uniform(0, 1, 1e6).astype(np.float32)
            c = np.zeros(1e6)
            nrep = 1000 * self.scale
        
        for i in range(nrep):
            c = a    
        return "Moving 1 million elements %d times" % nrep
    def bench_5(self, watch):
        # Test 5 - indirect addressing 1 million elements 100*self.scale times
        # this test does not seem right, should really create a random index
        # which accesses the elements
        with watch.paused():
            a = 1e6 * np.random.uniform(0, 1, 1e6).astype(np.float32)
            b = np.zeros(1e6).astype(np.uint32)
            nrep = 100 * self.scale
        
        for i in range(nrep):
            c = a[b]
        return 'Indirect addressing 1 million elements %d times' % nrep
    def bench_6(self, watch):
        #Test 6 - Shift 512 by 512 byte and store
        with watch.paused():
            a = np.random.randint(2**32, size=1e6)
            c = np.zeros(1e6, dtype=np.uint32)
            nrep = 1000 * self.scale
        
        for i in range(nrep):
            c = np.roll(a, 12, axis=0)
        return 'Shifting 1 million elements %d times' % nrep
    def bench_7(self, watch):
        #Test 7 - Cosine 1 million elements 100*self.scale times
        with watch.paused():
            nrep = 100 * self.scale
            a = 1e6 * np.random.uniform(0, 1, 1e6).astype(np.float32)
            c = np.zeros(1e6, dtype=np.float32)
        
        for i in range(nrep * 10):
            c = np.cos(a)
        return 'Cosine 1 million elements %d times' % nrep
    def bench_8(self, watch):
        #Test 8 - Alog 1 million elements 100*self.scale
        nrep = 100 * self.scale
        for i in range(nrep):
            np.log(a)
        return ' Alog 1 million elements %d times' % nrep
    def bench_9(self, watch):
        #Test 9 - Alog 1 million elements 100*self.scale
        with watch.paused():
            a = np.arange(1e6).astype(np.uint8)
            b = np.zeros(1e6).astype(np.uint8)
            nrep = 1000 * self.scale
    
        for i in range(nrep):
            np.save('test.npy', a)
            b = np.load('test.npy')
        
        os.remove('test.npy')
        return ' Writing and reading bytarr(1e6) %d times' % nrep

    
if __name__ == '__main__':
    bench = BenchmarkSuite('TIME_TESTR')
    sys.exit(bench.main(sys.argv))
