#!/usr/bin/env python
#-*- coding:utf-8 -*-
#
# <License info will go here...>
#
# Author: Keith Hughitt <keith.hughitt@nasa.gov>
# Author: Steven Christe <steven.d.christe@nasa.gov>
#
#pylint: disable=F0401,W0612
import sys
import math
import numpy as np
import benchmark
import pycuda.autoinit
#import pycuda.driver as cuda
import pycuda.gpuarray as gpuarray
import pycuda.curandom as curandom
import scikits.cuda.linalg
import scikits.cuda.fft
"""
Notes:
   2011/04/04
   scikits.cuda.linalg.transpose doesn't currently support int32, may need to 
   find another way to do this
"""
def main():
    """Main application"""
    timer = benchmark.BenchmarkTimer()
    
    options = timer.parse_arguments()
    timer.print_header("TIME_TEST3_CUDA")
    
    run_tests(timer, options.scale_factor)
    
    timer.print_summary()

def run_tests(timer, scale_factor):
    """PyCUDA port of time_test3.pro"""
    #nofileio = True
    
    # Initialize linear algebra extensions to PyCUDA
    scikits.cuda.linalg.init()

    #initialize time
    timer.reset()   

    #
    # khughitt (2011/04/04): Non-CUDA tests from above will go here...
    #
    
    #
    # Begin CUDA tests
    #
    siz = int(384 * math.sqrt(scale_factor))

    # a = curandom.rand((siz,siz), dtype=np.int32)
    a = curandom.rand((siz,siz))

    timer.reset()

    #Test 17 - Transpose byte array, TRANSPOSE function
    for i in range(100):
        b = scikits.cuda.linalg.transpose(a, pycuda.autoinit.device)
    timer.log('Transpose %d^2 byte, TRANSPOSE function x 100' % siz)
    
    n = 2**(17 * scale_factor)
    a  = gpuarray.arange(n, dtype=np.float32)
    timer.reset()
    
    #Test 20 - Forward and inverse FFT
    b = scikits.cuda.fft.fft(a)
    b = scikits.cuda.fft.ifft(b)
    timer.log('%d point forward plus inverse FFT' % n)
    
if __name__ == '__main__':
    sys.exit(main())
