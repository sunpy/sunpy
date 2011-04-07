#!/usr/bin/env python
#-*- coding:utf-8 -*-

import sys
import benchmark
import pycuda.autoinit
import pycuda.driver as cuda
import pycuda.gpuarray as gpuarray
import pycuda.curandom as curandom
import scikits.cuda.linalg

def main(argv):
    """Main application"""
    timer = benchmark.BenchmarkTimer()
    
    options = timer.parse_arguments()
    timer.print_header("TIME_TEST3_CUDA")
    
    run_tests(timer, options.scale_factor)
    
    timer.print_summary()

def run_tests(timer, scale_factor):
    """PyCUDA port of time_test3.pro"""
    nofileio = True

    #initialize time
    timer.reset()   

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
    siz = int(384 * math.sqrt(scale_factor))

    # Hmm... scikits.cuda.linalg.transpose doesn't currently support int32
    # May need to find another way to do this
    # a = curandom.rand((siz,siz), dtype=np.int32)
    a = curandom.rand((siz,siz))

    for i in xrange(100):
        b = scikits.cuda.linalg.transpose(a, pycuda.autoinit.device)

    time_test_timer('Transpose %d^2 byte, TRANSPOSE function x 100')
    
if __name__ == '__main__':
    sys.exit(main(sys.argv))
