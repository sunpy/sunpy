#-*- coding:utf-8 -*-
#
# Author: Steven Christe <steven.d.christe@nasa.gov>
# Written: 2011/03/30
#
# <License info will go here...>
#
"""Equivalent of time_test3.pro in IDL """

import math
import numpy
import time as tick

timer_common = 0.0
time = 0.0
lunno = 0.0
total_time = 0.0
geom_time = 0.0 
ntest = 0.0
demomode = True

#
# Notes:
# Steven Christe (30-Mar-2011)
#
# This module is not yet finished. The code in the comments below are what remains to be implemented.
#
def time_test_timer(name=None, mode=True):
    global timer_common, time, lunno, total_time, geom_time, ntest, demomode

    #Get current time
    t = tick.time()
    ntest = ntest + 1
    tt = t - time
    total_time = total_time + tt
    geom_time = geom_time
 
    if mode: 
        print ntest, float(tt), ' ', name
    else:
        print lunno, ntest, float(tt), ' ', name 
    time = tick.time()


def time_test3(fact=1):
    
    #initialize time
    time = tick.time()
    
    # Empty For loop
    nrep = 2000000 * fact
    nrepArray = xrange(nrep)
    for i in xrange(nrep): 
        pass

    time_test_timer("Empty For loop %d times." % nrep)
    
    def time_test_dummy(n=None):
        return 0
     
    nrep = 1000000 * fact
    for i in xrange(nrep): 
        time_test_dummy(1)
     
    time_test_timer("Call empty procedure (1 param) %d times." % nrep)
    
    # Add 200000 scalar ints:...
    nrep = 2000000 * fact
    for i in range(nrep):
        a = i + 1
    
    time_test_timer("Add %d integer scalars and store" % nrep)
    
    # Scalar arithmetic loop:
    # Add 200000 scalar ints:...
    nrep = 50000 * fact
    for i in xrange(1000):
        a = i + i - 2
        b = a / 2 + 1
        if b !== i:
            print "You screwed up", i, a, b
    
    time_test_timer("%d scalar loops each of 5 ops, 2 =, 1 if" % nrep)


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
    
    # Using numpy
    # numpy.random.randn(siz,siz).astype(numpy.int32)
    
    # Hmm... scikits.cuda.linalg.transpose doesn't currently support int32
    # May need to find another way to do this
    # a = curandom.rand((siz,siz), dtype=numpy.int32)
    a = curandom.rand((siz,siz))

    for i in xrange(100):
        b = scikits.cuda.linalg.transpose(a, pycuda.autoinit.device)

    time_test_timer('Transpose %d^2 byte, TRANSPOSE function x 100')
    

# ;    Scalar arithmetic loop:
# nrep = long(fact * 50000)
# for i=1L, nrep do begin
#     a = i + i -2
#     b = a / 2 + 1
#     if b ne i then print,'You screwed up',i,a,b
#     endfor
# time_test_timer, strtrim(nrep,2) + ' scalar loops each of 5 ops, 2 =, 1 if)'
# 
# a=replicate(2b,512,512)
# time_test_reset
# nrep = long(30L*fact)
# for i=1,nrep do b=a*2b
# time_test_timer,'Mult 512 by 512 byte by constant and store, 
# '+strtrim(nrep,2)+' times'
# nrep = long(300L*fact)
# for i=1,nrep do c = shift(b,10,10)
# time_test_timer,'Shift 512 by 512 byte and store, '+strtrim(nrep,2)+' times'
# 
# nrep = long(100L*fact)
# for i=1,nrep do b=a+3b
# time_test_timer,'Add constant to 512x512 byte array, '+strtrim(nrep,2)+' 
# times'
# 
# nrep = long(80L*fact)
# for i=1, nrep do b=a+b
# time_test_timer,'Add two 512 by 512 byte arrays and store, 
# '+strtrim(nrep,2)+' times'
# 
# a = randomu(seed, 512,512)
# time_test_reset
# nrep = long(30L*fact)
# for i=1, nrep do b=a*2b
# time_test_timer,'Mult 512 by 512 floating by constant, 
# '+strtrim(nrep,2)+' times'
# 
# nrep = long(60L*fact)
# for i=1,nrep do c = shift(b,10,10)
# time_test_timer,'Shift 512 x 512 array, '+strtrim(nrep,2)+' times'
# 
# nrep = long(40L*fact)
# for i=1, nrep do b=a+b
# time_test_timer,'Add two 512 by 512 floating images, '+strtrim(nrep,2)+' 
# times'
# 
# time_test_reset
# nrep = long(10L*fact)
# for i=1, nrep do a=randomu(qqq, 100000L)    ;Random number matrix
# time_test_timer, 'Generate '+strtrim(100000L*nrep,2)+' random numbers'
# 
# siz = long(sqrt(fact) * 192)
# a = randomu(seed, siz, siz)
# time_test_reset
# b = invert(a)
# time_test_timer,'Invert a '+strtrim(siz,2)+'^2 random matrix'
# 
# time_test_reset
# ludc, a, index
# time_test_timer, 'LU Decomposition of a '+strtrim(siz,2)+'^2 random matrix'
# 
# siz = long(384 * sqrt(fact))
# a=bindgen(siz,siz) & b=a
# time_test_reset
# for i=0,(siz-1) do for j=0,(siz-1) do b[j,i]=a[i,j]
# time_test_timer,'Transpose '+strtrim(siz,2)+'^2 byte, FOR loops'
# for j=1,10 do for i=0,(siz-1) do begin
#     b[0,i] = transpose(a[i,*])
#     end
# time_test_timer,'Transpose '+strtrim(siz,2)+'^2 byte, row and column ops 
# x 10'
# for i=1,100 do b=transpose(a)
# time_test_timer,'Transpose '+strtrim(siz,2)+'^2 byte, TRANSPOSE function 
# x 100'
# 
# siz = long(100000L*fact)
# a=findgen(siz)+1
# c=a
# b = a
# time_test_reset
# for i=0L, n_elements(a)-1 do b[i] = alog(a[i])
# time_test_timer,'Log of '+strtrim(siz,2)+' numbers, FOR loop'
# for i=1,10 do b = alog(a)
# time_test_timer,'Log of '+strtrim(siz,2)+' numbers, vector ops 10 times'
# 
# n = 2L^long(17*fact)
# a = findgen(n)
# time_test_reset
# b = fft(a,1)
# b = fft(b,-1)
# time_test_timer,strtrim(n,2) + ' point forward plus inverse FFT'
# 
# nrep = long(10L*fact)
# a=bytarr(512,512)
# a[200:250,200:250]=10b
# time_test_reset
# for i=1,nrep do b=smooth(a,5)
# time_test_timer,'Smooth 512 by 512 byte array, 5x5 boxcar, 
# '+strtrim(nrep,2)+' times'
# 
# nrep = long(5L*fact)
# a=float(a)
# time_test_reset
# for i=1,nrep do b=smooth(a,5)
# time_test_timer,'Smooth 512 by 512 floating array, 5x5 boxcar, 
# '+strtrim(nrep,2)+' times'
# 
# a=bindgen(512,512)
# aa =assoc(1,a)
# time_test_reset
# nrep = long(40L*fact)
# 
# 
# IF ((NOT demomode) AND (NOT nofileio)) THEN BEGIN
#     openw, 1, FILEPATH('test.dat', /TMP), 512, $
#         initial = 512L*nrep ;Must be changed for vax
#     FOR i=0, nrep-1 DO aa[i] = a
#     FOR i=0, nrep-1 DO a=aa[i]
#     time_test_timer, 'Write and read 512 by 512 byte array x 
# '+strtrim(nrep, 2)
#     close, 1
# END ELSE BEGIN
#     IF (nofileio) AND (NOT demomode) THEN $
#           PRINT,'                      Skipped read/write test' $
#     ELSE $
#           PRINT,'                      Skipped read/write test in demo 
# mode'
# ENDELSE
# 
# IF (demomode) THEN $
#   print, float(total_time),'=Total Time, ', $
#     exp(geom_time / ntest), '=Geometric mean,',ntest,' tests.' $
# ELSE printf, lunno, float(total_time),'=Total Time, ', $
#     exp(geom_time / ntest), '=Geometric mean,',ntest,' tests.'
# 
# ;  Remove the data file
# IF ((NOT demomode) AND (NOT nofileio)) THEN BEGIN
#     openw, 2, FILEPATH('test.dat', /TMP), /DELETE
#     close, 2
# ENDIF
# if lunno gt 0 then free_lun,lunno
