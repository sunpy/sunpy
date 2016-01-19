========================
SSWIDL/SunPy Cheat Sheet
========================

`SolarSoft (SSWIDL) <http://sohowww.nascom.nasa.gov/solarsoft/>`_ is a  
popular IDL software library for solar data analysis, and in fact, many parts 
of SunPy are inspired by data structures and functions in SSWIDL. Though IDL and Python are very different it sometimes helps to consider how to translate
simple tasks between the two languages. The primary packages which provide much of the functionality for scientific data analysis in Python are NumPy and SciPy. In the following we assume that those packages are available to you and that you are imported Numpy is imported as np with the following import statement::

    import numpy as np
    import scipy as sp

In the following examples, a and b could be arrays. For python the arrays must be numpy arrays which can be created simply through::

    np.array(a)

where a is a python list of numbers.

**Relational Operators**

=========  ========
 IDL       Python  
=========  ========
a EQ b     a == b
a LT b     a < b
a GT b     a > b
a GE b     a >= b
a LE b     a <= b
a NE b     a != b
=========  ========

**Logical Operators**

=========  ========
 IDL       Python  
=========  ========
a and b    a and b
a or b     a or b
=========  ========

**Math Functions**

=========  ========
 IDL       Python  
=========  ========
cos(a)     np.cos(a)
alog(a)    np.log(a)
alog10(a)  np.alog(a)
exp(a)     np.exp(a)
=========  ========

**Math Constants**

=========  ========
 IDL       Python  
=========  ========
!pi        np.pi
exp(1)     np.e
=========  ========

**Arrays Sequences**

============  ========
 IDL          Python  
============  ========
indgen(10)    np.arange(0,10)
findgen(10)   np.arange(0,10,dtype=np.float)
============  ========

**Array Creation**

=============  =========
 IDL           Python  
=============  =========
dblarr(3,5)    np.zeros((3,5))
intarr(3,5)    np.zeros((3,5),dtype=np.int)
dblarr(3,5)+1  np.ones((3,5))
intarr(3,5)+9  np.zeros((3,5),dtype=np.int) + 9
boolarr(10)    np.zeros(10,dtype=bool)
identity(3)    np.identity(3)
=============  =========

Many more examples can be found on this `page <http://mathesaurus.sourceforge.net/idl-numpy.html>`_
